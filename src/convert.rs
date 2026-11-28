use crate::{
    bed::{
        BedFormat,
        rastair1::{BedRecordsFilterParams, Rastair1BedFormat, VcfKeys, VcfRecord},
        writer::BedWriter,
    },
    io::{
        formats::{FromFileExtension, OutputFormat},
        vcf_writer,
    },
    utils::cli,
};
use clio::ClioPath;
use color_eyre::{
    Section as _,
    eyre::{Result, WrapErr, eyre},
};
use seqair::vcf::{FullRecord, Reader, VcfFile, VcfStream, Writer};
use std::{io::BufWriter, num::NonZeroUsize};
use tracing::{info, warn};

/// Convert between different file formats that rastair supports
///
/// Supported input formats include:
/// - VCF (Variant Call Format)
/// - BCF (Binary Call Format)
///
/// Supported output formats include:
/// - The same as input formats
/// - BED (Browser Extensible Data)
#[derive(Debug, clap::Args)]
pub struct ConvertParams {
    /// Input file
    #[arg(short = 'i', long, default_value = "-")]
    #[arg(help_heading = cli::sections::INPUT, value_hint=clap::ValueHint::FilePath)]
    pub input: ClioPath,

    /// Input file format, guessed from file extension if not specified
    #[arg(short = 'f', long)]
    #[arg(help_heading = cli::sections::INPUT)]
    pub input_format: Option<vcf_writer::VcfFormat>,

    /// Output file
    #[arg(short = 'o', long, default_value = "-")]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub output: ClioPath,

    /// Output file format, guessed from file extension if not specified
    #[arg(short = 'F', long)]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub output_format: Option<OutputFormat>,

    /// BED-specific parameters
    #[command(flatten)]
    pub bed_params: BedRecordsFilterParams,

    /// Total number of threads to use (e.g. for parallel compression)
    #[arg(short='@', long = "threads", env = "RASTAIR_THREADS", default_value_t = std::thread::available_parallelism().map(|n|n.get()).unwrap_or(2).max(1))]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub total_threads: usize,
}

pub fn convert(params: &ConvertParams) -> Result<()> {
    let (input_format, output_format) =
        params.formats().wrap_err("Failed to determine input and output formats")?;

    match (input_format, output_format) {
        // same format: copy the file
        (input, OutputFormat::VcfLike(output)) if input == output => {
            warn!("Input and output formats are the same, no conversion will be performed");

            let mut from = params
                .input
                .clone()
                .open()
                .wrap_err_with(|| format!("Failed to open input `{}`", params.input))?;
            let mut to = params
                .output
                .clone()
                .create()
                .wrap_err_with(|| format!("Failed to create output `{}`", params.output))?;

            std::io::copy(&mut from, &mut to).wrap_err("Failed to copy input to output")?;

            info!("Copied input to output without conversion");
            Ok(())
        }
        // converting from vcf to bcf or vice versa
        (input, OutputFormat::VcfLike(output)) => {
            // The total without the main thread: a third reads and inflates,
            // the rest compresses, at least one thread each.
            let background_threads = params.total_threads.saturating_sub(1);
            let reader_threads = background_threads.div_ceil(3).max(1);
            let writer_threads = background_threads.saturating_sub(reader_threads).max(1);

            let mut reader = params.reader::<FullRecord>(reader_threads)?;

            info!(from=?input, to=?output, "Converting");

            let sink = params
                .output
                .clone()
                .create()
                .wrap_err_with(|| format!("Failed to create output `{}`", params.output))?;
            let mut writer = Writer::new(BufWriter::new(sink), output.into())
                .compression_threads(writer_threads)
                .write_header(reader.header())
                .wrap_err_with(|| {
                    format!("Failed to create VCF writer for output file `{}`", params.output)
                })?;

            while let Some(record) =
                reader.next_record().wrap_err("Failed to read record from VCF file")?
            {
                writer
                    .copy_record(record)
                    .and_then(|copy| copy.emit())
                    .wrap_err("Failed to write record to output file")?;
            }
            let (sink, _index) = writer.finish().wrap_err("Failed to finish the output file")?;
            sink.into_inner()
                .map_err(|error| error.into_error())
                .wrap_err("Failed to flush the output file")?;

            Ok(())
        }
        (_, OutputFormat::Bed(format)) => {
            vcf_to_bed(params, format).wrap_err("Failed to convert VCF to BED")
        }
    }
}

impl ConvertParams {
    fn formats(&self) -> Result<(vcf_writer::VcfFormat, OutputFormat)> {
        let input_format = match self.input_format {
            Some(format) => format,
            None => {
                if self.input.is_std() {
                    return Err(eyre!("Input is stdin but no input format was specified"))
                        .note("Please specify the input format with `--input-format`");
                } else {
                    vcf_writer::VcfFormat::guess_format(&self.input)
                        .wrap_err("Failed to guess input format from file extension")?
                }
            }
        };

        let output_format = match self.output_format {
            Some(format) => format,
            None => {
                if self.output.is_std() {
                    return Err(eyre!("Output is stdout but no output format was specified")
                        .note("Please specify the output format with `--output-format`"));
                } else {
                    OutputFormat::guess_format(&self.output)
                        .wrap_err("Failed to guess output format from file extension")?
                }
            }
        };

        Ok((input_format, output_format))
    }
}

impl ConvertParams {
    /// A reader of the input, a file or stdin, inflating on `threads` threads.
    fn reader<F: seqair::vcf::Fields>(&self, threads: usize) -> Result<Reader<F>> {
        let threads = NonZeroUsize::new(threads).unwrap_or(NonZeroUsize::MIN);
        if self.input.is_std() {
            VcfStream::stdin()
                .reader::<F>()
                .threads(threads)
                .open()
                .wrap_err("Failed to open stdin to read VCF file")
        } else if self.input.is_file() {
            VcfFile::open(self.input.path())
                .and_then(|file| file.reader::<F>().threads(threads).open())
                .wrap_err_with(|| format!("Failed to open VCF file `{}`", self.input))
        } else {
            Err(eyre!("VCF input can only be file or stdin").note("If you need to convert from other source, please open an issue in the rastair2 repository"))
        }
    }
}

fn vcf_to_bed(params: &ConvertParams, format: BedFormat) -> Result<()> {
    let mut reader = params.reader::<VcfRecord>(2)?;
    let keys = VcfKeys::of(&reader.keys()).wrap_err("Failed to read the VCF header")?;

    let mut writer = BedWriter::new(&params.output, format).wrap_err_with(|| {
        format!("Failed to create BED writer for output file `{}`", params.output)
    })?;

    // `call` writes a position's main record first, then under `--all` its
    // rejected alts and its indels at the same position; the BED has one row
    // per position, as `call --bed` writes it.
    let mut last_row = None;
    while let Some(record) = reader.next_record().wrap_err("Failed to read record from VCF file")? {
        let locus = &record.0;
        let at = (locus.contig().index(), locus.start());
        if last_row == Some(at) {
            continue;
        }

        let Some(row) = Rastair1BedFormat::from_vcf(&keys, record, &params.bed_params)
            .wrap_err("Failed to convert record to BED format")?
        else {
            continue;
        };
        writer.write_record(&row).wrap_err("Failed to write record")?;
        last_row = Some(at);
    }

    Ok(())
}
