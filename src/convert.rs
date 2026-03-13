use crate::{
    bed::{
        BedFormat,
        rastair1::{BedRecordsFilterParams, Rastair1BedFormat},
        writer::BedWriter,
    },
    call::{RecordFilters, variant_calling::ErrorModel},
    io::{
        formats::{FromFileExtension, InputFormat, OutputFormat},
        mpk::{MessagePackReader, MpkEntry},
        vcf_writer,
    },
    utils::{cli, logging::ThisIsABug},
};
use clio::ClioPath;
use color_eyre::{
    Section as _,
    eyre::{ContextCompat, Result, WrapErr, bail, eyre},
};
use rastair_types::Probability;
use rust_htslib::bcf::{self, Read as _};
use std::num::NonZeroUsize;
use tracing::{debug, info, warn};

/// Convert between different file formats that rastair supports
///
/// Supported input formats include:
/// - VCF (Variant Call Format)
/// - BCF (Binary Call Format)
/// - Message Pack (rastair's internal format)
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
    pub input_format: Option<InputFormat>,

    /// Output file
    #[arg(short = 'o', long, default_value = "-")]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub output: ClioPath,

    /// Output file format, guessed from file extension if not specified
    #[arg(short = 'F', long)]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub output_format: Option<OutputFormat>,

    #[arg(long, default_value = "novaseq6000", value_parser = ErrorModel::value_parser())]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub error_model: ErrorModel,

    /// VCF filters
    #[command(flatten)]
    pub vcf_filter: RecordFilters,

    /// BED-specific parameters
    #[command(flatten)]
    pub bed_params: BedRecordsFilterParams,

    /// Minimum ML score to consider a position as variant
    ///
    /// This does nothing if the input data does not contain ML scores.
    #[arg(long = "bed-ml", default_value_t = Probability::new_panicky(0.5))]
    #[arg(help_heading = cli::sections::FILTER)]
    pub ml_threshold: Probability,

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
        (InputFormat::VcfLike(input), OutputFormat::VcfLike(output)) if input == output => {
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
        // converting from vcf to bcf or vice versa using htslib directly
        (InputFormat::VcfLike(vcf_writer::Format::Vcf(input)), OutputFormat::VcfLike(output)) => {
            use vcf_writer::*;

            // htslib allows setting the number of background compression
            // threads, so we take the total (without the main thread) and give
            // 1/3 to the reader and 2/3 to the writer, with a minimum of 1
            // thread for each.
            let background_threads = params.total_threads.saturating_sub(1);
            let reader_threads = background_threads.div_ceil(3).max(1);
            let writer_threads = background_threads.saturating_sub(reader_threads).max(1);

            let mut reader = bcf::Reader::from_path(params.input.path())
                .wrap_err_with(|| format!("Failed to open VCF file `{}`", params.input))?;

            reader.set_threads(reader_threads).wrap_err("Failed to set VCF reader threads")?;

            let header = bcf::Header::from_template(reader.header());

            let (format, uncompressed) = match output {
                Format::Vcf(VcfFormat::Bcf) => (bcf::Format::Bcf, false),
                Format::Vcf(VcfFormat::Vcf) => (bcf::Format::Vcf, true),
                Format::Vcf(VcfFormat::VcfCompressed) => (bcf::Format::Vcf, false),
                Format::MessagePack => {
                    bail!("Cannot convert {input:?} to MessagePack format")
                }
            };

            info!(from=?input, to=?output, "Converting using htslib");

            let mut writer =
                bcf::Writer::from_path(params.output.path(), &header, uncompressed, format)
                    .wrap_err_with(|| {
                        format!("Failed to create VCF writer for output file `{}`", params.output)
                    })?;
            writer.set_threads(writer_threads).wrap_err("Failed to set VCF writer threads")?;

            for result in reader.records() {
                let record = result.wrap_err("Failed to read record from VCF file")?;
                writer.write(&record).wrap_err("Failed to write record to output file")?;
            }

            Ok(())
        }
        (InputFormat::VcfLike(vcf_writer::Format::Vcf(_vcf)), OutputFormat::Bed(format)) => {
            vcf_to_bed(params, format).wrap_err("Failed to convert VCF to BED")
        }
        (InputFormat::VcfLike(vcf_writer::Format::MessagePack), OutputFormat::Bed(format)) => {
            mpk_to_bed(params, format).wrap_err("Failed to convert MessagePack to BED")
        }
        (
            InputFormat::VcfLike(vcf_writer::Format::MessagePack),
            OutputFormat::VcfLike(vcf_writer::Format::Vcf(format)),
        ) => mpk_to_vcf(params, format).wrap_err("Failed to convert MessagePack to VCF"),
        _ => {
            bail!("Unsupported conversion from {:?} to {:?}", input_format, output_format);
        }
    }
}

impl ConvertParams {
    fn formats(&self) -> Result<(InputFormat, OutputFormat)> {
        let input_format = match self.input_format {
            Some(format) => format,
            None => {
                if self.input.is_std() {
                    return Err(eyre!("Input is stdin but no input format was specified"))
                        .note("Please specify the input format with `--input-format`");
                } else {
                    InputFormat::guess_format(&self.input)
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

fn vcf_to_bed(params: &ConvertParams, format: BedFormat) -> Result<()> {
    let mut reader = if params.input.is_std() {
        bcf::Reader::from_stdin().wrap_err("Failed to open stdin to read VCF file")?
    } else if params.input.is_file() {
        bcf::Reader::from_path(params.input.path())
            .wrap_err_with(|| format!("Failed to open VCF file `{}`", params.input))?
    } else {
        return Err(eyre!("VCF input can only be file or stdin").note("If you need to convert from other source, please open an issue in the rastair2 repository"));
    };
    reader.set_threads(2).wrap_err("Failed to set VCF reader threads")?;

    let mut writer = BedWriter::new(&params.output, format).wrap_err_with(|| {
        format!("Failed to create BED writer for output file `{}`", params.output)
    })?;

    let mut record = reader.empty_record();
    while let Some(res) = reader.read(&mut record) {
        if let Err(error) = res {
            debug!(%error, "Skipping invalid record in VCF file");
            continue;
        };

        let params = crate::bed::rastair1::BedRecordsConvertParams {
            ml_threshold: params.ml_threshold,
            filters: params.bed_params.clone(),
        };
        let Some(record) = Rastair1BedFormat::from_vcf(&record, &params)
            .wrap_err("Failed to convert record to BED format")?
        else {
            continue;
        };
        writer.write_record(&record).wrap_err("Failed to write record")?;
    }

    Ok(())
}

fn mpk_to_vcf(params: &ConvertParams, format: vcf_writer::VcfFormat) -> Result<()> {
    let r = MessagePackReader::new(&params.input)
        .wrap_err("Failed to create MessagePack reader")
        .and_then(|reader| reader.read().wrap_err("Failed to read file header"))
        .wrap_err_with(|| format!("Failed to read MessagePack from `{}`", params.input))?;
    debug!(header=?r.header, "opened mpk file");
    let Some(meta) = r.vcf_header else {
        bail!("MessagePack file does not contain a VCF header");
    };

    let vcf_params = vcf_writer::VcfParams {
        vcf: Some(params.output.clone()),
        vcf_threads: NonZeroUsize::new(4).expect("valid number"),
        vcf_info_fields: Vec::new(),
        vcf_format_fields: Vec::new(),
        vcf_all_fields: false,
    };

    let (format, compression) = format.into();
    let mut writer = vcf_params
        .vcf_writer(&meta.contigs, &meta.samples, &meta.metadata, format, compression)
        .wrap_err_with(|| {
            format!("Failed to create VCF writer for output file `{}`", params.output)
        })?
        .wrap_err("No writer requested")?;

    for entry in r.entries {
        match entry {
            Ok(MpkEntry::Record(record)) => {
                for vcf_record in record
                    .to_vcf_records(Some(params.ml_threshold), &params.error_model)
                    .wrap_err("Failed to convert record to VCF format")
                    .this_is_a_bug()?
                    .to_vec(&params.vcf_filter)
                {
                    writer.add(vcf_record).wrap_err("Failed to write record")?;
                }
            }
            Ok(x) => {
                warn!(?x, "Skipping unsupported entry type in MessagePack file");
                continue;
            }
            Err(e) => Err(e)?,
        }
    }

    Ok(())
}

fn mpk_to_bed(params: &ConvertParams, format: BedFormat) -> Result<()> {
    let r = MessagePackReader::new(&params.input)
        .wrap_err("Failed to create MessagePack reader")
        .and_then(|reader| reader.read().wrap_err("Failed to read file header"))
        .wrap_err_with(|| format!("Failed to read MessagePack file `{}`", params.input))?;

    let mut writer = BedWriter::new(&params.output, format).wrap_err_with(|| {
        format!("Failed to create BED writer for output file `{}`", params.output)
    })?;

    for entry in r.entries {
        match entry {
            Ok(MpkEntry::Record(record)) => {
                let params = crate::bed::rastair1::BedRecordsConvertParams {
                    ml_threshold: params.ml_threshold,
                    filters: params.bed_params.clone(),
                };
                let Some(record) = Rastair1BedFormat::from_metrics(&record, &params)
                    .wrap_err("Failed to convert record to BED format")?
                else {
                    continue;
                };
                writer.write_record(&record).wrap_err("Failed to write record")?;
            }
            Ok(x) => {
                warn!(?x, "Skipping unsupported entry type in MessagePack file");
                continue;
            }
            Err(e) => Err(e)?,
        }
    }

    Ok(())
}
