use crate::{
    bed::{BedWriter, Rastair1BedFormat},
    io::{
        formats::{FromFileExtension, InputFormat, OutputFormat},
        mpk::{MessagePackReader, MpkEntry},
        vcf_writer,
    },
};
use clap::value_parser;
use clio::ClioPath;
use color_eyre::{
    Section as _,
    eyre::{Result, WrapErr, bail, eyre},
};
use rust_htslib::bcf::Read as _;
use std::num::NonZeroUsize;
use tracing::{debug, info, warn};

/// Convert between different file formats that rastair2 supports
///
/// Supported input formats include:
/// - VCF (Variant Call Format)
/// - BCF (Binary Call Format)
/// - Message Pack (rastair2's internal format)
///
/// Supported output formats include:
/// - The same as input formats
/// - BED (Browser Extensible Data)
#[derive(Debug, clap::Args)]
pub struct ConvertParams {
    /// Input file
    #[arg(long, value_parser=value_parser!(ClioPath).exists().is_file())]
    pub input: ClioPath,

    /// Input file format, guessed from file extension if not specified
    #[arg(long)]
    pub input_format: Option<InputFormat>,

    /// Output file
    #[arg(short = 'o', long)]
    pub output: ClioPath,

    /// Output file format, guessed from file extension if not specified
    #[arg(long)]
    pub output_format: Option<OutputFormat>,
}

pub fn convert(params: &ConvertParams) -> Result<()> {
    let (input_format, output_format) =
        params.formats().wrap_err("Failed to determine input and output formats")?;

    match (input_format, output_format) {
        (InputFormat::VcfLike(input), OutputFormat::VcfLike(output)) if input == output => {
            warn!("Input and output formats are the same, no conversion will be performed");
            std::fs::copy(params.input.path(), params.output.path())
                .wrap_err("Failed to copy input file to output file")?;
            info!("Copied input file to output file without conversion");
            Ok(())
        }
        (InputFormat::VcfLike(vcf_writer::Format::Vcf(_vcf)), OutputFormat::Bed) => {
            vcf_to_bed(params).wrap_err("Failed to convert VCF to BED")
        }
        (InputFormat::VcfLike(vcf_writer::Format::MessagePack), OutputFormat::Bed) => {
            mpk_to_bed(params).wrap_err("Failed to convert MessagePack to BED")
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

fn vcf_to_bed(params: &ConvertParams) -> Result<()> {
    let mut reader = rust_htslib::bcf::Reader::from_path(params.input.path())
        .wrap_err_with(|| format!("Failed to open VCF file `{}`", params.input))?;

    let mut writer = BedWriter::new(&params.output).wrap_err_with(|| {
        format!("Failed to create BED writer for output file `{}`", params.output)
    })?;

    for record in reader.records() {
        let Ok(record) = record else {
            warn!("Skipping invalid record in VCF file");
            continue;
        };
        let record = Rastair1BedFormat::try_from(&record)
            .wrap_err("Failed to convert record to BED format")?;
        writer.write_record(&record).wrap_err("Failed to write record")?;
    }

    Ok(())
}

fn mpk_to_vcf(params: &ConvertParams, format: vcf_writer::VcfFormat) -> Result<()> {
    let r = MessagePackReader::new(&params.input)
        .wrap_err("Failed to create MessagePack reader")
        .and_then(|reader| reader.read().wrap_err("Failed to read file header"))
        .wrap_err_with(|| format!("Failed to read MessagePack file `{}`", params.input))?;
    debug!(header=?r.header, "opened mpk file");
    let Some(meta) = r.vcf_header else {
        bail!("MessagePack file does not contain a VCF header");
    };

    let params = vcf_writer::Params {
        vcf_output: params.output.clone(),
        vcf_threads: NonZeroUsize::new(4).expect("valid number"),
    };

    let (format, compression) = format.into();
    let mut writer = params
        .vcf_writer(&meta.contigs, &meta.samples, &meta.metadata, format, compression)
        .wrap_err_with(|| {
            format!("Failed to create VCF writer for output file `{}`", params.vcf_output)
        })?;

    for entry in r.entries {
        match entry {
            Ok(MpkEntry::Record(record)) => {
                writer.add(&record).wrap_err("Failed to write record")?;
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

fn mpk_to_bed(params: &ConvertParams) -> Result<()> {
    let r = MessagePackReader::new(&params.input)
        .wrap_err("Failed to create MessagePack reader")
        .and_then(|reader| reader.read().wrap_err("Failed to read file header"))
        .wrap_err_with(|| format!("Failed to read MessagePack file `{}`", params.input))?;

    let mut writer = BedWriter::new(&params.output).wrap_err_with(|| {
        format!("Failed to create BED writer for output file `{}`", params.output)
    })?;

    for entry in r.entries {
        match entry {
            Ok(MpkEntry::Record(record)) => {
                let record = Rastair1BedFormat::try_from(record.as_ref())
                    .wrap_err("Failed to convert record to BED format")?;
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
