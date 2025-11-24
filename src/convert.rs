use crate::{
    bed::{
        BedFormat,
        rastair1::{BedRecordsFilterParams, Rastair1BedFormat},
        writer::BedWriter,
    },
    io::{
        formats::{FromFileExtension, InputFormat, OutputFormat},
        mpk::{MessagePackReader, MpkEntry},
        vcf_writer,
    },
    utils::{cli, logging::ThisIsABug},
    vcf::{DeNovoCpGCandidate, InCpG},
};
use clio::ClioPath;
use color_eyre::{
    Section as _,
    eyre::{ContextCompat, Result, WrapErr, bail, eyre},
};
use rastair_types::Probability;
use rastair_vcf::VcfField;
use rust_htslib::bcf::Read as _;
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

    /// BED-specific parameters
    #[command(flatten)]
    pub bed_params: BedRecordsFilterParams,

    /// Write tabix index for the BED output file
    #[arg(long)]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub bed_index: bool,

    /// Minimum ML score to consider a position as variant
    ///
    /// This does nothing if the input data does not contain ML scores.
    #[arg(long = "bed-ml", default_value_t = Probability::new_panicky(0.8))]
    #[arg(help_heading = cli::sections::FILTER)]
    pub ml_threshold: Probability,
}

pub fn convert(params: &ConvertParams) -> Result<()> {
    let (input_format, output_format) =
        params.formats().wrap_err("Failed to determine input and output formats")?;

    match (input_format, output_format) {
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
        rust_htslib::bcf::Reader::from_stdin().wrap_err("Failed to open stdin to read VCF file")?
    } else if params.input.is_file() {
        rust_htslib::bcf::Reader::from_path(params.input.path())
            .wrap_err_with(|| format!("Failed to open VCF file `{}`", params.input))?
    } else {
        return Err(eyre!("VCF input can only be file or stdin").note("If you need to convert from other source, please open an issue in the rastair2 repository"));
    };
    reader.set_threads(2).wrap_err("Failed to set VCF reader threads")?;

    let mut writer =
        BedWriter::new(&params.output, format, params.bed_index).wrap_err_with(|| {
            format!("Failed to create BED writer for output file `{}`", params.output)
        })?;

    let mut record = reader.empty_record();
    while let Some(res) = reader.read(&mut record) {
        if let Err(error) = res {
            debug!(%error, "Skipping invalid record in VCF file");
            continue;
        };

        let cpg = record.info(InCpG::ID.as_bytes()).flag().unwrap_or_default();
        let dn_cpg = record.info(DeNovoCpGCandidate::ID.as_bytes()).flag().unwrap_or_default();
        if !cpg && !dn_cpg {
            continue;
        }

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
                    .to_vcf_records(Some(params.ml_threshold))
                    .wrap_err("Failed to convert record to VCF format")
                    .this_is_a_bug()?
                {
                    writer.add(&vcf_record).wrap_err("Failed to write record")?;
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

    let mut writer =
        BedWriter::new(&params.output, format, params.bed_index).wrap_err_with(|| {
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
