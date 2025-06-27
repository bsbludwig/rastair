use crate::io::{
    formats::{FromFileExtension, InputFormat, OutputFormat},
    mpk::MessagePackReader,
    vcf_writer,
};
use clap::value_parser;
use clio::ClioPath;
use color_eyre::{
    Section as _,
    eyre::{ContextCompat, Result, WrapErr, bail, eyre},
};
use tracing::{info, warn};

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
    let input_format = match params.input_format {
        Some(format) => format,
        None => {
            if params.input.is_std() {
                return Err(eyre!("Input is stdin but no input format was specified")
                    .note("Please specify the input format with `--input-format`"));
            } else {
                guess_format::<InputFormat>(&params.input)
                    .wrap_err("Failed to guess input format from file extension")?
            }
        }
    };

    let output_format = match params.output_format {
        Some(format) => format,
        None => {
            if params.output.is_std() {
                return Err(eyre!("Output is stdout but no output format was specified")
                    .note("Please specify the output format with `--output-format`"));
            } else {
                guess_format::<OutputFormat>(&params.output)
                    .wrap_err("Failed to guess output format from file extension")?
            }
        }
    };

    match (input_format, output_format) {
        (InputFormat::VcfLike(input), OutputFormat::VcfLike(output)) if input == output => {
            warn!("Input and output formats are the same, no conversion will be performed");
            std::fs::copy(params.input.path(), params.output.path())
                .wrap_err("Failed to copy input file to output file")?;
            info!("Copied input file to output file without conversion");
            Ok(())
        }
        (InputFormat::VcfLike(vcf_writer::Format::Vcf(_vcf)), OutputFormat::Bed) => {
            todo!("Implement VCF->BED conversion")
        }
        (InputFormat::VcfLike(vcf_writer::Format::MessagePack), OutputFormat::Bed) => {
            bail!(
                "Cannot convert MessagePack to BED format directly. Please convert to VCF or BCF first."
            )
        }
        (InputFormat::VcfLike(vcf_writer::Format::MessagePack), OutputFormat::VcfLike(_)) => {
            // FIXME: This is just a placeholder for debugging purposes
            let r = MessagePackReader::new(&params.input)
                .wrap_err("Failed to create MessagePack reader")
                .and_then(|reader| reader.read().wrap_err("Failed to read file header"))
                .wrap_err_with(|| format!("Failed to read MessagePack file `{}`", params.input))?;
            info!(header=?r.header, "opened mpk file");
            dbg!(&r.vcf_header);
            dbg!(r.entries.count());
            Ok(())
        }
        _ => {
            bail!("Unsupported conversion from {:?} to {:?}", input_format, output_format);
        }
    }
}

fn guess_format<T: FromFileExtension>(path: &ClioPath) -> Result<T> {
    let Some(name) = path.path().file_name().and_then(|x| x.to_str()) else {
        bail!("No file name found in path `{path}`");
    };

    T::from_file_extension(name).wrap_err_with(|| {
        eyre!("Could not determine format from file extension `{name}`").suggestion(
            "You can specify the format explicitly with `--input-format` or `--output-format`",
        )
    })
}

// fn mpack_to_vcf(input: &ClioPath, output: &ClioPath) -> Result<()> {
//     todo!()
// }
