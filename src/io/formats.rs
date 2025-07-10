use crate::io::vcf_writer::{self, VcfFormat};
use clio::ClioPath;
use color_eyre::{
    Result, Section as _,
    eyre::{ContextCompat as _, bail, eyre},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
    VcfLike(vcf_writer::Format),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    VcfLike(vcf_writer::Format),
    /// BED format
    Bed,
}

impl clap::ValueEnum for InputFormat {
    fn value_variants<'a>() -> &'a [Self] {
        &[
            InputFormat::VcfLike(vcf_writer::Format::Vcf(VcfFormat::Vcf)),
            InputFormat::VcfLike(vcf_writer::Format::Vcf(VcfFormat::Bcf)),
            InputFormat::VcfLike(vcf_writer::Format::Vcf(VcfFormat::VcfCompressed)),
            InputFormat::VcfLike(vcf_writer::Format::MessagePack),
        ]
    }

    fn to_possible_value(&self) -> Option<clap::builder::PossibleValue> {
        match self {
            InputFormat::VcfLike(format) => format.to_possible_value(),
        }
    }
}

impl clap::ValueEnum for OutputFormat {
    fn value_variants<'a>() -> &'a [Self] {
        &[
            OutputFormat::VcfLike(vcf_writer::Format::Vcf(VcfFormat::Vcf)),
            OutputFormat::VcfLike(vcf_writer::Format::Vcf(VcfFormat::Bcf)),
            OutputFormat::VcfLike(vcf_writer::Format::Vcf(VcfFormat::VcfCompressed)),
            OutputFormat::VcfLike(vcf_writer::Format::MessagePack),
            OutputFormat::Bed,
        ]
    }

    fn to_possible_value(&self) -> Option<clap::builder::PossibleValue> {
        match self {
            OutputFormat::VcfLike(format) => format.to_possible_value(),
            OutputFormat::Bed => Some(clap::builder::PossibleValue::new("bed")),
        }
    }
}

pub trait FromFileExtension: Sized {
    fn from_file_extension(path: &str) -> Option<Self>;

    fn guess_format(path: &ClioPath) -> Result<Self> {
        let Some(filename) = path.path().file_name().and_then(|x| x.to_str()) else {
            bail!("No file name found in path `{path}`");
        };

        Self::from_file_extension(filename).wrap_err_with(|| {
            eyre!("Could not determine format from file extension `{filename}`").suggestion(
                "You can specify the format explicitly with `--input-format` or `--output-format`",
            )
        })
    }
}

impl FromFileExtension for vcf_writer::Format {
    fn from_file_extension(p: &str) -> Option<Self> {
        if p.ends_with(".bcf") {
            Some(vcf_writer::Format::Vcf(VcfFormat::Bcf))
        } else if p.ends_with(".vcf.gz") {
            Some(vcf_writer::Format::Vcf(VcfFormat::VcfCompressed))
        } else if p.ends_with(".vcf") {
            Some(vcf_writer::Format::Vcf(VcfFormat::Vcf))
        } else if p.ends_with("mpk.lz4") {
            Some(vcf_writer::Format::MessagePack)
        } else {
            None
        }
    }
}

impl FromFileExtension for InputFormat {
    fn from_file_extension(p: &str) -> Option<Self> {
        vcf_writer::Format::from_file_extension(p).map(InputFormat::VcfLike)
    }
}

impl FromFileExtension for OutputFormat {
    fn from_file_extension(p: &str) -> Option<Self> {
        vcf_writer::Format::from_file_extension(p)
            .map(OutputFormat::VcfLike)
            .or_else(|| if p.ends_with(".bed") { Some(OutputFormat::Bed) } else { None })
    }
}
