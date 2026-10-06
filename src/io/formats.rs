use crate::{bed::BedFormat, io::vcf_writer::VcfFormat};
use clio::ClioPath;
use color_eyre::{
    Result, Section as _,
    eyre::{ContextCompat as _, bail, eyre},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
    VcfLike(VcfFormat),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    VcfLike(VcfFormat),
    Bed(BedFormat),
}

impl clap::ValueEnum for InputFormat {
    fn value_variants<'a>() -> &'a [Self] {
        &[
            InputFormat::VcfLike(VcfFormat::Vcf),
            InputFormat::VcfLike(VcfFormat::Bcf),
            InputFormat::VcfLike(VcfFormat::VcfCompressed),
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
            OutputFormat::VcfLike(VcfFormat::Vcf),
            OutputFormat::VcfLike(VcfFormat::Bcf),
            OutputFormat::VcfLike(VcfFormat::VcfCompressed),
            OutputFormat::Bed(BedFormat::Bed),
            OutputFormat::Bed(BedFormat::BedGz),
        ]
    }

    fn to_possible_value(&self) -> Option<clap::builder::PossibleValue> {
        match self {
            OutputFormat::VcfLike(format) => format.to_possible_value(),
            OutputFormat::Bed(format) => format.to_possible_value(),
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

impl FromFileExtension for VcfFormat {
    fn from_file_extension(p: &str) -> Option<Self> {
        if p.ends_with(".bcf") {
            Some(VcfFormat::Bcf)
        } else if p.ends_with(".vcf.gz") {
            Some(VcfFormat::VcfCompressed)
        } else if p.ends_with(".vcf") {
            Some(VcfFormat::Vcf)
        } else {
            None
        }
    }
}

impl FromFileExtension for BedFormat {
    fn from_file_extension(p: &str) -> Option<Self> {
        if p.ends_with(".bed.gz") {
            Some(BedFormat::BedGz)
        } else if p.ends_with(".bed") {
            Some(BedFormat::Bed)
        } else {
            None
        }
    }
}

impl FromFileExtension for InputFormat {
    fn from_file_extension(p: &str) -> Option<Self> {
        VcfFormat::from_file_extension(p).map(InputFormat::VcfLike)
    }
}

impl FromFileExtension for OutputFormat {
    fn from_file_extension(p: &str) -> Option<Self> {
        VcfFormat::from_file_extension(p)
            .map(OutputFormat::VcfLike)
            .or_else(|| BedFormat::from_file_extension(p).map(OutputFormat::Bed))
    }
}
