use crate::io::vcf_writer::{self, VcfFormat};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
    VcfLike(vcf_writer::Format),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    VcfLike(vcf_writer::Format),
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

pub trait FromFileExtension {
    fn from_file_extension(p: &str) -> Option<Self>
    where
        Self: Sized;
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
        vcf_writer::Format::from_file_extension(p).map(OutputFormat::VcfLike)
    }
}
