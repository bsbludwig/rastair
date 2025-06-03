use std::fmt;

mod format_field;
mod info_field;

pub use format_field::{FormatField, FormatFieldValue};
pub use info_field::{InfoField, InfoFieldValue};

/// A field that can be used in the header section.
pub trait HeaderField {
    /// Description of the field, used in the VCF header.
    const DESCRIPTION: &'static str;
}

/// A field that can be used in VCF.
pub trait VcfField: HeaderField {
    /// The ID of the field, used in the VCF header and record.
    const ID: &'static str;
    /// The number of values that can be included with the field.
    const NUMBER: FieldNumber;
}

/// The number of values that can be included with the INFO field
pub enum FieldNumber {
    /// A flag, no values. Written as "0" in the header.
    Flag,
    /// A fixed number of values, must be non-zero. Written as a number in the header.
    Num(u32),
    /// One value per alternative allele. Written as "A" in the header.
    OneValPerAlt,
    /// One value per alternative allele and reference allele. Written as "R" in the header.
    OnePerAltAndRef,
    /// One value per genotype. Written as "G" in the header.
    OnePerGenotype,
    /// Variable number of values or unknown or unbounded. Written as "." in the header.
    Dot,
}

impl fmt::Display for FieldNumber {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FieldNumber::Flag => write!(f, "0"),
            FieldNumber::Num(n) => write!(f, "{n}"),
            FieldNumber::OneValPerAlt => write!(f, "A"),
            FieldNumber::OnePerAltAndRef => write!(f, "R"),
            FieldNumber::OnePerGenotype => write!(f, "G"),
            FieldNumber::Dot => write!(f, "."),
        }
    }
}
