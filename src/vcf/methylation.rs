use std::fmt;

use color_eyre::{Result, eyre::WrapErr};
use rastair_vcf::{FormatFieldNumber, FormatFieldValue, VcfField as _};
use rust_htslib::bcf::Record;

/// Methylation information
#[derive(Clone, Default, serde::Serialize, serde::Deserialize)]
pub enum Methylated {
    /// Unknown methylation status, i.e., no processing was done
    #[default]
    Unknown,
    /// No evidence of methylation at this site
    NoEvidence,
    /// Original CpG site
    OriginalCpG { beta: f64 },
    /// De-novo CpG site
    DeNovoCpG { beta: f64 },
}

impl Methylated {
    pub fn beta(&self) -> Option<f64> {
        match self {
            Methylated::Unknown => None,
            Methylated::NoEvidence => Some(0.0),
            Methylated::OriginalCpG { beta } => Some(*beta),
            Methylated::DeNovoCpG { beta } => Some(*beta),
        }
    }
}

impl fmt::Debug for Methylated {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Methylated::Unknown => f.debug_tuple("Methylated::Unknown").finish(),
            Methylated::NoEvidence => f.debug_tuple("Methylated::NoEvidence").finish(),
            Methylated::OriginalCpG { beta } => {
                f.debug_tuple("Methylated::OriginalCpG").field(beta).finish()
            }
            Methylated::DeNovoCpG { beta } => {
                f.debug_tuple("Methylated::DeNovoCpG").field(beta).finish()
            }
        }
    }
}

impl rastair_vcf::VcfField for Methylated {
    const ID: &'static cstr8::CStr8 = cstr8::cstr8!("M5mC");
}

impl rastair_vcf::HeaderField for Methylated {
    const DESCRIPTION: &'static str = "Methylation level at CpG sites";
}

impl rastair_vcf::FormatField for Methylated {
    type Type = Option<f64>;
    const NUMBER: FormatFieldNumber = FormatFieldNumber::Num(1);

    fn write(&self, record: &mut Record) -> Result<()> {
        <Option<f64> as FormatFieldValue>::write(record, Self::ID, &[self.beta()]).wrap_err_with(
            || {
                format!(
                    "Failed to write format field {} (type {})",
                    Self::ID,
                    <Self::Type as FormatFieldValue>::TYPE_NAME
                )
            },
        )
    }
}
