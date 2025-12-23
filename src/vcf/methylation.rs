use std::fmt;

use color_eyre::{Result, eyre::WrapErr};
use rastair_types::{Probability, SmallVec, smallvec::smallvec};
use rastair_vcf::{FormatFieldNumber, FormatFieldValue, VcfField as _};
use rust_htslib::bcf::Record;

use crate::utils::IntoF64;

/// Methylation information
#[derive(Clone, Default, serde::Serialize, serde::Deserialize)]
#[must_use]
pub enum Methylated {
    /// Unknown methylation status, i.e., no processing was done
    #[default]
    Unknown,
    /// No evidence of methylation at this site
    NoEvidence,
    /// Original CpG site
    OriginalCpG { beta: Probability },
    /// De-novo CpG site
    DeNovoCpG { beta: Probability },
    /// Both original CpG and de-novo CpG at this position
    Both { original_beta: Probability, denovo_beta: Probability },
}

impl Methylated {
    /// Returns the beta value(s) as a vector.
    /// For positions with both original and de-novo CpG, returns both values.
    pub fn betas(&self) -> SmallVec<Probability, 2> {
        match self {
            Methylated::Unknown => smallvec![],
            Methylated::NoEvidence => smallvec![Probability::ZERO],
            Methylated::OriginalCpG { beta } | Methylated::DeNovoCpG { beta } => smallvec![*beta],
            Methylated::Both { original_beta, denovo_beta } => {
                smallvec![*original_beta, *denovo_beta]
            }
        }
    }

    /// Returns the first beta value, for backwards compatibility.
    pub fn beta(&self) -> Option<Probability> {
        match self {
            Methylated::Unknown => None,
            Methylated::NoEvidence => Some(Probability::ZERO),
            Methylated::OriginalCpG { beta } | Methylated::DeNovoCpG { beta } => Some(*beta),
            Methylated::Both { original_beta, .. } => Some(*original_beta),
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
            Methylated::Both { original_beta, denovo_beta } => {
                f.debug_tuple("Methylated::Both").field(original_beta).field(denovo_beta).finish()
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
    const NUMBER: FormatFieldNumber = FormatFieldNumber::OnePerPossibleBaseModification;

    fn write(&self, record: &mut Record) -> Result<()> {
        let betas = self.betas();
        let values: Vec<Option<f64>> =
            if betas.is_empty() { vec![None] } else { betas.iter().map(|b| Some(b.f())).collect() };

        <Option<f64> as FormatFieldValue>::write(record, Self::ID, &values).wrap_err_with(|| {
            format!(
                "Failed to write format field {} (type {})",
                Self::ID,
                <Self::Type as FormatFieldValue>::TYPE_NAME
            )
        })
    }
}
