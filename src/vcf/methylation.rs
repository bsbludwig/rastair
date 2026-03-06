use crate::utils::IntoF64;
use color_eyre::{Result, eyre::WrapErr};
use rastair_vcf::{FormatFieldNumber, FormatFieldValue, VcfField as _};
use rust_htslib::bcf::Record;
use seqair_types::{Probability, SmallVec, smallvec::smallvec};
use std::fmt;

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
    OriginalCpG { beta: Probability, mod_count: u32, total_count: u32 },
    /// De-novo CpG site
    DeNovoCpG { beta: Probability, mod_count: u32, total_count: u32 },
    /// Both original CpG and de-novo CpG at this position
    Both {
        original_beta: Probability,
        original_mod_count: u32,
        original_total_count: u32,
        denovo_beta: Probability,
        denovo_mod_count: u32,
        denovo_total_count: u32,
    },
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for Methylated {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Methylated::Unknown => f.debug_tuple("Methylated::Unknown").finish(),
            Methylated::NoEvidence => f.debug_tuple("Methylated::NoEvidence").finish(),
            Methylated::OriginalCpG { beta, .. } => {
                f.debug_tuple("Methylated::OriginalCpG").field(beta).finish()
            }
            Methylated::DeNovoCpG { beta, .. } => {
                f.debug_tuple("Methylated::DeNovoCpG").field(beta).finish()
            }
            Methylated::Both { original_beta, denovo_beta, .. } => {
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
        let values: &[Option<f64>] = match self {
            // Unknown: no processing was done, write None
            Methylated::Unknown => &[None],
            // NoEvidence: we checked and found no methylation, write 0.0
            Methylated::NoEvidence => &[Some(0.0)],
            // Single context: write one beta value
            Methylated::OriginalCpG { beta, .. } | Methylated::DeNovoCpG { beta, .. } => {
                &[Some(beta.f())]
            }
            // Dual context: write both beta values (original first, de-novo second)
            Methylated::Both { original_beta, denovo_beta, .. } => {
                &[Some(original_beta.f()), Some(denovo_beta.f())]
            }
        };

        <Option<f64> as FormatFieldValue>::write(record, Self::ID, values).wrap_err_with(|| {
            format!(
                "Failed to write format field {} (type {})",
                Self::ID,
                <Self::Type as FormatFieldValue>::TYPE_NAME
            )
        })
    }
}

/// Total read depth for 5-methylcytosine detection (mod_count + unmod_count), before het adjustment.
#[derive(Clone, Default, Debug, serde::Serialize, serde::Deserialize)]
pub struct MethylationDepth(pub SmallVec<Option<u32>, 2>);

impl From<&Methylated> for MethylationDepth {
    fn from(m: &Methylated) -> Self {
        let values: SmallVec<Option<u32>, 2> = match m {
            Methylated::Unknown => smallvec![None],
            Methylated::NoEvidence => smallvec![Some(0)],
            Methylated::OriginalCpG { total_count, .. }
            | Methylated::DeNovoCpG { total_count, .. } => smallvec![Some(*total_count)],
            Methylated::Both { original_total_count, denovo_total_count, .. } => {
                smallvec![Some(*original_total_count), Some(*denovo_total_count)]
            }
        };
        Self(values)
    }
}

impl rastair_vcf::VcfField for MethylationDepth {
    const ID: &'static cstr8::CStr8 = cstr8::cstr8!("DPM5mC");
}

impl rastair_vcf::HeaderField for MethylationDepth {
    const DESCRIPTION: &'static str = "Total read depth for 5-methylcytosine detection";
}

impl rastair_vcf::FormatField for MethylationDepth {
    type Type = Option<u32>;
    const NUMBER: FormatFieldNumber = FormatFieldNumber::OnePerPossibleBaseModification;

    fn write(&self, record: &mut Record) -> Result<()> {
        <Option<u32> as FormatFieldValue>::write(record, Self::ID, &self.0)
            .wrap_err("Failed to write DPM5mC")
    }
}

/// Read depth supporting 5-methylcytosine modification (mod_count only), before het adjustment.
#[derive(Clone, Default, Debug, serde::Serialize, serde::Deserialize)]
pub struct MethylationAltDepth(pub SmallVec<Option<u32>, 2>);

impl From<&Methylated> for MethylationAltDepth {
    fn from(m: &Methylated) -> Self {
        let values: SmallVec<Option<u32>, 2> = match m {
            Methylated::Unknown => smallvec![None],
            Methylated::NoEvidence => smallvec![Some(0)],
            Methylated::OriginalCpG { mod_count, .. } | Methylated::DeNovoCpG { mod_count, .. } => {
                smallvec![Some(*mod_count)]
            }
            Methylated::Both { original_mod_count, denovo_mod_count, .. } => {
                smallvec![Some(*original_mod_count), Some(*denovo_mod_count)]
            }
        };
        Self(values)
    }
}

impl rastair_vcf::VcfField for MethylationAltDepth {
    const ID: &'static cstr8::CStr8 = cstr8::cstr8!("ADM5mC");
}

impl rastair_vcf::HeaderField for MethylationAltDepth {
    const DESCRIPTION: &'static str = "Read depth supporting 5-methylcytosine modification";
}

impl rastair_vcf::FormatField for MethylationAltDepth {
    type Type = Option<u32>;
    const NUMBER: FormatFieldNumber = FormatFieldNumber::OnePerPossibleBaseModification;

    fn write(&self, record: &mut Record) -> Result<()> {
        <Option<u32> as FormatFieldValue>::write(record, Self::ID, &self.0)
            .wrap_err("Failed to write ADM5mC")
    }
}
