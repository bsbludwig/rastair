use crate::utils::IntoF64;
use color_eyre::{Result, eyre::WrapErr};
use rastair_vcf::{FormatFieldNumber, FormatFieldValue, VcfField as _};
use rust_htslib::bcf::Record;
use seqair_types::{Probability, SmallVec, smallvec::smallvec};
use std::fmt;

/// Whether a CpG is from the reference sequence or created by a de-novo variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum CpgOrigin {
    Original,
    DeNovo,
}

/// Methylation measurement for one CpG allele at a position.
#[derive(Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct CpgBeta {
    pub origin: CpgOrigin,
    pub beta: Probability,
    pub mod_count: u32,
    pub total_count: u32,
}

impl CpgBeta {
    pub fn has_evidence(&self) -> bool {
        self.total_count > 0
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for CpgBeta {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "CpgBeta({:?}, beta={:.3}, {}/{})",
            self.origin, *self.beta, self.mod_count, self.total_count
        )
    }
}

/// Methylation information for a position.
///
/// Empty means the position is not in a CpG context. One or two entries
/// represent original and/or de-novo CpG measurements.
#[derive(Clone, Default, Debug, serde::Serialize, serde::Deserialize)]
#[must_use]
pub struct Methylated(pub SmallVec<CpgBeta, 2>);

impl Methylated {
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn has_evidence(&self) -> bool {
        self.0.iter().any(CpgBeta::has_evidence)
    }

    pub fn original(&self) -> Option<&CpgBeta> {
        self.0.iter().find(|b| b.origin == CpgOrigin::Original)
    }

    pub fn denovo(&self) -> Option<&CpgBeta> {
        self.0.iter().find(|b| b.origin == CpgOrigin::DeNovo)
    }

    pub fn iter(&self) -> impl Iterator<Item = &CpgBeta> {
        self.0.iter()
    }

    /// Extract values in canonical order (Original first, then de-novo) for VCF output.
    fn ordered_values<T>(&self, f: impl Fn(&CpgBeta) -> T) -> SmallVec<Option<T>, 2> {
        if self.is_empty() {
            return smallvec![None];
        }
        let mut out = SmallVec::new();
        if let Some(b) = self.original() {
            out.push(Some(f(b)));
        }
        if let Some(b) = self.denovo() {
            out.push(Some(f(b)));
        }
        out
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
        let values: SmallVec<Option<f64>, 2> = self.ordered_values(|b| b.beta.f());

        <Option<f64> as FormatFieldValue>::write(record, Self::ID, &values).wrap_err_with(|| {
            format!(
                "Failed to write format field {} (type {})",
                Self::ID,
                <Self::Type as FormatFieldValue>::TYPE_NAME
            )
        })
    }
}

/// Total read depth for 5-methylcytosine detection (`mod_count` + `unmod_count`), before het adjustment.
#[derive(Clone, Default, Debug, serde::Serialize, serde::Deserialize)]
pub struct MethylationDepth(pub SmallVec<Option<u32>, 2>);

impl From<&Methylated> for MethylationDepth {
    fn from(m: &Methylated) -> Self {
        Self(m.ordered_values(|b| b.total_count))
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

/// Read depth supporting 5-methylcytosine modification (`mod_count` only), before het adjustment.
#[derive(Clone, Default, Debug, serde::Serialize, serde::Deserialize)]
pub struct MethylationAltDepth(pub SmallVec<Option<u32>, 2>);

impl From<&Methylated> for MethylationAltDepth {
    fn from(m: &Methylated) -> Self {
        Self(m.ordered_values(|b| b.mod_count))
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
