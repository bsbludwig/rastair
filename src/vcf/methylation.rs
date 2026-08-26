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
    pub(crate) fn ordered_values<T>(&self, f: impl Fn(&CpgBeta) -> T) -> SmallVec<Option<T>, 2> {
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

/// Total read depth for 5-methylcytosine detection (`mod_count` + `unmod_count`), before het adjustment.
#[derive(Clone, Default, Debug, serde::Serialize, serde::Deserialize)]
pub struct MethylationDepth(pub SmallVec<Option<u32>, 2>);

impl From<&Methylated> for MethylationDepth {
    fn from(m: &Methylated) -> Self {
        Self(m.ordered_values(|b| b.total_count))
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
