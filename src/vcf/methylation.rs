use seqair_types::Probability;
use std::{fmt, num::NonZeroU32};

/// Whether a CpG is from the reference sequence or created by a de-novo variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CpgOrigin {
    Original,
    DeNovo,
}

/// Methylation measured at one CpG of a position.
#[derive(Clone, Copy)]
pub struct CpgBeta {
    pub beta: Probability,
    pub mod_count: u32,
    /// A CpG without reads has no measurement.
    pub total_count: NonZeroU32,
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for CpgBeta {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "CpgBeta(beta={:.3}, {}/{})", *self.beta, self.mod_count, self.total_count)
    }
}

/// The methylation measured at a position, for the reference's CpG and for
/// one a called alt makes; neither for a position in no CpG or without reads.
#[derive(Clone, Copy, Default, Debug)]
#[must_use]
pub struct Methylated {
    pub original: Option<CpgBeta>,
    pub denovo: Option<CpgBeta>,
}

impl Methylated {
    pub fn get(&self, origin: CpgOrigin) -> Option<&CpgBeta> {
        match origin {
            CpgOrigin::Original => self.original.as_ref(),
            CpgOrigin::DeNovo => self.denovo.as_ref(),
        }
    }

    pub(crate) fn get_mut(&mut self, origin: CpgOrigin) -> &mut Option<CpgBeta> {
        match origin {
            CpgOrigin::Original => &mut self.original,
            CpgOrigin::DeNovo => &mut self.denovo,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.original.is_none() && self.denovo.is_none()
    }
}
