//! Methylation detection for variant candidates.
//!
//! Options:
//! 1. methylated
//! 2. snip
//! 3. both

use crate::{call::variants::VariantCandidatePileup, utils::Base};

pub mod threshold;

impl VariantCandidatePileup {
    /// Is this a C->G variant candidate?
    pub fn is_cpg(&self) -> bool {
        self.reference_base == Base::C && self.sequence_after().first() == Some(&Base::G)
    }

    /// With TAPS, methylation is detected by observing a C->T transition
    pub fn could_be_methylation_event(&self) -> bool {
        self.bases.iter().any(|b| b.base == Base::T)
    }
}
