use crate::utils::Base;

use super::{
    scores::VariantCandidatePileupMetrics,
    variants::{Counter, VariantCandidatePileup},
};

impl VariantCandidatePileup {
    /// With TAPS, methylation is detected by observing a C->T transition
    pub fn could_be_methylation_event(&self) -> bool {
        self.is_cpg() && self.bases.iter().any(|b| b.base == Base::T)
    }

    /// arbitrary filters
    pub fn likely_methylation_event(&self, metrics: &VariantCandidatePileupMetrics) -> bool {
        if !self.could_be_methylation_event() {
            return false;
        }
        if metrics.binomial > 0.08 {
            return false;
        }
        if metrics.vaf < 0.1 {
            return false;
        }

        let counter = Counter::from_iter(self.bases.iter().map(|b| b.base));
        let methylation_evidence = counter.t + counter.c;
        // there are more A and G than C and T
        if methylation_evidence < (self.bases.len() - methylation_evidence) {
            return false;
        }

        true
    }

    pub fn beta(&self) -> f64 {
        let total_evidence = self.bases.len();
        if total_evidence == 0 {
            return 0.0;
        }
        let methylation_evidence = self.bases.iter().filter(|b| b.base == Base::T).count() as f64;
        methylation_evidence / total_evidence as f64
    }
}
