use crate::call::{
    scores::VariantCandidatePileupMetrics,
    variants::{Counter, VariantCandidatePileup},
};

impl VariantCandidatePileup {
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
}
