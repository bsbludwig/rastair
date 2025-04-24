use super::variants::VariantCandidatePileup;
use crate::utils::Base;

impl VariantCandidatePileup {
    /// With TAPS, methylation is detected by observing a C->T transition
    pub fn could_be_methylation_event(&self) -> bool {
        self.is_cpg() && self.bases.iter().any(|b| b.base == Base::T)
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
