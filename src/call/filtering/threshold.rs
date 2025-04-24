use crate::call::{
    scores::VariantCandidatePileupMetrics,
    variants::{Counter, VariantCandidatePileup},
};

#[derive(Debug, clap::Args)]
pub struct ThresholdConfig {
    /// The minimum VAF to call a variant
    #[clap(long, default_value_t = 0.1)]
    pub min_vaf: f64,

    /// The maximum binomial p-value to call a variant
    #[clap(long, default_value_t = 0.08)]
    pub max_binomial: f64,

    /// The minimum number of reads to call a variant
    #[clap(long, default_value_t = 5)]
    pub min_reads: usize,
}

impl VariantCandidatePileup {
    /// arbitrary filters
    pub fn likely_methylation_event(
        &self,
        metrics: &VariantCandidatePileupMetrics,
        config: &ThresholdConfig,
    ) -> bool {
        if !self.could_be_methylation_event() {
            return false;
        }
        if metrics.binomial > config.max_binomial {
            return false;
        }
        if metrics.vaf < config.min_vaf {
            return false;
        }
        let num_reads = self.bases.len();
        if num_reads < config.min_reads {
            return false;
        }

        // if the methylation evidence is less than half of the total evidence, we don't call it
        let counter = Counter::from_iter(self.bases.iter().map(|b| b.base));
        let methylation_evidence = counter.t + counter.c;
        // there are more A and G than C and T
        if methylation_evidence < (num_reads - methylation_evidence) {
            return false;
        }

        true
    }
}
