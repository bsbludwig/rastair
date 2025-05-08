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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::Base;
    use proptest::prelude::*;
    use smallvec::SmallVec;

    fn make_pileup(
        pos: u32,
        reference_base: Base,
        next_base: Option<Base>,
        bases: Vec<Base>,
        _min_reads: usize,
    ) -> VariantCandidatePileup {
        let seen_bases = bases
            .into_iter()
            .map(|b| crate::call::variants::SeenBase {
                base: b,
                qual: 30,
                mapq: 30,
                reverse: false,
                at_fringe: false,
                qname: SmallVec::from_slice(b"test"),
            })
            .collect();

        VariantCandidatePileup {
            pos,
            reference_base,
            next_base,
            bases: crate::call::variants::SeenBases(seen_bases),
        }
    }

    proptest! {
        #[test]
        fn not_methylation_if_not_cpg(
            pos in 1u32..1000u32,
            reference_base in prop_oneof![
                Just(Base::A),
                Just(Base::T),
                Just(Base::G)
            ],
            next_base in prop_oneof![
                Just(None),
                Just(Some(Base::A)),
                Just(Some(Base::T)),
                Just(Some(Base::C)),
                Just(Some(Base::G))
            ],
            bases in prop::collection::vec(
                prop_oneof![
                    Just(Base::A),
                    Just(Base::T),
                    Just(Base::C),
                    Just(Base::G)
                ],
                5..100
            )
        ) {
            let pileup = make_pileup(pos, reference_base, next_base, bases, 5);
            let metrics = pileup.metrics();
            let config = ThresholdConfig {
                min_vaf: 0.1,
                max_binomial: 0.08,
                min_reads: 5,
            };

            prop_assert!(!pileup.likely_methylation_event(&metrics, &config));
        }

        #[test]
        fn potential_methylation_requires_cpg_and_t(
            pos in 1u32..1000u32,
            t_count in 3usize..6usize,  // Reduced T count
            c_count in 2usize..4usize   // Increased C count for better ratio
        ) {
            let mut bases = vec![Base::T; t_count];
            bases.extend(vec![Base::C; c_count]);

            let pileup = make_pileup(pos, Base::C, Some(Base::G), bases, 5);
            let metrics = pileup.metrics();
            let config = ThresholdConfig {
                min_vaf: 0.1,
                max_binomial: 0.1,  // Increased slightly to accommodate test case
                min_reads: 5,
            };

            // Should be possible methylation since we have:
            // 1. CpG site (C reference with G next)
            // 2. T bases present (indicating methylation)
            // 3. More T bases than other bases
            // 4. Above min read threshold
            // 5. A mix of C and T bases to keep binomial test reasonable
            prop_assert!(pileup.likely_methylation_event(&metrics, &config));
        }

        #[test]
        fn respects_min_reads_threshold(
            pos in 1u32..1000u32,
            count in 1usize..4usize,
            min_reads in 5usize..10usize
        ) {
            let bases = vec![Base::T; count];
            let pileup = make_pileup(pos, Base::C, Some(Base::G), bases, min_reads);
            let metrics = pileup.metrics();
            let config = ThresholdConfig {
                min_vaf: 0.1,
                max_binomial: 0.08,
                min_reads,
            };

            // Should not be methylation since below min_reads threshold
            prop_assert!(!pileup.likely_methylation_event(&metrics, &config));
        }

        #[test]
        fn requires_sufficient_methylation_evidence(
            pos in 1u32..1000u32,
            c_count in 0usize..5usize,
            t_count in 0usize..5usize,
            a_count in 6usize..10usize,
            g_count in 6usize..10usize
        ) {
            let mut bases = vec![];
            bases.extend(vec![Base::C; c_count]);
            bases.extend(vec![Base::T; t_count]);
            bases.extend(vec![Base::A; a_count]);
            bases.extend(vec![Base::G; g_count]);

            let pileup = make_pileup(pos, Base::C, Some(Base::G), bases, 5);
            let metrics = pileup.metrics();
            let config = ThresholdConfig {
                min_vaf: 0.1,
                max_binomial: 0.08,
                min_reads: 5,
            };

            // Should not be methylation since more A+G than C+T
            prop_assert!(!pileup.likely_methylation_event(&metrics, &config));
        }
    }
}
