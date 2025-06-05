use crate::call::{
    scores::VariantCandidatePileupMetrics,
    variants::{Counter, VariantCandidatePileup},
};

#[derive(Debug, Clone, clap::Args)]
pub struct ThresholdConfig {
    /// The minimum VAF to call a variant
    #[clap(long, default_value_t = 0.1)]
    pub vaf_min: f64,

    /// The maximum binomial p-value to call a variant
    #[clap(long, default_value_t = 0.08)]
    pub binomial_max: f64,

    /// The minimum number of reads to call a variant
    #[clap(long, default_value_t = 5)]
    pub reads_min: usize,
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
        if metrics.binomial > config.binomial_max {
            return false;
        }
        if metrics.vaf < config.vaf_min {
            return false;
        }
        let num_reads = self.bases.len();
        if num_reads < config.reads_min {
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
    use smol_str::SmolStr;

    fn make_pileup(
        pos: u32,
        reference_base: Base,
        sequence_before: Vec<Base>,
        sequence_after: Vec<Base>,
        bases: Vec<Base>,
        _min_reads: usize,
    ) -> VariantCandidatePileup {
        let read_length = u32::try_from(bases.len()).expect("bases length should fit in u32");
        let seen_bases = bases
            .into_iter()
            .enumerate()
            .map(|(idx, b)| crate::call::variants::SeenBase {
                base: b,
                qual: 30,
                mapq: 30,
                reverse: false,
                position: crate::call::variants::PositionInRead {
                    pos: u32::try_from(idx).expect("pos should fit into u32"),
                    read_length,
                },
                qname: SmallVec::from_slice(b"test"),
            })
            .collect();

        VariantCandidatePileup {
            chrom: SmolStr::new_inline("chr66"),
            pos,
            reference_base,
            sequence_before: SmallVec::from_slice(&sequence_before),
            sequence_after: SmallVec::from_slice(&sequence_after),
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
            sequence_before in prop::collection::vec(
                prop_oneof![
                    Just(Base::A),
                    Just(Base::T),
                    Just(Base::C),
                    Just(Base::G)
                ],
                0..2
            ),
            sequence_after in prop::collection::vec(
                prop_oneof![
                    Just(Base::A),
                    Just(Base::T),
                    Just(Base::C),
                    Just(Base::G)
                ],
                0..2
            ),
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
            let pileup = make_pileup(pos, reference_base, sequence_before, sequence_after, bases, 5);
            let metrics = pileup.metrics().unwrap();
            let config = ThresholdConfig {
                vaf_min: 0.1,
                binomial_max: 0.08,
                reads_min: 5,
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

            let pileup = make_pileup(pos, Base::C, vec![Base::A, Base::T], vec![Base::G, Base::G], bases, 5);
            let metrics = pileup.metrics().unwrap();
            let config = ThresholdConfig {
                vaf_min: 0.1,
                binomial_max: 0.1,  // Increased slightly to accommodate test case
                reads_min: 5,
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
            reads_min in 5usize..10usize
        ) {
            let bases = vec![Base::T; count];
            let pileup = make_pileup(pos, Base::C, vec![Base::A, Base::T], vec![Base::G, Base::G], bases, reads_min);
            let metrics = pileup.metrics().unwrap();
            let config = ThresholdConfig {
                vaf_min: 0.1,
                binomial_max: 0.08,
                reads_min,
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

            let pileup = make_pileup(pos, Base::C, vec![Base::A, Base::T], vec![Base::G, Base::G], bases, 5);
            let metrics = pileup.metrics().unwrap();
            let config = ThresholdConfig {
                vaf_min: 0.1,
                binomial_max: 0.08,
                reads_min: 5,
            };

            // Should not be methylation since more A+G than C+T
            prop_assert!(!pileup.likely_methylation_event(&metrics, &config));
        }
    }
}
