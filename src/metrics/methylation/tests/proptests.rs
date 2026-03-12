use std::num::NonZeroU8;

use super::*;
use crate::{
    call::{
        pileup::{Pileup, SimpleRead, SimpleReads},
        variant_calling::{EstimatedGenotype, GenotypeTag},
    },
    metrics::AltCall,
    sequence::{ChunkRegion, Region},
    utils::default,
    vcf::{CpgOrigin, SequenceContext},
};
use proptest::prelude::*;
use rastair_types::{Base::*, Probability};

// ---------------------------------------------------------------------------
// Scenario: the *desired outcome* that we generate first, then build reads for
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct MethylationScenario {
    side: CpgSide,
    origin: CpgOrigin,
    /// Number of methylated-looking reads on the informative strand.
    mod_reads: u8,
    /// Number of unmethylated-looking reads on the informative strand.
    unmod_reads: u8,
    adjustment: Adjustment,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Adjustment {
    None,
    HetConfounded,
    /// Only valid for Original CpGs.
    HomAlt,
}

impl MethylationScenario {
    /// What `call()` should return for this scenario.
    fn expected_beta(&self) -> Option<f64> {
        let m = self.mod_reads as f64;
        let u = self.unmod_reads as f64;
        if m + u == 0.0 {
            return None;
        }
        Some(match self.adjustment {
            Adjustment::HomAlt => 0.0,
            Adjustment::HetConfounded => {
                let total = m + u;
                let excess = (m - total / 2.0).max(0.0);
                excess / (u + excess)
            }
            Adjustment::None => m / (m + u),
        })
    }

    /// Construct a `PileupMetrics` that embodies this scenario.
    fn build(&self) -> PileupMetrics {
        let (ref_base, context) = match (self.side, self.origin) {
            (CpgSide::C, CpgOrigin::Original) => (C, ctx_after_g(C)),
            (CpgSide::G, CpgOrigin::Original) => (G, ctx_before_c(G)),
            (CpgSide::C, CpgOrigin::DeNovo) => (T, ctx_after_g(T)),
            (CpgSide::G, CpgOrigin::DeNovo) => (A, ctx_before_c(A)),
        };

        let strand = self.side.strand();
        let mod_base = self.side.mod_base();
        let unmod_base = self.side.unmod_base();
        let (before_base, after_base) = match self.side {
            CpgSide::C => (None, Some(G)),
            CpgSide::G => (Some(C), None),
        };

        let make_read = |base| SimpleRead {
            base,
            strand,
            before_base,
            after_base,
            ..default()
        };

        let mut reads: Vec<SimpleRead> = Vec::new();
        reads.extend((0..self.mod_reads).map(|_| make_read(mod_base)));
        reads.extend((0..self.unmod_reads).map(|_| make_read(unmod_base)));

        let pileup = Pileup {
            region: dummy_region(),
            context,
            pos: 1000,
            reads: SimpleReads(reads.into()),
            reference_base: ref_base,
        };

        let mut metrics = PileupMetrics::new(pileup).unwrap();

        // For de-novo: mark the new CpG-forming alt as a real variant.
        if self.origin == CpgOrigin::DeNovo {
            let denovo_base = self.side.unmod_base();
            if let Some(alt) = metrics.alts.iter_mut().find(|a| a.base == denovo_base) {
                alt.call = AltCall::RealVariant;
            }
        }

        metrics.pos_metrics.extended.genotype = self.build_genotype();
        metrics
    }

    fn build_genotype(&self) -> Option<EstimatedGenotype> {
        let gt = match self.adjustment {
            Adjustment::None => return None,
            Adjustment::HetConfounded => GenotypeTag::RefHet(NonZeroU8::new(1).unwrap()),
            Adjustment::HomAlt => GenotypeTag::HomAlt(NonZeroU8::new(1).unwrap()),
        };
        Some(EstimatedGenotype {
            genotype: gt,
            likelihood: Probability::new(0.99).unwrap(),
            confidence: Probability::new(0.99).unwrap(),
        })
    }

    /// Mirror this scenario to the opposite CpgSide.
    fn mirror(&self) -> Self {
        Self {
            side: match self.side {
                CpgSide::C => CpgSide::G,
                CpgSide::G => CpgSide::C,
            },
            ..self.clone()
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers for building minimal Pileup / PileupMetrics by hand
// ---------------------------------------------------------------------------

fn dummy_region() -> ChunkRegion {
    ChunkRegion {
        region: Region { contig: "chr_test".into(), start: 1000, end: 1002 },
        last_position: 1002,
        overlap_start: 0,
        overlap_end: 0,
    }
}

fn ctx_after_g(me: Base) -> SequenceContext {
    SequenceContext { me, after_1: Some(G), ..default() }
}

fn ctx_before_c(me: Base) -> SequenceContext {
    SequenceContext { me, before_1: Some(C), ..default() }
}

// ---------------------------------------------------------------------------
// Strategy
// ---------------------------------------------------------------------------

fn scenario_strategy() -> impl Strategy<Value = MethylationScenario> {
    (
        prop_oneof![Just(CpgSide::C), Just(CpgSide::G)],
        prop_oneof![Just(CpgOrigin::Original), Just(CpgOrigin::DeNovo)],
        0u8..50,
        0u8..50,
    )
        .prop_filter("DeNovo needs the CpG-forming alt to exist", |&(_, origin, _, unmod)| {
            origin == CpgOrigin::Original || unmod >= 1
        })
        .prop_flat_map(|(side, origin, mod_reads, unmod_reads)| {
            let mut adjustments = vec![Adjustment::None];
            // HetConfounded needs the confounding base to exist as an alt
            // (for Original) or as the ref (for DeNovo). Require mod_reads >= 1
            // so the alt is present in either case.
            if mod_reads >= 1 {
                adjustments.push(Adjustment::HetConfounded);
                if origin == CpgOrigin::Original {
                    adjustments.push(Adjustment::HomAlt);
                }
            }
            (
                Just(side),
                Just(origin),
                Just(mod_reads),
                Just(unmod_reads),
                proptest::sample::select(adjustments),
            )
        })
        .prop_map(|(side, origin, mod_reads, unmod_reads, adjustment)| {
            MethylationScenario { side, origin, mod_reads, unmod_reads, adjustment }
        })
}

/// Strategy that only produces scenarios with evidence (mod + unmod > 0).
fn scenario_with_evidence() -> impl Strategy<Value = MethylationScenario> {
    scenario_strategy().prop_filter("need at least one read", |s| {
        s.mod_reads + s.unmod_reads > 0
    })
}

// ---------------------------------------------------------------------------
// Property tests
// ---------------------------------------------------------------------------

proptest! {
    /// The core roundtrip: build reads from a scenario, call methylation,
    /// and verify the result matches the expected beta.
    #[test]
    fn roundtrip(scenario in scenario_strategy()) {
        let metrics = scenario.build();
        let result = call(&metrics).unwrap();

        match scenario.expected_beta() {
            None => {
                prop_assert!(
                    result.as_ref().is_none_or(|m| !m.has_evidence()),
                    "Expected no evidence for {scenario:?}, got {result:?}"
                );
            }
            Some(expected) => {
                let methylated = result.as_ref().expect(&format!(
                    "Expected Some(Methylated) for {scenario:?}"
                ));
                let cpg = methylated
                    .iter()
                    .find(|b| b.origin == scenario.origin)
                    .expect(&format!(
                        "Expected {:?} CpG in {methylated:?} for {scenario:?}",
                        scenario.origin
                    ));
                let actual = *cpg.beta;
                prop_assert!(
                    (actual - expected).abs() < 1e-9,
                    "Beta mismatch for {scenario:?}: expected {expected}, got {actual}"
                );
            }
        }
    }

    /// Beta is always a valid probability [0, 1].
    /// (This is enforced by the Probability type, but let's verify the math
    /// never tries to produce an out-of-range value.)
    #[test]
    fn beta_always_valid_probability(scenario in scenario_with_evidence()) {
        let metrics = scenario.build();
        let result = call(&metrics).unwrap();

        if let Some(methylated) = &result {
            for cpg in methylated.iter() {
                let b = *cpg.beta;
                prop_assert!(
                    (0.0..=1.0).contains(&b),
                    "Beta {b} out of [0,1] for {scenario:?}"
                );
            }
        }
    }

    /// C-side and G-side are symmetric: mirroring a scenario should produce
    /// the same beta value.
    #[test]
    fn c_g_symmetry(scenario in scenario_with_evidence()) {
        let mirror = scenario.mirror();

        let result_a = call(&scenario.build()).unwrap();
        let result_b = call(&mirror.build()).unwrap();

        let beta_a = result_a.as_ref()
            .and_then(|m| m.iter().find(|b| b.origin == scenario.origin))
            .map(|b| *b.beta);
        let beta_b = result_b.as_ref()
            .and_then(|m| m.iter().find(|b| b.origin == mirror.origin))
            .map(|b| *b.beta);

        prop_assert_eq!(
            beta_a.map(|b| (b * 1e9).round()),
            beta_b.map(|b| (b * 1e9).round()),
            "Symmetry broken: {:?} → {:?}, mirror {:?} → {:?}",
            scenario, beta_a, mirror, beta_b
        );
    }

    /// Het-confounded adjustment always produces a beta ≤ the unadjusted beta.
    #[test]
    fn het_confounded_reduces_beta(
        side in prop_oneof![Just(CpgSide::C), Just(CpgSide::G)],
        origin in prop_oneof![Just(CpgOrigin::Original), Just(CpgOrigin::DeNovo)],
        mod_reads in 1u8..50,
        unmod_reads in 1u8..50,
    ) {
        let unadjusted = MethylationScenario {
            side, origin, mod_reads, unmod_reads,
            adjustment: Adjustment::None,
        };
        let confounded = MethylationScenario {
            adjustment: Adjustment::HetConfounded,
            ..unadjusted.clone()
        };

        let beta_plain = call(&unadjusted.build()).unwrap()
            .and_then(|m| m.iter().find(|b| b.origin == origin).map(|b| *b.beta));
        let beta_het = call(&confounded.build()).unwrap()
            .and_then(|m| m.iter().find(|b| b.origin == origin).map(|b| *b.beta));

        if let (Some(plain), Some(het)) = (beta_plain, beta_het) {
            prop_assert!(
                het <= plain + 1e-9,
                "HetConfounded beta ({het}) > unadjusted ({plain}) for {confounded:?}"
            );
        }
    }

    /// HomAlt on an original CpG always produces beta = 0.0.
    #[test]
    fn hom_alt_original_is_zero(
        side in prop_oneof![Just(CpgSide::C), Just(CpgSide::G)],
        mod_reads in 1u8..50,
        unmod_reads in 0u8..50,
    ) {
        let scenario = MethylationScenario {
            side,
            origin: CpgOrigin::Original,
            mod_reads,
            unmod_reads,
            adjustment: Adjustment::HomAlt,
        };

        let result = call(&scenario.build()).unwrap();
        let methylated = result.as_ref().expect("HomAlt with reads should produce result");
        let cpg = methylated.iter().find(|b| b.origin == CpgOrigin::Original)
            .expect("should have Original CpG");
        prop_assert!(
            *cpg.beta < 1e-9,
            "HomAlt Original should be 0.0, got {} for {scenario:?}", *cpg.beta
        );
    }

    /// More methylated reads (with fixed unmethylated count) should give
    /// a higher or equal beta.
    #[test]
    fn monotonic_in_mod_reads(
        side in prop_oneof![Just(CpgSide::C), Just(CpgSide::G)],
        origin in prop_oneof![Just(CpgOrigin::Original), Just(CpgOrigin::DeNovo)],
        mod_lo in 0u8..25,
        mod_extra in 1u8..25,
        // DeNovo requires the CpG-forming alt (unmod_base) to exist
        unmod_reads in 1u8..50,
    ) {
        let mod_hi = mod_lo + mod_extra;

        let lo = MethylationScenario {
            side, origin, mod_reads: mod_lo, unmod_reads,
            adjustment: Adjustment::None,
        };
        let hi = MethylationScenario {
            mod_reads: mod_hi,
            ..lo.clone()
        };

        let beta_lo = call(&lo.build()).unwrap()
            .and_then(|m| m.iter().find(|b| b.origin == origin).map(|b| *b.beta))
            .unwrap_or(0.0);
        let beta_hi = call(&hi.build()).unwrap()
            .and_then(|m| m.iter().find(|b| b.origin == origin).map(|b| *b.beta))
            .unwrap_or(0.0);

        prop_assert!(
            beta_hi >= beta_lo - 1e-9,
            "More mod reads should give higher beta: {mod_lo}→{beta_lo}, {mod_hi}→{beta_hi}"
        );
    }

    /// mod_count and total_count in the result should match the raw read counts.
    #[test]
    fn counts_match_reads(scenario in scenario_with_evidence()) {
        let metrics = scenario.build();
        let result = call(&metrics).unwrap();

        if let Some(methylated) = &result {
            if let Some(cpg) = methylated.iter().find(|b| b.origin == scenario.origin) {
                prop_assert_eq!(
                    cpg.mod_count, u32::from(scenario.mod_reads),
                    "mod_count mismatch for {:?}", scenario
                );
                prop_assert_eq!(
                    cpg.total_count,
                    u32::from(scenario.mod_reads) + u32::from(scenario.unmod_reads),
                    "total_count mismatch for {:?}", scenario
                );
            }
        }
    }
}
