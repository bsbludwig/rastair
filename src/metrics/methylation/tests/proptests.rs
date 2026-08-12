use std::num::NonZeroU8;

use super::*;
use crate::{
    call::{
        pileup::{Pileup, SimpleRead, SimpleReads},
        variant_calling::{EstimatedGenotype, GenotypeTag},
    },
    metrics::{AlleleMetrics, Alt, AltCall, AltFilters, DenovoAdjecent},
    sequence::{ChunkRegion, Region},
    utils::default,
    vcf::{CpgOrigin, SequenceContext},
};
use proptest::prelude::*;
use seqair_types::{Base, Base::*, Probability};

const PROPTEST_CASES: u32 = 2048;

// ---------------------------------------------------------------------------
// Scenario: the *desired outcome* that we generate first, then build reads for
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct MethylationScenario {
    side: CpgSide,
    origin: CpgOrigin,
    /// For `DeNovo`: the original ref base that was mutated.
    /// For Original: ignored (ref is determined by side).
    denovo_ref: Base,
    /// Number of methylated-looking reads on the informative strand.
    mod_reads: u8,
    /// Number of unmethylated-looking reads on the informative strand.
    unmod_reads: u8,
    /// Reads that should be ignored by `read_counts()` — tests the filtering.
    noise: NoiseReads,
    genotype: GenotypeScenario,
}

#[derive(Debug, Clone, Copy, Default)]
struct NoiseReads {
    /// Mod-base reads on the wrong strand (right adj).
    wrong_strand: u8,
    /// Mod-base reads on right strand with wrong adjacent base.
    wrong_adj: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GenotypeScenario {
    /// No genotype set.
    None,
    /// Explicit `HomRef` — should behave identically to None.
    HomRef,
    /// Het with the confounding base.
    HetConfounded,
    /// Het with a non-confounding base — should behave identically to None.
    HetNonConfounding,
    /// `HomAlt` (Original only) — beta forced to 0.0.
    HomAlt,
}

fn empty_pileup() -> Pileup {
    Pileup {
        region: dummy_region(),
        context: SequenceContext::default(),
        pos: 1000,
        reads: SimpleReads(std::sync::Arc::new([])),
        reference_base: Base::Unknown,
        indel_observations: SmallVec::new(),
        depth_offset: 0,
        homopolymer_run: 0,
        dinucleotide_run: 0,
        soft_clip_count: 0,
        indel_ref_window: SmallVec::new(),
        indel_ref_anchor: 0,
    }
}

impl MethylationScenario {
    /// What `call()` should return for this scenario.
    ///
    /// NOTE: This oracle mirrors the production formula. The formula itself is
    /// verified against domain-derived analytical values in the `formula_*`
    /// tests at the bottom of this file. This proptest verifies the *wiring*:
    /// that the right reads reach the right formula with the right adjustment.
    fn expected_beta(&self) -> Option<f64> {
        let m = f64::from(self.mod_reads);
        let u = f64::from(self.unmod_reads);
        if m + u == 0.0 {
            return None;
        }
        Some(match self.effective_adjustment() {
            GenotypeAdjustment::HomAlt => 0.0,
            GenotypeAdjustment::HetConfounded => {
                let total = m + u;
                let excess = (m - total / 2.0).max(0.0);
                excess / (u + excess)
            }
            GenotypeAdjustment::None => m / (m + u),
        })
    }

    fn effective_adjustment(&self) -> GenotypeAdjustment {
        match self.genotype {
            GenotypeScenario::HomAlt if self.origin == CpgOrigin::Original => {
                GenotypeAdjustment::HomAlt
            }
            GenotypeScenario::HetConfounded => GenotypeAdjustment::HetConfounded,
            _ => GenotypeAdjustment::None,
        }
    }

    /// Construct a `PileupMetrics` that embodies this scenario.
    fn build(&self) -> PileupMetrics {
        let (ref_base, context) = self.ref_and_context();

        let strand = self.side.strand();
        let wrong_strand = match strand {
            Strand::OT => Strand::OB,
            Strand::OB => Strand::OT,
            Strand::Unknown => unreachable!(),
        };
        let mod_base = self.side.mod_base();
        let unmod_base = self.side.unmod_base();
        let (before_base, after_base) = self.adj_bases();
        let wrong_adj = Some(T); // T is never the correct adj (C or G)

        let make_read = |base, s, bb, ab| SimpleRead {
            base,
            strand: s,
            before_base: bb,
            after_base: ab,
            ..default()
        };

        let mut reads: Vec<SimpleRead> = Vec::new();

        // Real signal reads
        reads.extend(
            (0..self.mod_reads).map(|_| make_read(mod_base, strand, before_base, after_base)),
        );
        reads.extend(
            (0..self.unmod_reads).map(|_| make_read(unmod_base, strand, before_base, after_base)),
        );

        // Noise: mod-base reads on the wrong strand (should be ignored)
        reads.extend(
            (0..self.noise.wrong_strand)
                .map(|_| make_read(mod_base, wrong_strand, before_base, after_base)),
        );

        // Noise: mod-base reads on right strand but wrong adjacent base
        let (wrong_bb, wrong_ab) = match self.side {
            CpgSide::C => (before_base, wrong_adj),
            CpgSide::G => (wrong_adj, after_base),
        };
        reads.extend(
            (0..self.noise.wrong_adj).map(|_| make_read(mod_base, strand, wrong_bb, wrong_ab)),
        );

        let pileup = Pileup {
            region: dummy_region(),
            context,
            pos: 1000,
            reads: SimpleReads(reads.into()),
            reference_base: ref_base,
            ..empty_pileup()
        };

        let mut metrics = PileupMetrics::new(pileup).unwrap();

        // For de-novo: mark the new CpG-forming alt as a real variant.
        if self.origin == CpgOrigin::DeNovo {
            let denovo_base = self.side.unmod_base();
            if let Some(alt) = metrics.alts.iter_mut().find(|a| a.base == denovo_base) {
                alt.call = AltCall::RealVariant;
            }
        }

        let gt = self.build_genotype(&mut metrics);
        metrics.pos_metrics.extended.genotype = gt;
        metrics
    }

    fn ref_and_context(&self) -> (Base, SequenceContext) {
        match (self.side, self.origin) {
            (CpgSide::C, CpgOrigin::Original) => (C, ctx_after_g(C)),
            (CpgSide::G, CpgOrigin::Original) => (G, ctx_before_c(G)),
            (CpgSide::C, CpgOrigin::DeNovo) => (self.denovo_ref, ctx_after_g(self.denovo_ref)),
            (CpgSide::G, CpgOrigin::DeNovo) => (self.denovo_ref, ctx_before_c(self.denovo_ref)),
        }
    }

    fn adj_bases(&self) -> (Option<Base>, Option<Base>) {
        match self.side {
            CpgSide::C => (Option::None, Some(G)),
            CpgSide::G => (Some(C), Option::None),
        }
    }

    fn build_genotype(&self, metrics: &mut PileupMetrics) -> Option<EstimatedGenotype> {
        let gt = match self.genotype {
            GenotypeScenario::None => return None,
            GenotypeScenario::HomRef => GenotypeTag::HomRef,
            GenotypeScenario::HetConfounded => {
                // For denovo where ref == confounding base, any het triggers
                // confounding via ref_base() check. Use RefHet(1) pointing at
                // the first alt (the CpG-forming base for denovo, or the
                // confounding base for original).
                GenotypeTag::RefHet(NonZeroU8::new(1).unwrap())
            }
            GenotypeScenario::HetNonConfounding => self.build_het_non_confounding(metrics),
            GenotypeScenario::HomAlt => GenotypeTag::HomAlt(NonZeroU8::new(1).unwrap()),
        };
        Some(EstimatedGenotype {
            genotype: gt,
            likelihood: Probability::new(0.99).unwrap(),
            confidence: Probability::new(0.99).unwrap(),
        })
    }

    /// Build a het genotype that does NOT involve the confounding base.
    fn build_het_non_confounding(&self, metrics: &mut PileupMetrics) -> GenotypeTag {
        let confounding = self.side.mod_base();

        // Find an existing non-confounding alt.
        let existing_idx = metrics.alts.iter().position(|a| a.base != confounding);

        let idx = if let Some(idx) = existing_idx {
            idx
        } else {
            // For Original CpGs, unmod reads are ref reads (no alt created).
            // Push a dummy non-confounding alt (e.g., G at a C-ref site).
            let dummy_base = match self.side {
                CpgSide::C => G, // confounding is T, dummy is G
                CpgSide::G => T, // confounding is A, dummy is T
            };
            push_dummy_alt(metrics, dummy_base);
            metrics.alts.len() - 1
        };

        let allele = NonZeroU8::new((idx + 1) as u8).unwrap();
        GenotypeTag::RefHet(allele)
    }

    /// Mirror this scenario to the opposite `CpgSide`.
    fn mirror(&self) -> Self {
        let mirrored_side = match self.side {
            CpgSide::C => CpgSide::G,
            CpgSide::G => CpgSide::C,
        };
        // Mirror the denovo ref to maintain the same relationship.
        // C-side T→C mirrors to G-side A→G, etc.
        let mirrored_denovo_ref = match (self.side, self.denovo_ref) {
            (CpgSide::C, T) => A,       // T→C mirrors to A→G
            (CpgSide::C, A) => T,       // A→C mirrors to T→G
            (CpgSide::C, G) => Base::C, // G→C mirrors to C→G
            (CpgSide::G, A) => T,       // A→G mirrors to T→C
            (CpgSide::G, T) => A,       // T→G mirrors to A→C
            (CpgSide::G, Base::C) => G, // C→G mirrors to G→C
            _ => self.denovo_ref,
        };
        Self { side: mirrored_side, denovo_ref: mirrored_denovo_ref, ..self.clone() }
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

/// Push a dummy alt onto metrics so genotype indices resolve.
fn push_dummy_alt(metrics: &mut PileupMetrics, base: Base) {
    metrics.alts.push(Alt {
        base,
        metrics: AlleleMetrics { base, ..default() },
        filters: AltFilters::default(),
        call: default(),
    });
}

fn extract_beta(result: &Option<Methylated>, origin: CpgOrigin) -> Option<f64> {
    result.as_ref().and_then(|m| m.iter().find(|b| b.origin == origin)).map(|b| *b.beta)
}

// ---------------------------------------------------------------------------
// Strategy
// ---------------------------------------------------------------------------

/// Valid denovo ref bases for a given side (anything except the CpG base itself).
fn denovo_ref_strategy(side: CpgSide) -> impl Strategy<Value = Base> {
    match side {
        // C-side: ref was T, A, or G — mutated to C
        CpgSide::C => prop_oneof![Just(T), Just(A), Just(G)].boxed(),
        // G-side: ref was A, T, or C — mutated to G
        CpgSide::G => prop_oneof![Just(A), Just(T), Just(C)].boxed(),
    }
}

/// Weight boundaries (0, 1) more heavily for better edge case coverage.
fn read_count_strategy(lo: u8, hi: u8) -> impl Strategy<Value = u8> {
    prop_oneof![
        3 => lo..hi,               // uniform bulk
        1 => Just(lo),             // boundary: minimum
        1 => Just(lo.saturating_add(1).min(hi.saturating_sub(1))), // boundary: min+1
    ]
}

fn scenario_strategy() -> impl Strategy<Value = MethylationScenario> {
    let side = prop_oneof![Just(CpgSide::C), Just(CpgSide::G)];
    let origin = prop_oneof![Just(CpgOrigin::Original), Just(CpgOrigin::DeNovo)];

    (side, origin)
        .prop_flat_map(|(side, origin)| {
            let denovo_ref = if origin == CpgOrigin::DeNovo {
                denovo_ref_strategy(side).boxed()
            } else {
                Just(A).boxed() // placeholder, unused for Original
            };

            // For DeNovo, unmod_reads >= 1 (CpG-forming alt must exist as reads).
            let unmod_lo = if origin == CpgOrigin::DeNovo { 1 } else { 0 };

            (
                Just(side),
                Just(origin),
                denovo_ref,
                read_count_strategy(0, 50),
                read_count_strategy(unmod_lo, 50),
                read_count_strategy(0, 10), // noise: wrong strand
                read_count_strategy(0, 10), // noise: wrong adj
            )
        })
        .prop_flat_map(|(side, origin, denovo_ref, mod_reads, unmod_reads, nws, nwa)| {
            let confounding = side.mod_base();
            let ref_base = match origin {
                CpgOrigin::Original => side.unmod_base(),
                CpgOrigin::DeNovo => denovo_ref,
            };

            let mut genotypes = vec![GenotypeScenario::None, GenotypeScenario::HomRef];

            // HetConfounded: for Original, need confounding base as alt (mod_reads >= 1).
            // For DeNovo where ref == confounding, any het is confounded (mod_reads can be 0
            // if unmod_reads provides the alt). For DeNovo where ref != confounding, need
            // confounding base as alt (mod_reads >= 1).
            let can_confound = if origin == CpgOrigin::DeNovo && ref_base == confounding {
                true // ref IS the confounding base
            } else {
                mod_reads >= 1 // confounding base exists as alt via mod reads
            };
            if can_confound {
                genotypes.push(GenotypeScenario::HetConfounded);
            }

            // HetNonConfounding: need a non-confounding alt, AND for DeNovo
            // the ref must not be the confounding base (otherwise any het is
            // confounded via the ref_base() == confounding check).
            let ref_is_confounding = origin == CpgOrigin::DeNovo && ref_base == confounding;
            if unmod_reads >= 1 && !ref_is_confounding {
                genotypes.push(GenotypeScenario::HetNonConfounding);
            }

            // HomAlt: only for Original, need any alt.
            if origin == CpgOrigin::Original && (mod_reads >= 1 || unmod_reads >= 1) {
                genotypes.push(GenotypeScenario::HomAlt);
            }
            // DeNovo + HomAlt: both chromosomes carry the variant → normal beta.
            // This is a distinct code path worth testing.
            if origin == CpgOrigin::DeNovo {
                genotypes.push(GenotypeScenario::HomAlt);
            }

            (
                Just(side),
                Just(origin),
                Just(denovo_ref),
                Just(mod_reads),
                Just(unmod_reads),
                Just(NoiseReads { wrong_strand: nws, wrong_adj: nwa }),
                proptest::sample::select(genotypes),
            )
        })
        .prop_map(|(side, origin, denovo_ref, mod_reads, unmod_reads, noise, genotype)| {
            MethylationScenario {
                side,
                origin,
                denovo_ref,
                mod_reads,
                unmod_reads,
                noise,
                genotype,
            }
        })
}

fn scenario_with_evidence() -> impl Strategy<Value = MethylationScenario> {
    scenario_strategy().prop_filter("need at least one read", |s| s.mod_reads + s.unmod_reads > 0)
}

// ---------------------------------------------------------------------------
// Property tests
// ---------------------------------------------------------------------------

proptest! {
    #![proptest_config(ProptestConfig::with_cases(PROPTEST_CASES))]

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
                    "Expected no evidence for {:?}, got {:?}", scenario, result
                );
            }
            Some(expected) => {
                let beta = extract_beta(&result, scenario.origin);
                let actual = match beta {
                    Some(b) => b,
                    None => {
                        prop_assert!(false,
                            "Expected beta={} for {:?}, got {:?}", expected, scenario, result
                        );
                        unreachable!()
                    }
                };
                prop_assert!(
                    (actual - expected).abs() < 1e-9,
                    "Beta mismatch for {:?}: expected {}, got {}", scenario, expected, actual
                );
            }
        }
    }

    /// Beta is always a valid probability [0, 1].
    #[test]
    fn beta_always_valid_probability(scenario in scenario_with_evidence()) {
        let metrics = scenario.build();
        let result = call(&metrics).unwrap();

        if let Some(methylated) = &result {
            for cpg in methylated.iter() {
                let b = *cpg.beta;
                prop_assert!(
                    (0.0..=1.0).contains(&b),
                    "Beta {} out of [0,1] for {:?}", b, scenario
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

        let beta_a = extract_beta(&result_a, scenario.origin);
        let beta_b = extract_beta(&result_b, mirror.origin);

        match (beta_a, beta_b) {
            (Some(a), Some(b)) => {
                prop_assert!(
                    (a - b).abs() < 1e-9,
                    "Symmetry broken: {:?} → {}, mirror {:?} → {}", scenario, a, mirror, b
                );
            }
            (None, None) => {} // both None is fine
            _ => {
                prop_assert!(false,
                    "Symmetry broken (one None): {:?} → {:?}, mirror {:?} → {:?}",
                    scenario, beta_a, mirror, beta_b
                );
            }
        }
    }

    /// Het-confounded adjustment always produces a beta <= the unadjusted beta.
    #[test]
    fn het_confounded_reduces_beta(
        side in prop_oneof![Just(CpgSide::C), Just(CpgSide::G)],
        origin in prop_oneof![Just(CpgOrigin::Original), Just(CpgOrigin::DeNovo)],
        mod_reads in 1u8..50,
        unmod_reads in 1u8..50,
    ) {
        let unadjusted = MethylationScenario {
            side, origin,
            denovo_ref: if origin == CpgOrigin::DeNovo { side.mod_base() } else { A },
            mod_reads, unmod_reads,
            noise: default(),
            genotype: GenotypeScenario::None,
        };
        let confounded = MethylationScenario {
            genotype: GenotypeScenario::HetConfounded,
            ..unadjusted.clone()
        };

        let beta_plain = extract_beta(&call(&unadjusted.build()).unwrap(), origin);
        let beta_het = extract_beta(&call(&confounded.build()).unwrap(), origin);

        if let (Some(plain), Some(het)) = (beta_plain, beta_het) {
            prop_assert!(
                het <= plain + 1e-9,
                "HetConfounded beta ({}) > unadjusted ({}) for {:?}", het, plain, confounded
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
            denovo_ref: A,
            mod_reads,
            unmod_reads,
            noise: default(),
            genotype: GenotypeScenario::HomAlt,
        };

        let beta = extract_beta(&call(&scenario.build()).unwrap(), CpgOrigin::Original);
        let beta = beta.expect("HomAlt with reads should produce result");
        prop_assert!(
            beta < 1e-9,
            "HomAlt Original should be 0.0, got {} for {:?}", beta, scenario
        );
    }

    /// HomAlt on a *DeNovo* CpG should give normal beta (both chromosomes
    /// carry the variant → CpG exists on both → no adjustment).
    #[test]
    fn hom_alt_denovo_is_normal_beta(
        side in prop_oneof![Just(CpgSide::C), Just(CpgSide::G)],
        mod_reads in 1u8..50,
        unmod_reads in 1u8..50,
    ) {
        let normal = MethylationScenario {
            side,
            origin: CpgOrigin::DeNovo,
            denovo_ref: side.mod_base(),
            mod_reads,
            unmod_reads,
            noise: default(),
            genotype: GenotypeScenario::None,
        };
        let hom_alt = MethylationScenario {
            genotype: GenotypeScenario::HomAlt,
            ..normal.clone()
        };

        let beta_normal = extract_beta(&call(&normal.build()).unwrap(), CpgOrigin::DeNovo);
        let beta_hom = extract_beta(&call(&hom_alt.build()).unwrap(), CpgOrigin::DeNovo);

        if let (Some(n), Some(h)) = (beta_normal, beta_hom) {
            prop_assert!(
                (n - h).abs() < 1e-9,
                "DeNovo HomAlt should match normal: {} vs {} for {:?}", n, h, hom_alt
            );
        }
    }

    /// More methylated reads (with fixed unmethylated count) should give
    /// a higher or equal beta.
    #[test]
    fn monotonic_in_mod_reads(
        side in prop_oneof![Just(CpgSide::C), Just(CpgSide::G)],
        origin in prop_oneof![Just(CpgOrigin::Original), Just(CpgOrigin::DeNovo)],
        mod_lo in 0u8..25,
        mod_extra in 1u8..25,
        unmod_reads in 1u8..50,
    ) {
        let mod_hi = mod_lo + mod_extra;

        let lo = MethylationScenario {
            side, origin,
            denovo_ref: side.mod_base(),
            mod_reads: mod_lo, unmod_reads,
            noise: default(),
            genotype: GenotypeScenario::None,
        };
        let hi = MethylationScenario { mod_reads: mod_hi, ..lo.clone() };

        let beta_lo = extract_beta(&call(&lo.build()).unwrap(), origin).unwrap_or(0.0);
        let beta_hi = extract_beta(&call(&hi.build()).unwrap(), origin).unwrap_or(0.0);

        prop_assert!(
            beta_hi >= beta_lo - 1e-9,
            "More mod reads should give higher beta: {}→{}, {}→{}", mod_lo, beta_lo, mod_hi, beta_hi
        );
    }

    /// mod_count and total_count in the result should match the raw read counts
    /// (noise reads must NOT be included).
    #[test]
    fn counts_match_reads(scenario in scenario_with_evidence()) {
        let metrics = scenario.build();
        let result = call(&metrics).unwrap();

        if let Some(methylated) = &result {
            let cpg = methylated.iter().find(|b| b.origin == scenario.origin);
            // This MUST find the CpG — if it doesn't, the build is broken.
            prop_assert!(
                cpg.is_some(),
                "Expected {:?} CpG in result {:?} for {:?}",
                scenario.origin, methylated, scenario
            );
            let cpg = cpg.unwrap();
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

    /// Noise reads (wrong strand, wrong adjacent base) must not affect the
    /// beta value. Compare a scenario with noise to the same scenario without.
    #[test]
    fn noise_reads_are_ignored(
        side in prop_oneof![Just(CpgSide::C), Just(CpgSide::G)],
        origin in prop_oneof![Just(CpgOrigin::Original), Just(CpgOrigin::DeNovo)],
        mod_reads in 1u8..30,
        unmod_reads in 1u8..30,
        wrong_strand in 1u8..20,
        wrong_adj in 1u8..20,
    ) {
        let clean = MethylationScenario {
            side, origin,
            denovo_ref: side.mod_base(),
            mod_reads, unmod_reads,
            noise: NoiseReads::default(),
            genotype: GenotypeScenario::None,
        };
        let noisy = MethylationScenario {
            noise: NoiseReads { wrong_strand, wrong_adj },
            ..clean.clone()
        };

        let beta_clean = extract_beta(&call(&clean.build()).unwrap(), origin);
        let beta_noisy = extract_beta(&call(&noisy.build()).unwrap(), origin);

        match (beta_clean, beta_noisy) {
            (Some(c), Some(n)) => {
                prop_assert!(
                    (c - n).abs() < 1e-9,
                    "Noise affected beta: clean={}, noisy={} for {:?}", c, n, noisy
                );
            }
            (None, None) => {}
            _ => {
                prop_assert!(false,
                    "Noise changed presence: clean={:?}, noisy={:?} for {:?}",
                    beta_clean, beta_noisy, noisy
                );
            }
        }
    }

    /// HomRef and HetNonConfounding should produce the same beta as no genotype.
    #[test]
    fn non_confounding_genotypes_are_transparent(
        side in prop_oneof![Just(CpgSide::C), Just(CpgSide::G)],
        origin in prop_oneof![Just(CpgOrigin::Original), Just(CpgOrigin::DeNovo)],
        mod_reads in 1u8..30,
        unmod_reads in 1u8..30,
        genotype in prop_oneof![
            Just(GenotypeScenario::HomRef),
            Just(GenotypeScenario::HetNonConfounding),
        ],
    ) {
        // For DeNovo + HetNonConfounding, ref must NOT be the confounding base.
        // Use a non-confounding ref (the CpG-forming base's adjacent counterpart).
        let denovo_ref = match side {
            CpgSide::C => A, // not T (confounding)
            CpgSide::G => T, // not A (confounding)
        };
        // Skip impossible combo: DeNovo + HetNonConfounding where ref==confounding
        // (already handled by using non-confounding denovo_ref above)

        let baseline = MethylationScenario {
            side, origin,
            denovo_ref,
            mod_reads, unmod_reads,
            noise: default(),
            genotype: GenotypeScenario::None,
        };
        let with_gt = MethylationScenario { genotype, ..baseline.clone() };

        let beta_none = extract_beta(&call(&baseline.build()).unwrap(), origin);
        let beta_gt = extract_beta(&call(&with_gt.build()).unwrap(), origin);

        match (beta_none, beta_gt) {
            (Some(a), Some(b)) => {
                prop_assert!(
                    (a - b).abs() < 1e-9,
                    "{:?} changed beta: {} vs {} for {:?}", genotype, a, b, with_gt
                );
            }
            (None, None) => {}
            _ => {
                prop_assert!(false,
                    "{:?} changed presence: {:?} vs {:?}", genotype, beta_none, beta_gt
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Negative tests: verify call() returns None when it should
// ---------------------------------------------------------------------------

proptest! {
    #![proptest_config(ProptestConfig::with_cases(PROPTEST_CASES))]

    /// Non-CpG context should never produce methylation calls.
    #[test]
    fn non_cpg_context_returns_none(
        ref_base in prop_oneof![Just(A), Just(T)],
        read_base in prop_oneof![Just(A), Just(T), Just(C), Just(G)],
        strand in prop_oneof![Just(Strand::OT), Just(Strand::OB)],
        n_reads in 1u8..20,
    ) {
        let reads: Vec<SimpleRead> = (0..n_reads)
            .map(|_| SimpleRead { base: read_base, strand, ..default() })
            .collect();

        let pileup = Pileup {
            region: dummy_region(),
            // No CpG context: ref is A or T, no adjacent C or G
            context: SequenceContext { me: ref_base, after_1: Some(A), before_1: Some(T), ..default() },
            pos: 1000,
            reads: SimpleReads(reads.into()),
            reference_base: ref_base,
            ..empty_pileup()
        };
        let metrics = PileupMetrics::new(pileup).unwrap();
        let result = call(&metrics).unwrap();

        prop_assert!(
            result.is_none(),
            "Non-CpG context produced methylation: {:?}", result
        );
    }

    /// All reads on the wrong strand should produce no evidence.
    #[test]
    fn wrong_strand_only_returns_no_evidence(
        side in prop_oneof![Just(CpgSide::C), Just(CpgSide::G)],
        n_reads in 1u8..20,
    ) {
        let (ref_base, context) = match side {
            CpgSide::C => (C, ctx_after_g(C)),
            CpgSide::G => (G, ctx_before_c(G)),
        };
        let wrong_strand = match side {
            CpgSide::C => Strand::OB,
            CpgSide::G => Strand::OT,
        };
        let mod_base = side.mod_base();
        let (bb, ab) = match side {
            CpgSide::C => (Option::None, Some(G)),
            CpgSide::G => (Some(C), Option::None),
        };

        let reads: Vec<SimpleRead> = (0..n_reads)
            .map(|_| SimpleRead {
                base: mod_base,
                strand: wrong_strand,
                before_base: bb,
                after_base: ab,
                ..default()
            })
            .collect();

        let pileup = Pileup {
            region: dummy_region(),
            context,
            pos: 1000,
            reads: SimpleReads(reads.into()),
            reference_base: ref_base,
            ..empty_pileup()
        };
        let metrics = PileupMetrics::new(pileup).unwrap();
        let result = call(&metrics).unwrap();

        prop_assert!(
            result.as_ref().is_none_or(|m| !m.has_evidence()),
            "Wrong-strand-only produced evidence: {:?}", result
        );
    }

    /// DeNovo alt that is NOT marked as RealVariant should not be detected.
    #[test]
    fn uncalled_denovo_returns_none(
        side in prop_oneof![Just(CpgSide::C), Just(CpgSide::G)],
        n_mod in 1u8..10,
        n_unmod in 1u8..10,
    ) {
        let denovo_ref = side.mod_base();
        let scenario = MethylationScenario {
            side,
            origin: CpgOrigin::DeNovo,
            denovo_ref,
            mod_reads: n_mod,
            unmod_reads: n_unmod,
            noise: default(),
            genotype: GenotypeScenario::None,
        };

        let mut metrics = scenario.build();
        // UN-mark the denovo alt — revert to Uncalled
        let denovo_base = side.unmod_base();
        if let Some(alt) = metrics.alts.iter_mut().find(|a| a.base == denovo_base) {
            alt.call = AltCall::Uncalled;
        }

        let result = call(&metrics).unwrap();
        let denovo_beta = extract_beta(&result, CpgOrigin::DeNovo);
        prop_assert!(
            denovo_beta.is_none(),
            "Uncalled denovo alt produced DeNovo beta: {:?}", result
        );
    }
}

// ---------------------------------------------------------------------------
// Targeted tests for specific code paths
// ---------------------------------------------------------------------------

/// `DenovoAdjecent` path: `cpg_origin` detects Original CpG via `denovo_adj` flag
/// (partner of a denovo CpG) rather than via `InCpG` from ref context.
#[test]
fn denovo_adjacent_produces_original_cpg() {
    for side in [CpgSide::C, CpgSide::G] {
        // Use a ref base that does NOT form a natural CpG context.
        let ref_base = side.unmod_base(); // C for C-side, G for G-side
        // Deliberately set the wrong adjacent base so InCpG::from would say No.
        let context = match side {
            CpgSide::C => SequenceContext { me: ref_base, after_1: Some(A), ..default() },
            CpgSide::G => SequenceContext { me: ref_base, before_1: Some(T), ..default() },
        };

        let strand = side.strand();
        let (bb, ab) = match side {
            // Reads still need the correct adj for PairedCounts
            CpgSide::C => (Option::None, Some(G)),
            CpgSide::G => (Some(C), Option::None),
        };
        let reads: Vec<SimpleRead> = (0..5)
            .map(|_| SimpleRead {
                base: side.mod_base(),
                strand,
                before_base: bb,
                after_base: ab,
                ..default()
            })
            .chain((0..5).map(|_| SimpleRead {
                base: side.unmod_base(),
                strand,
                before_base: bb,
                after_base: ab,
                ..default()
            }))
            .collect();

        let pileup = Pileup {
            region: dummy_region(),
            context,
            pos: 1000,
            reads: SimpleReads(reads.into()),
            reference_base: ref_base,
            ..empty_pileup()
        };
        let mut metrics = PileupMetrics::new(pileup).unwrap();

        // Set the denovo_adj flag — this is what makes cpg_origin detect it.
        metrics.pos_metrics.extended.denovo_adj = match side {
            CpgSide::C => DenovoAdjecent::ThisIsTheMatchingC,
            CpgSide::G => DenovoAdjecent::ThisIsTheMatchingG,
        };

        let result = call(&metrics).unwrap();
        let beta = extract_beta(&result, CpgOrigin::Original);
        assert!(
            beta.is_some(),
            "DenovoAdjecent should produce Original CpG for {:?}, got {:?}",
            side,
            result
        );
        let beta = beta.unwrap();
        assert!((beta - 0.5).abs() < 1e-9, "Expected beta=0.5, got {} for {:?}", beta, side);
    }
}

/// Both C-side and G-side can produce independent betas at the same position
/// (e.g., a position that is both the C of one CpG and the G of another).
#[test]
fn both_sides_independent() {
    // ref=C, before_1=C (makes InCpG::G from neighbor), after_1=G (makes InCpG::C)
    // This position is the C-side of a CpG (C followed by G) AND
    // we set denovo_adj to also make it the G-side partner.
    let context = SequenceContext { me: C, before_1: Some(C), after_1: Some(G), ..default() };

    // C-side reads (OT strand, adj G): 3 mod (T), 2 unmod (C)
    let c_reads: Vec<SimpleRead> = (0..3)
        .map(|_| SimpleRead {
            base: T,
            strand: Strand::OT,
            before_base: Option::None,
            after_base: Some(G),
            ..default()
        })
        .chain((0..2).map(|_| SimpleRead {
            base: C,
            strand: Strand::OT,
            before_base: Option::None,
            after_base: Some(G),
            ..default()
        }))
        .collect();

    // G-side reads (OB strand, adj C): 1 mod (A), 4 unmod (G)
    // Note: these use before_base=C because G-side looks at before_counts
    let g_reads: Vec<SimpleRead> = std::iter::once(SimpleRead {
        base: A,
        strand: Strand::OB,
        before_base: Some(C),
        after_base: Option::None,
        ..default()
    })
    .chain((0..4).map(|_| SimpleRead {
        base: G,
        strand: Strand::OB,
        before_base: Some(C),
        after_base: Option::None,
        ..default()
    }))
    .collect();

    let mut all_reads = c_reads;
    all_reads.extend(g_reads);

    let pileup = Pileup {
        region: dummy_region(),
        context,
        pos: 1000,
        reads: SimpleReads(all_reads.into()),
        reference_base: C,
        ..empty_pileup()
    };

    let mut metrics = PileupMetrics::new(pileup).unwrap();
    // Make this also the G-side partner of a denovo CpG
    metrics.pos_metrics.extended.denovo_adj = DenovoAdjecent::ThisIsTheMatchingG;

    let result = call(&metrics).unwrap();
    let methylated = result.as_ref().expect("should have methylation");

    // Verify both sides produced output
    assert!(
        methylated.0.len() == 2,
        "Expected 2 CpG entries (C-side + G-side), got {}: {:?}",
        methylated.0.len(),
        methylated
    );

    // C-side: 3 mod / 5 total = 0.6
    let c_beta = *methylated.0.first().expect("C-side entry").beta;
    assert!((c_beta - 0.6).abs() < 1e-9, "C-side beta: expected 0.6, got {}", c_beta);
    // G-side: 1 mod / 5 total = 0.2
    let g_beta = *methylated.0.get(1).expect("G-side entry").beta;
    assert!((g_beta - 0.2).abs() < 1e-9, "G-side beta: expected 0.2, got {}", g_beta);
}

// ---------------------------------------------------------------------------
// Formula verification: domain-derived analytical values, NOT copied from code
// ---------------------------------------------------------------------------
//
// These test adjusted_beta against values derived from first principles:
// - In TAPS, methylated C reads as T, unmethylated C reads as C
// - At a het C/T site, ~half of T reads come from the SNP chromosome
// - HomAlt means the reference base is gone on both chromosomes

#[test]
fn formula_no_adjustment_basic_ratios() {
    // All methylated → beta = 1.0
    assert_eq!(adjusted_beta(10.0, 0.0, GenotypeAdjustment::None), 1.0);
    // All unmethylated → beta = 0.0
    assert_eq!(adjusted_beta(0.0, 10.0, GenotypeAdjustment::None), 0.0);
    // Equal → beta = 0.5
    assert_eq!(adjusted_beta(5.0, 5.0, GenotypeAdjustment::None), 0.5);
    // 3:1 ratio
    assert!((adjusted_beta(3.0, 1.0, GenotypeAdjustment::None) - 0.75).abs() < 1e-9);
}

#[test]
fn formula_hom_alt_always_zero() {
    // Regardless of read counts, HomAlt = 0.0
    assert_eq!(adjusted_beta(100.0, 0.0, GenotypeAdjustment::HomAlt), 0.0);
    assert_eq!(adjusted_beta(0.0, 100.0, GenotypeAdjustment::HomAlt), 0.0);
    assert_eq!(adjusted_beta(50.0, 50.0, GenotypeAdjustment::HomAlt), 0.0);
}

#[test]
fn formula_het_confounded_known_values() {
    // 50/50 split: all "mod" reads could come from the SNP → beta = 0
    assert_eq!(adjusted_beta(5.0, 5.0, GenotypeAdjustment::HetConfounded), 0.0);

    // 6 mod, 4 unmod (total=10): SNP contributes ~5, excess = 1
    // beta = 1 / (4 + 1) = 0.2
    let b = adjusted_beta(6.0, 4.0, GenotypeAdjustment::HetConfounded);
    assert!((b - 0.2).abs() < 1e-9, "got {b}");

    // All mod reads (10 mod, 0 unmod): SNP contributes ~5, excess = 5
    // beta = 5 / (0 + 5) = 1.0
    let b = adjusted_beta(10.0, 0.0, GenotypeAdjustment::HetConfounded);
    assert!((b - 1.0).abs() < 1e-9, "got {b}");

    // 7 mod, 3 unmod (total=10): SNP contributes ~5, excess = 2
    // beta = 2 / (3 + 2) = 0.4
    let b = adjusted_beta(7.0, 3.0, GenotypeAdjustment::HetConfounded);
    assert!((b - 0.4).abs() < 1e-9, "got {b}");

    // Fewer mod than half (3 mod, 7 unmod): excess = max(3-5, 0) = 0
    // beta = 0 / (7 + 0) = 0
    let b = adjusted_beta(3.0, 7.0, GenotypeAdjustment::HetConfounded);
    assert!(b < 1e-9, "got {b}");
}

#[test]
fn formula_het_confounded_never_exceeds_unadjusted() {
    for m in 0..=20 {
        for u in 1..=20 {
            let plain = adjusted_beta(f64::from(m), f64::from(u), GenotypeAdjustment::None);
            let het = adjusted_beta(f64::from(m), f64::from(u), GenotypeAdjustment::HetConfounded);
            assert!(het <= plain + 1e-9, "het ({het}) > plain ({plain}) for m={m}, u={u}");
        }
    }
}
