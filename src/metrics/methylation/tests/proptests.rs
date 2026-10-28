use std::num::NonZeroU8;

use super::*;
use crate::{
    call::{
        tests::pileup::{Pileup, SimpleRead, SimpleReads},
        variant_calling::{EstimatedGenotype, GenotypeSource, GenotypeTag},
    },
    metrics::{AlleleMetrics, Alt},
    sequence::{ChunkRegion, Region, Segment},
    utils::default,
    vcf::{CpgOrigin, SequenceContext},
};
use hegel::{TestCase, generators as gs};
use seqair_types::{Base, Phred, Strand};

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

impl MethylationScenario {
    /// What `call()` should return for this scenario.
    ///
    /// NOTE: This oracle mirrors the production formula. The formula itself is
    /// verified against domain-derived analytical values in the `formula_*`
    /// tests at the bottom of this file. This property verifies the *wiring*:
    /// that the right reads reach the right formula with the right adjustment.
    fn expected_beta(&self) -> Option<f64> {
        let m = f64::from(self.mod_reads);
        let u = f64::from(self.unmod_reads);
        if m + u == 0.0 {
            return None;
        }
        Some(match self.effective_adjustment() {
            GenotypeAdjustment::NoReference => 0.0,
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
                GenotypeAdjustment::NoReference
            }
            GenotypeScenario::HetConfounded => GenotypeAdjustment::HetConfounded,
            _ => GenotypeAdjustment::None,
        }
    }

    /// Construct a `PileupMetrics` that embodies this scenario.
    fn build(&self) -> PileupMetrics {
        let context = self.context();

        let strand = self.side.strand();
        let wrong_strand = match strand {
            Strand::OT => Strand::OB,
            Strand::OB => Strand::OT,
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

        let mut metrics = column(context, reads);

        // For de-novo: mark the new CpG-forming alt as a real variant.
        if self.origin == CpgOrigin::DeNovo {
            let denovo_base = self.side.unmod_base();
            if let Some(alt) = metrics.alts.iter_mut().find(|a| a.base == denovo_base) {
                alt.verdict = crate::metrics::AltVerdict::Called(seqair_types::Probability::ONE);
            }
        }

        let gt = self.build_genotype(&mut metrics);
        metrics.pos_metrics.genotype = gt;
        metrics
    }

    fn context(&self) -> SequenceContext {
        match (self.side, self.origin) {
            (CpgSide::C, CpgOrigin::Original) => ctx_after_g(C),
            (CpgSide::G, CpgOrigin::Original) => ctx_before_c(G),
            (CpgSide::C, CpgOrigin::DeNovo) => ctx_after_g(self.denovo_ref),
            (CpgSide::G, CpgOrigin::DeNovo) => ctx_before_c(self.denovo_ref),
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
                // Point the het at the confounding base where it is an alt.
                // For denovo where ref == confounding base it is not, and any
                // het triggers confounding via the ref_base() check, so the
                // first alt will do. Alts are ordered by depth, not by the
                // order the reads were pushed, so the index has to be looked up.
                let confounding = self.side.mod_base();
                let idx = metrics.alts.iter().position(|a| a.base == confounding).unwrap_or(0);
                GenotypeTag::RefHet(NonZeroU8::new(idx as u8 + 1).unwrap())
            }
            GenotypeScenario::HetNonConfounding => self.build_het_non_confounding(metrics),
            GenotypeScenario::HomAlt => GenotypeTag::HomAlt(NonZeroU8::new(1).unwrap()),
        };
        Some(EstimatedGenotype {
            genotype: gt,
            quality: Phred::from_phred(20),
            source: GenotypeSource::Binomial,
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

    /// Mirror this scenario to the opposite CpG side.
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

/// The column at the `me` of `context`, on a segment holding exactly the
/// context's known bases.
fn column(context: SequenceContext, reads: Vec<SimpleRead>) -> PileupMetrics {
    let before: Vec<Base> = [context.before_2, context.before_1].into_iter().flatten().collect();
    let after = [context.after_1, context.after_2].into_iter().flatten();
    let sequence: Vec<u8> =
        before.iter().copied().chain([context.me]).chain(after).map(|b| b as u8).collect();
    let start = 1000;
    let end = start + sequence.len() as u64;
    let segment = Segment {
        range: std::sync::Arc::new(ChunkRegion {
            region: Region { contig: "chr_test".into(), start, end },
            last_position: end,
            overlap_start: 0,
            overlap_end: 0,
        }),
        sequence,
    };
    let pileup = Pileup {
        region: segment.range.clone(),
        pos: u32::try_from(start + before.len() as u64).unwrap(),
        reads: SimpleReads(reads.into()),
        reference_base: context.me,
        context,
    };
    pileup.metrics(&segment).unwrap()
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
        verdict: crate::metrics::AltVerdict::Unscored,
    });
}

fn extract_beta(result: &Methylated, origin: CpgOrigin) -> Option<f64> {
    result.get(origin).map(|b| *b.beta)
}

// ---------------------------------------------------------------------------
// Strategy
// ---------------------------------------------------------------------------

fn sides() -> impl gs::PrintableGenerator<CpgSide> {
    gs::sampled_from([CpgSide::C, CpgSide::G].as_slice())
}

fn origins() -> impl gs::PrintableGenerator<CpgOrigin> {
    gs::sampled_from([CpgOrigin::Original, CpgOrigin::DeNovo].as_slice())
}

fn read_counts(lo: u8, below: u8) -> impl gs::PrintableGenerator<u8> {
    gs::integers::<u8>().min_value(lo).max_value(below - 1)
}

/// A scenario; with `evidence`, one with at least one informative read,
/// drawn that way rather than filtered.
#[hegel::composite]
fn scenarios(tc: &TestCase, evidence: bool) -> MethylationScenario {
    let side = tc.draw(sides());
    let origin = tc.draw(origins());
    // Valid denovo ref bases: anything except the CpG base itself. Unused for
    // Original.
    let denovo_ref = match (origin, side) {
        (CpgOrigin::Original, _) => A,
        // C-side: ref was T, A, or G — mutated to C
        (CpgOrigin::DeNovo, CpgSide::C) => tc.draw(gs::sampled_from([T, A, G].as_slice())),
        // G-side: ref was A, T, or C — mutated to G
        (CpgOrigin::DeNovo, CpgSide::G) => tc.draw(gs::sampled_from([A, T, C].as_slice())),
    };
    let mod_reads = tc.draw(read_counts(0, 50));
    // For DeNovo, unmod_reads >= 1 (CpG-forming alt must exist as reads).
    let unmod_lo = u8::from(origin == CpgOrigin::DeNovo || (evidence && mod_reads == 0));
    let unmod_reads = tc.draw(read_counts(unmod_lo, 50));
    let noise = NoiseReads {
        wrong_strand: tc.draw(read_counts(0, 10)),
        wrong_adj: tc.draw(read_counts(0, 10)),
    };

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
    let genotype = tc.draw(gs::sampled_from(genotypes));

    MethylationScenario { side, origin, denovo_ref, mod_reads, unmod_reads, noise, genotype }
}

// ---------------------------------------------------------------------------
// Property tests
// ---------------------------------------------------------------------------

/// The core roundtrip: build reads from a scenario, call methylation,
/// and verify the result matches the expected beta.
#[hegel::test(test_cases = crate::pinned::cases(2048))]
fn roundtrip(tc: TestCase) {
    let scenario = tc.draw(scenarios(false));
    let metrics = scenario.build();
    let result = call(&metrics).unwrap();

    match scenario.expected_beta() {
        None => {
            assert!(result.is_empty(), "Expected no evidence for {:?}, got {:?}", scenario, result);
        }
        Some(expected) => {
            let beta = extract_beta(&result, scenario.origin);
            let actual = match beta {
                Some(b) => b,
                None => {
                    panic!("Expected beta={} for {:?}, got {:?}", expected, scenario, result)
                }
            };
            assert!(
                (actual - expected).abs() < 1e-9,
                "Beta mismatch for {:?}: expected {}, got {}",
                scenario,
                expected,
                actual
            );
        }
    }
}

/// Beta is always a valid probability [0, 1].
#[hegel::test(test_cases = crate::pinned::cases(2048))]
fn beta_always_valid_probability(tc: TestCase) {
    let scenario = tc.draw(scenarios(true));
    let metrics = scenario.build();
    let result = call(&metrics).unwrap();

    for cpg in [result.original, result.denovo].iter().flatten() {
        let b = *cpg.beta;
        assert!((0.0..=1.0).contains(&b), "Beta {} out of [0,1] for {:?}", b, scenario);
    }
}

/// `M5mC_Strands` and the BED's `unmod` and `mod` columns count the reads β
/// is computed from, so β can be recomputed from what is written. Reads on
/// the converting strand whose neighbour is not the CpG partner (`noise`'s
/// `wrong_adj`) are not in β.
#[hegel::test]
fn written_counts_are_the_reads_of_beta(tc: TestCase) {
    let scenario = tc.draw(scenarios(true));
    let metrics = scenario.build();
    let cpg = *call(&metrics).unwrap().get(scenario.origin).expect("evidence");
    let written = crate::metrics::MethylationEvidenceStrandInfo::from_pileup(&metrics);
    assert_eq!(
        (written.modified, written.unmod + written.modified),
        (cpg.mod_count, cpg.total_count.get()),
        "{scenario:?}"
    );
}

/// C-side and G-side are symmetric: mirroring a scenario should produce
/// the same beta value.
#[hegel::test(test_cases = crate::pinned::cases(2048))]
fn c_g_symmetry(tc: TestCase) {
    let scenario = tc.draw(scenarios(true));
    let mirror = scenario.mirror();

    let result_a = call(&scenario.build()).unwrap();
    let result_b = call(&mirror.build()).unwrap();

    let beta_a = extract_beta(&result_a, scenario.origin);
    let beta_b = extract_beta(&result_b, mirror.origin);

    match (beta_a, beta_b) {
        (Some(a), Some(b)) => {
            assert!(
                (a - b).abs() < 1e-9,
                "Symmetry broken: {:?} → {}, mirror {:?} → {}",
                scenario,
                a,
                mirror,
                b
            );
        }
        (None, None) => {} // both None is fine
        _ => {
            panic!(
                "Symmetry broken (one None): {:?} → {:?}, mirror {:?} → {:?}",
                scenario, beta_a, mirror, beta_b
            );
        }
    }
}

/// Het-confounded adjustment always produces a beta <= the unadjusted beta.
#[hegel::test(test_cases = crate::pinned::cases(2048))]
fn het_confounded_reduces_beta(tc: TestCase) {
    let side = tc.draw(sides());
    let origin = tc.draw(origins());
    let mod_reads = tc.draw(read_counts(1, 50));
    let unmod_reads = tc.draw(read_counts(1, 50));
    let unadjusted = MethylationScenario {
        side,
        origin,
        denovo_ref: if origin == CpgOrigin::DeNovo { side.mod_base() } else { A },
        mod_reads,
        unmod_reads,
        noise: default(),
        genotype: GenotypeScenario::None,
    };
    let confounded =
        MethylationScenario { genotype: GenotypeScenario::HetConfounded, ..unadjusted.clone() };

    let beta_plain = extract_beta(&call(&unadjusted.build()).unwrap(), origin);
    let beta_het = extract_beta(&call(&confounded.build()).unwrap(), origin);

    if let (Some(plain), Some(het)) = (beta_plain, beta_het) {
        assert!(
            het <= plain + 1e-9,
            "HetConfounded beta ({}) > unadjusted ({}) for {:?}",
            het,
            plain,
            confounded
        );
    }
}

/// `HomAlt` on an original CpG always produces beta = 0.0.
#[hegel::test(test_cases = crate::pinned::cases(2048))]
fn hom_alt_original_is_zero(tc: TestCase) {
    let side = tc.draw(sides());
    let mod_reads = tc.draw(read_counts(1, 50));
    let unmod_reads = tc.draw(read_counts(0, 50));
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
    assert!(beta < 1e-9, "HomAlt Original should be 0.0, got {} for {:?}", beta, scenario);
}

/// `HomAlt` on a `DeNovo` CpG should give normal beta (both chromosomes
/// carry the variant → CpG exists on both → no adjustment).
#[hegel::test(test_cases = crate::pinned::cases(2048))]
fn hom_alt_denovo_is_normal_beta(tc: TestCase) {
    let side = tc.draw(sides());
    let mod_reads = tc.draw(read_counts(1, 50));
    let unmod_reads = tc.draw(read_counts(1, 50));
    let normal = MethylationScenario {
        side,
        origin: CpgOrigin::DeNovo,
        denovo_ref: side.mod_base(),
        mod_reads,
        unmod_reads,
        noise: default(),
        genotype: GenotypeScenario::None,
    };
    let hom_alt = MethylationScenario { genotype: GenotypeScenario::HomAlt, ..normal.clone() };

    let beta_normal = extract_beta(&call(&normal.build()).unwrap(), CpgOrigin::DeNovo);
    let beta_hom = extract_beta(&call(&hom_alt.build()).unwrap(), CpgOrigin::DeNovo);

    if let (Some(n), Some(h)) = (beta_normal, beta_hom) {
        assert!(
            (n - h).abs() < 1e-9,
            "DeNovo HomAlt should match normal: {} vs {} for {:?}",
            n,
            h,
            hom_alt
        );
    }
}

/// More methylated reads (with fixed unmethylated count) should give
/// a higher or equal beta.
#[hegel::test(test_cases = crate::pinned::cases(2048))]
fn monotonic_in_mod_reads(tc: TestCase) {
    let side = tc.draw(sides());
    let origin = tc.draw(origins());
    let mod_lo = tc.draw(read_counts(0, 25));
    let mod_extra = tc.draw(read_counts(1, 25));
    let unmod_reads = tc.draw(read_counts(1, 50));
    let mod_hi = mod_lo + mod_extra;

    let lo = MethylationScenario {
        side,
        origin,
        denovo_ref: side.mod_base(),
        mod_reads: mod_lo,
        unmod_reads,
        noise: default(),
        genotype: GenotypeScenario::None,
    };
    let hi = MethylationScenario { mod_reads: mod_hi, ..lo.clone() };

    let beta_lo = extract_beta(&call(&lo.build()).unwrap(), origin).unwrap_or(0.0);
    let beta_hi = extract_beta(&call(&hi.build()).unwrap(), origin).unwrap_or(0.0);

    assert!(
        beta_hi >= beta_lo - 1e-9,
        "More mod reads should give higher beta: {}→{}, {}→{}",
        mod_lo,
        beta_lo,
        mod_hi,
        beta_hi
    );
}

/// `mod_count` and `total_count` in the result should match the raw read counts
/// (noise reads must NOT be included).
#[hegel::test(test_cases = crate::pinned::cases(2048))]
fn counts_match_reads(tc: TestCase) {
    let scenario = tc.draw(scenarios(true));
    let metrics = scenario.build();
    let result = call(&metrics).unwrap();

    let cpg = result.get(scenario.origin);
    // This MUST find the CpG — if it doesn't, the build is broken.
    assert!(
        cpg.is_some(),
        "Expected {:?} CpG in result {:?} for {:?}",
        scenario.origin,
        result,
        scenario
    );
    let cpg = cpg.unwrap();
    assert_eq!(
        cpg.mod_count,
        u32::from(scenario.mod_reads),
        "mod_count mismatch for {:?}",
        scenario
    );
    assert_eq!(
        cpg.total_count.get(),
        u32::from(scenario.mod_reads) + u32::from(scenario.unmod_reads),
        "total_count mismatch for {:?}",
        scenario
    );
}

/// Noise reads (wrong strand, wrong adjacent base) must not affect the
/// beta value. Compare a scenario with noise to the same scenario without.
#[hegel::test(test_cases = crate::pinned::cases(2048))]
fn noise_reads_are_ignored(tc: TestCase) {
    let side = tc.draw(sides());
    let origin = tc.draw(origins());
    let mod_reads = tc.draw(read_counts(1, 30));
    let unmod_reads = tc.draw(read_counts(1, 30));
    let wrong_strand = tc.draw(read_counts(1, 20));
    let wrong_adj = tc.draw(read_counts(1, 20));
    let clean = MethylationScenario {
        side,
        origin,
        denovo_ref: side.mod_base(),
        mod_reads,
        unmod_reads,
        noise: NoiseReads::default(),
        genotype: GenotypeScenario::None,
    };
    let noisy =
        MethylationScenario { noise: NoiseReads { wrong_strand, wrong_adj }, ..clean.clone() };

    let beta_clean = extract_beta(&call(&clean.build()).unwrap(), origin);
    let beta_noisy = extract_beta(&call(&noisy.build()).unwrap(), origin);

    match (beta_clean, beta_noisy) {
        (Some(c), Some(n)) => {
            assert!(
                (c - n).abs() < 1e-9,
                "Noise affected beta: clean={}, noisy={} for {:?}",
                c,
                n,
                noisy
            );
        }
        (None, None) => {}
        _ => {
            panic!(
                "Noise changed presence: clean={:?}, noisy={:?} for {:?}",
                beta_clean, beta_noisy, noisy
            );
        }
    }
}

/// `HomRef` and `HetNonConfounding` should produce the same beta as no genotype.
#[hegel::test(test_cases = crate::pinned::cases(2048))]
fn non_confounding_genotypes_are_transparent(tc: TestCase) {
    let side = tc.draw(sides());
    let origin = tc.draw(origins());
    let mod_reads = tc.draw(read_counts(1, 30));
    let unmod_reads = tc.draw(read_counts(1, 30));
    let genotype = tc.draw(gs::sampled_from(
        [GenotypeScenario::HomRef, GenotypeScenario::HetNonConfounding].as_slice(),
    ));
    // For DeNovo + HetNonConfounding, ref must NOT be the confounding base.
    // Use a non-confounding ref (the CpG-forming base's adjacent counterpart).
    let denovo_ref = match side {
        CpgSide::C => A, // not T (confounding)
        CpgSide::G => T, // not A (confounding)
    };
    // Skip impossible combo: DeNovo + HetNonConfounding where ref==confounding
    // (already handled by using non-confounding denovo_ref above)

    let baseline = MethylationScenario {
        side,
        origin,
        denovo_ref,
        mod_reads,
        unmod_reads,
        noise: default(),
        genotype: GenotypeScenario::None,
    };
    let with_gt = MethylationScenario { genotype, ..baseline.clone() };

    let beta_none = extract_beta(&call(&baseline.build()).unwrap(), origin);
    let beta_gt = extract_beta(&call(&with_gt.build()).unwrap(), origin);

    match (beta_none, beta_gt) {
        (Some(a), Some(b)) => {
            assert!(
                (a - b).abs() < 1e-9,
                "{:?} changed beta: {} vs {} for {:?}",
                genotype,
                a,
                b,
                with_gt
            );
        }
        (None, None) => {}
        _ => {
            panic!("{:?} changed presence: {:?} vs {:?}", genotype, beta_none, beta_gt);
        }
    }
}

// ---------------------------------------------------------------------------
// Negative tests: verify call() returns None when it should
// ---------------------------------------------------------------------------

/// Non-CpG context should never produce methylation calls.
#[hegel::test(test_cases = crate::pinned::cases(2048))]
fn non_cpg_context_returns_none(tc: TestCase) {
    let ref_base = tc.draw(gs::sampled_from([A, T].as_slice()));
    let read_base = tc.draw(gs::sampled_from([A, T, C, G].as_slice()));
    let strand = tc.draw(gs::sampled_from([Strand::OT, Strand::OB].as_slice()));
    let n_reads = tc.draw(read_counts(1, 20));
    let reads: Vec<SimpleRead> =
        (0..n_reads).map(|_| SimpleRead { base: read_base, strand, ..default() }).collect();

    let metrics = column(
        SequenceContext { me: ref_base, after_1: Some(A), before_1: Some(T), ..default() },
        reads,
    );
    let result = call(&metrics).unwrap();

    assert!(result.is_empty(), "Non-CpG context produced methylation: {:?}", result);
}

/// All reads on the wrong strand should produce no evidence.
#[hegel::test(test_cases = crate::pinned::cases(2048))]
fn wrong_strand_only_returns_no_evidence(tc: TestCase) {
    let side = tc.draw(sides());
    let n_reads = tc.draw(read_counts(1, 20));
    let context = match side {
        CpgSide::C => ctx_after_g(C),
        CpgSide::G => ctx_before_c(G),
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

    let metrics = column(context, reads);
    let result = call(&metrics).unwrap();

    assert!(result.is_empty(), "Wrong-strand-only produced evidence: {:?}", result);
}

/// `DeNovo` alt that is NOT marked as `RealVariant` should not be detected.
#[hegel::test(test_cases = crate::pinned::cases(2048))]
fn uncalled_denovo_returns_none(tc: TestCase) {
    let side = tc.draw(sides());
    let n_mod = tc.draw(read_counts(1, 10));
    let n_unmod = tc.draw(read_counts(1, 10));
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
        alt.verdict = crate::metrics::AltVerdict::Unscored;
    }

    let result = call(&metrics).unwrap();
    let denovo_beta = extract_beta(&result, CpgOrigin::DeNovo);
    assert!(denovo_beta.is_none(), "Uncalled denovo alt produced DeNovo beta: {:?}", result);
}

// ---------------------------------------------------------------------------
// Targeted tests for specific code paths
// ---------------------------------------------------------------------------

#[test]
fn denovo_adjacent_without_called_partner_is_not_a_cpg() {
    for side in [CpgSide::C, CpgSide::G] {
        let ref_base = side.unmod_base();
        let context = match side {
            CpgSide::C => SequenceContext { me: ref_base, after_1: Some(A), ..default() },
            CpgSide::G => SequenceContext { me: ref_base, before_1: Some(T), ..default() },
        };

        let strand = side.strand();
        let (bb, ab) = match side {
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

        let metrics = column(context, reads);

        let result = call(&metrics).unwrap();
        assert!(
            result.is_empty(),
            "adjacency to a rejected denovo candidate produced methylation for {:?}: {:?}",
            side,
            result
        );
    }
}

/// The partner of a called de-novo alt: `cpg_origin` detects an Original CpG
/// beside the called alt rather than from the reference context.
#[test]
fn denovo_adjacent_produces_original_cpg() {
    for side in [CpgSide::C, CpgSide::G] {
        // Use a ref base that does NOT form a natural CpG context.
        let ref_base = side.unmod_base(); // C for C-side, G for G-side
        // Deliberately set the wrong adjacent base so the reference context is no CpG.
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

        let mut metrics = column(context, reads);

        metrics.beside_called_alts = enumset::EnumSet::only(side);

        let result = call(&metrics).unwrap();
        let beta = extract_beta(&result, CpgOrigin::Original);
        assert!(
            beta.is_some(),
            "a called de-novo partner should produce an Original CpG for {:?}, got {:?}",
            side,
            result
        );
        let beta = beta.unwrap();
        assert!((beta - 0.5).abs() < 1e-9, "Expected beta=0.5, got {} for {:?}", beta, side);
    }
}

/// Both C-side and G-side can produce independent betas at the same position:
/// the `C` of a reference CpG that a called `G` alt makes the `G` of another.
#[test]
fn both_sides_independent() {
    // ref=C, after_1=G (a reference CpG C); before_1=C, so a G here forms a CpG with it.
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

    let mut metrics = column(context, all_reads);
    if let Some(alt) = metrics.alts.iter_mut().find(|a| a.base == G) {
        alt.verdict = crate::metrics::AltVerdict::Called(seqair_types::Probability::ONE);
    }

    let methylated = call(&metrics).unwrap();

    // C-side: 3 mod / 5 total = 0.6
    let c_beta = *methylated.original.expect("C-side (reference CpG) entry").beta;
    assert!((c_beta - 0.6).abs() < 1e-9, "C-side beta: expected 0.6, got {}", c_beta);
    // G-side: 1 mod / 5 total = 0.2
    let g_beta = *methylated.denovo.expect("G-side (de-novo CpG) entry").beta;
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
    assert_eq!(adjusted_beta(100.0, 0.0, GenotypeAdjustment::NoReference), 0.0);
    assert_eq!(adjusted_beta(0.0, 100.0, GenotypeAdjustment::NoReference), 0.0);
    assert_eq!(adjusted_beta(50.0, 50.0, GenotypeAdjustment::NoReference), 0.0);
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
