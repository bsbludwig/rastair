//! Non-ML hard-filter indel calling.
//!
//! The pathway used by `rastair call --experimental-indels` (and by
//! `--experimental-indels-ml --no-ml`, which degrades to it): candidate indels are
//! accepted by a fixed chain of hard filters — minimum filtered depth, minimum
//! alternate observations, an OT/OB strand-bias test, and a binomial genotype
//! test — with no ML scoring.
//!
//! Notes:
//!   * Genotyping excludes noisy fragments — a terminal tandem repeat or a
//!     soft-clip — under `--indel-noise-exclusion`, whose three settings differ in
//!     whether a noisy *supporting* fragment counts toward the min-AO gate and
//!     toward the ratio's numerator; see
//!     [`IndelNoiseExclusion`](super::IndelNoiseExclusion). The default
//!     (`symmetric`) drops it from both, which leaves VAF unbiased. All of this is
//!     kept separate from the ML-facing `depth_offset`, which stays one-sided so
//!     ML feature distributions do not move. Every count involved is accumulated
//!     over exactly the reads that end up in `Pileup::reads`.
//!   * "Noisy" is decided without looking at the observed base, so — unlike the
//!     read-level mismatch filter, which is TAPS-aware — the exclusion cannot vary
//!     with local methylation.
//!   * The default strand gate is **both-strand concordance**
//!     (`IndelAlleleCounts::on_both_strands`): an allele confined to one strand
//!     fails. This is a *presence* criterion — genuine TAPS indels are frequently
//!     strand-asymmetric (alt-allele reference bias, not methylation), so a *degree*
//!     test rejects real calls, whereas single-strandedness is ~9x enriched in
//!     artefacts. Its known cost is the low-AO regime: a genuine 2/0 split is a coin
//!     flip yet fails here. Trade-off measured genotype-aware (hap.py) vs
//!     genotype-blind (`verify`); the two metrics rank concordance vs the
//!     significance test differently, which is a live discussion.
//!   * The significance test (`strand_bias_p_value`) stays available as an **opt-in
//!     additional** gate via `--indel-strand-bias-alpha` (default 0 = off, since a
//!     p-value is never < 0). Its null is the locus' own strand mix, taken over all
//!     supporting fragments, noisy included, so both sides are on the same footing.
//!     At the 0.05 it once defaulted to it rejected 4,288 true chr12 indels to
//!     remove 252 false ones; the p-value is still informative as an ML feature.
//!   * Strand is judged on OT/OB, not the alignment reverse flag — see
//!     [`crate::call::pileup::indels::IndelAlleleCounts`] for why the reverse flag
//!     cannot work under per-fragment deduplication.
//!   * Failing alleles are emitted with a FILTER tag rather than dropped, matching
//!     rastair's "always emit indels" convention. `IndelVerdict` travels on the call
//!     to the VCF layer, so `rastair convert` renders the same FILTERs as a direct
//!     VCF run.
//!   * Tie-break: `super::binomial_genotype` resolves exact-probability ties toward
//!     hom-ref (`>=`); such exact ties do not occur with continuous binomial masses.
//!
//! To remove the whole hard-filter pathway, delete this file and:
//!   * in the parent `indel_calling.rs`: `pub mod hard_filters;`, the
//!     `hard_filter_verdict` field on [`IndelCall`], the `!use_ml` delegation, and
//!     the `experimental_indels` / `experimental_indels_ml` split on [`IndelParams`]
//!   * the verdict-rendering arm in `build_indel_records` (`src/vcf/metrics_to_vcf.rs`)
//!   * the `indel_strand` / `indel_hom_ref` FILTER ids (`src/vcf.rs`)
//!   * `Pileup::noisy_ref_count` (`src/call/pileup.rs`), `IndelCounts::noisy_ref_count`
//!     / `clean_depth` and `IndelAlleleCounts::{noisy, clean_total}`
//!     (`src/call/pileup/indels.rs`), and `IndelObservation::noisy`
//!
//! The OT/OB split on `IndelAlleleCounts` and the fragment-level indel accounting in
//! `from_hts.rs` are *not* part of this pathway — the ML path depends on them too.

use super::{IndelCall, IndelParams};
use crate::call::pileup::indels::IndelCounts;
use crate::call::variant_calling::GenotypeTag;
use tracing::trace;

/// Verdict of the non-ML hard-filter indel chain.
///
/// Every allele that reaches genotyping is emitted; failures carry the reason so
/// `build_indel_records` can render the matching VCF FILTER tag.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum IndelVerdict {
    /// Passed every hard filter.
    Pass,
    /// OT/OB split of the supporting fragments is significantly skewed against
    /// the strand mix of the rest of the locus.
    FailStrand,
    /// Genotyped as homozygous reference by the binomial model.
    FailHomRef,
}

/// Call indels with the non-ML hard-filter chain.
///
/// Shares the locus depth gate, per-allele min-AO and binomial genotype of
/// [`super::call_indels`], and adds an OT/OB strand-bias test. Unlike the ML
/// path, every genotyped allele is emitted with an [`IndelVerdict`] instead of
/// being dropped.
pub fn call_indels(indels: &IndelCounts, params: &IndelParams) -> Vec<IndelCall> {
    let mut calls = Vec::new();

    // How noisy fragments (terminal tandem repeat or soft-clip) are kept out of
    // the ratio is `--indel-noise-exclusion`; see [`IndelNoiseExclusion`]. All of
    // it is kept separate from the ML-facing `depth_offset`, which stays
    // one-sided because the shipped models were trained against it.
    let noise = params.indel_noise_exclusion;
    let filtered_depth = noise.depth(indels);

    if filtered_depth < params.min_indel_depth {
        trace!(
            depth = filtered_depth,
            min = params.min_indel_depth,
            noisy_ref = indels.noisy_ref_count,
            "Indel locus skipped: below min depth"
        );
        return calls;
    }

    for allele_counts in &indels.alleles {
        let (observations, alt_count) = noise.alt_counts(allele_counts);
        if observations < params.min_indel_ao {
            trace!(
                allele = ?allele_counts.allele,
                observations,
                noisy = allele_counts.noisy,
                min = params.min_indel_ao,
                "Indel skipped: below min AO"
            );
            continue;
        }

        let Some(genotype) = super::binomial_genotype(
            alt_count as usize,
            filtered_depth as usize,
            params.indel_error_rate,
            params.indel_het_vaf,
        ) else {
            trace!(allele = ?allele_counts.allele, "Indel skipped: not genotypable");
            continue;
        };

        // Both-strand concordance is the default strand gate: an allele confined to
        // a single strand fails. The significance test (below) is judged on *all*
        // supporting fragments, noisy included, against the locus' own strand mix,
        // and stays available as an opt-in additional gate via
        // `--indel-strand-bias-alpha` (default 0 leaves it off, since a p-value is
        // never < 0). See the module docs for why presence beats degree on TAPS.
        // FailStrand if the allele is single-stranded (concordance, default) or, when
        // the opt-in significance test is enabled, if its OT/OB split is significantly
        // biased against the locus' own strand mix (alpha default 0 = off).
        let strand_bias = allele_counts.strand_bias_p_value(indels.null_ot_fraction(allele_counts));
        let verdict = if !allele_counts.on_both_strands()
            || strand_bias < params.indel_strand_bias_alpha
        {
            IndelVerdict::FailStrand
        } else if matches!(genotype.tag, GenotypeTag::HomRef) {
            IndelVerdict::FailHomRef
        } else {
            IndelVerdict::Pass
        };

        calls.push(IndelCall {
            allele: allele_counts.allele.clone(),
            genotype: genotype.tag,
            quality: genotype.quality,
            ml: None,
            depth: filtered_depth,
            alt_count,
            hard_filter_verdict: Some(verdict),
        });
    }

    calls
}

#[cfg(test)]
mod tests {
    use super::super::IndelNoiseExclusion;
    use super::*;
    use crate::call::pileup::indels::{IndelAllele, IndelAlleleCounts};
    use seqair_types::{Base, SmallVec};

    fn insertion(seq: &str, ot: u32, ob: u32) -> IndelAlleleCounts {
        let bases: SmallVec<Base, 4> = seq.bytes().map(Base::from).collect();
        IndelAlleleCounts {
            allele: IndelAllele::Insertion(bases),
            ot,
            ob,
            unknown_strand: 0,
            noisy: 0,
        }
    }

    fn deletion(seq: &str, ot: u32, ob: u32) -> IndelAlleleCounts {
        let bases: SmallVec<Base, 4> = seq.bytes().map(Base::from).collect();
        IndelAlleleCounts {
            allele: IndelAllele::Deletion(bases),
            ot,
            ob,
            unknown_strand: 0,
            noisy: 0,
        }
    }

    /// Locus whose *background* coverage is strand balanced, so that a single
    /// allele's null works out to 0.5 and the expected p-values are the plain
    /// `2 * 0.5^n`. Balancing the total instead would leave the null tilted
    /// against whichever strand the allele sits on, since the null is taken over
    /// the non-supporting fragments. Cases that need a skewed background
    /// override `ot_depth`/`ob_depth` directly.
    fn counts(ref_count: u32, alleles: Vec<IndelAlleleCounts>) -> IndelCounts {
        let alleles: SmallVec<IndelAlleleCounts, 2> = alleles.into_iter().collect();
        let alt_ot: u32 = alleles.iter().map(|a| a.ot).sum();
        let alt_ob: u32 = alleles.iter().map(|a| a.ob).sum();
        IndelCounts {
            alleles,
            ot_depth: ref_count.div_ceil(2) + alt_ot,
            ob_depth: ref_count / 2 + alt_ob,
            ref_count,
            ..Default::default()
        }
    }

    /// The default strand gate is both-strand concordance; the significance test is
    /// opt-in, so tests that exercise it have to switch it on.
    /// `significance_test_is_off_by_default` pins the default.
    fn with_strand_test() -> IndelParams {
        IndelParams { indel_strand_bias_alpha: 0.05, ..IndelParams::default() }
    }

    fn verdicts(calls: &[IndelCall]) -> Vec<Option<IndelVerdict>> {
        calls.iter().map(|c| c.hard_filter_verdict).collect()
    }

    #[test]
    fn heterozygous_indel_on_both_strands_passes() {
        // 5 alt / 10 depth on both strands → binomial het → PASS.
        let c = counts(5, vec![insertion("A", 3, 2)]);
        let calls = call_indels(&c, &IndelParams::default());
        assert_eq!(verdicts(&calls), vec![Some(IndelVerdict::Pass)]);
        assert!(calls[0].genotype.is_heterozygous());
    }

    #[test]
    fn single_stranded_indel_fails_concordance() {
        // A single-stranded allele fails concordance by default; no significance
        // test needed.
        let c = counts(20, vec![insertion("A", 8, 0)]);
        let calls = call_indels(&c, &IndelParams::default());
        assert_eq!(verdicts(&calls), vec![Some(IndelVerdict::FailStrand)]);
    }

    /// Known cost of concordance at low AO: a 2/0 split is a coin flip under a
    /// balanced null, yet it fails here for being single-stranded. The significance
    /// test deliberately kept it; concordance does not. See the module docs.
    #[test]
    fn two_supporting_fragments_on_one_strand_fail_concordance() {
        let c = counts(10, vec![insertion("A", 2, 0)]);
        let calls = call_indels(&c, &IndelParams::default());
        assert_eq!(verdicts(&calls), vec![Some(IndelVerdict::FailStrand)]);
    }

    /// Under concordance every one-sided split fails, independent of AO or alpha —
    /// unlike the significance test, which only rejected 6/0 and up. Pinned so the
    /// behavioural change from the significance gate stays explicit.
    #[test]
    fn concordance_rejects_every_one_sided_split() {
        let verdict_for = |ot: u32| -> Option<IndelVerdict> {
            let calls = call_indels(&counts(30, vec![insertion("A", ot, 0)]), &IndelParams::default());
            calls.first()?.hard_filter_verdict
        };
        for ot in 2..=8 {
            assert_eq!(verdict_for(ot), Some(IndelVerdict::FailStrand), "{ot}/0 is single-stranded");
        }
    }

    /// Concordance depends only on the allele's own OT/OB support, not the locus
    /// background — so the same one-sided allele fails whether surrounding coverage
    /// is balanced or skewed. (The significance test used the background as its null
    /// and could pass the skewed case; concordance cannot.)
    #[test]
    fn concordance_ignores_locus_background() {
        let balanced = counts(10, vec![insertion("A", 8, 0)]);
        assert_eq!(
            verdicts(&call_indels(&balanced, &IndelParams::default())),
            vec![Some(IndelVerdict::FailStrand)]
        );

        let mut skewed = balanced.clone();
        skewed.ot_depth = 17;
        skewed.ob_depth = 1;
        assert_eq!(
            verdicts(&call_indels(&skewed, &IndelParams::default())),
            vec![Some(IndelVerdict::FailStrand)]
        );
    }

    /// Only concordance runs by default; the significance test is opt-in. A
    /// two-sided but skewed 18/2 allele satisfies concordance and passes by default,
    /// but the opt-in significance test rejects it.
    #[test]
    fn significance_test_is_off_by_default() {
        assert_eq!(IndelParams::default().indel_strand_bias_alpha, 0.0);

        let c = counts(20, vec![insertion("A", 18, 2)]);
        assert_eq!(
            verdicts(&call_indels(&c, &IndelParams::default())),
            vec![Some(IndelVerdict::Pass)]
        );
        assert_eq!(
            verdicts(&call_indels(&c, &with_strand_test())),
            vec![Some(IndelVerdict::FailStrand)]
        );
    }

    #[test]
    fn zero_alpha_disables_only_the_significance_test() {
        // A two-sided skewed allele the significance test would reject passes at
        // alpha 0; concordance is satisfied because it is present on both strands.
        let c = counts(20, vec![insertion("A", 12, 2)]);
        let params = IndelParams { indel_strand_bias_alpha: 0.0, ..IndelParams::default() };
        assert_ne!(verdicts(&call_indels(&c, &params)), vec![Some(IndelVerdict::FailStrand)]);
    }

    #[test]
    fn low_vaf_indel_genotyped_hom_ref_fails() {
        // 2 alt in a depth of 32 (both strands) reads as sequencing noise → hom ref.
        let c = counts(30, vec![insertion("A", 1, 1)]);
        let calls = call_indels(&c, &IndelParams::default());
        assert_eq!(verdicts(&calls), vec![Some(IndelVerdict::FailHomRef)]);
    }

    #[test]
    fn below_min_ao_is_dropped_not_emitted() {
        // Sub-threshold alt observations are never candidates.
        let c = counts(20, vec![insertion("A", 1, 0)]);
        let calls = call_indels(&c, &IndelParams::default());
        assert!(calls.is_empty());
    }

    #[test]
    fn below_min_depth_locus_is_dropped() {
        let c = counts(0, vec![insertion("A", 1, 0)]);
        let calls = call_indels(&c, &IndelParams::default());
        assert!(calls.is_empty());
    }

    #[test]
    fn noise_reduces_genotyping_depth() {
        // 1 ref + 2 alt = depth 3, but 2 fragments are noise (soft-clip/repeat)
        // → clean depth 1 < min_indel_depth, locus dropped.
        let mut c = counts(1, vec![insertion("A", 1, 1)]);
        c.noisy_ref_count = 1;
        c.alleles[0].noisy = 1;
        let calls = call_indels(&c, &IndelParams::default());
        assert!(calls.is_empty());
    }

    /// The showstopper this replaced. Noise is a property of the *read*, so it
    /// lands on supporting and non-supporting fragments alike; subtracting only
    /// the non-supporting ones raises VAF by roughly the noise rate. The binomial
    /// flips from hom-ref to het at VAF ≈ 0.218, so a locus anywhere below that
    /// can be walked across the boundary by the haircut alone.
    ///
    /// 6 alt in 40 fragments is VAF 0.15, hom-ref. Half of every fragment class is
    /// noise. Excluding both sides leaves 3-in-20 — still 0.15, still hom-ref.
    #[test]
    fn noise_exclusion_does_not_inflate_vaf() {
        let mut c = counts(34, vec![insertion("A", 3, 3)]);
        c.noisy_ref_count = 17;
        c.alleles[0].noisy = 3;

        let calls = call_indels(&c, &IndelParams::default());
        assert_eq!(verdicts(&calls), vec![Some(IndelVerdict::FailHomRef)]);
        assert_eq!(calls[0].depth, 20);
        assert_eq!(calls[0].alt_count, 3);

        // What the reference-only haircut produced from the same locus: the 17
        // noisy reference fragments leave, all 6 supporting ones stay, and 6-in-23
        // (VAF 0.26) is over the boundary.
        let one_sided = counts(17, vec![insertion("A", 3, 3)]);
        let calls = call_indels(&one_sided, &IndelParams::default());
        assert_eq!(calls[0].depth, 23);
        assert_eq!(verdicts(&calls), vec![Some(IndelVerdict::Pass)]);
        assert!(calls[0].genotype.is_heterozygous());
    }

    /// The same locus under all three noise settings, so the axes they differ on
    /// are visible side by side: 6 supporting fragments of which 5 are noisy, in a
    /// locus of 40 fragments half of which are noisy.
    ///
    /// `symmetric` charges the noise twice — 1 clean observation is below the
    /// min-AO gate, so there is no candidate at all. `ratio-only` charges it once:
    /// the allele is a candidate on its 6 raw observations but genotyped on 1-in-20.
    /// `depth-only` charges it to the denominator alone, which turns VAF 0.15 into
    /// 6-in-23 = 0.26 and flips the call to het.
    #[test]
    fn the_three_noise_settings_differ_on_gate_and_ratio() {
        let locus = || {
            let mut c = counts(34, vec![insertion("A", 3, 3)]);
            c.noisy_ref_count = 17;
            c.alleles[0].noisy = 5;
            c
        };
        let with = |mode| IndelParams { indel_noise_exclusion: mode, ..IndelParams::default() };

        let symmetric = call_indels(&locus(), &with(IndelNoiseExclusion::Symmetric));
        assert!(symmetric.is_empty(), "1 clean observation is below the min-AO gate of 2");

        let ratio_only = call_indels(&locus(), &with(IndelNoiseExclusion::RatioOnly));
        assert_eq!(verdicts(&ratio_only), vec![Some(IndelVerdict::FailHomRef)]);
        assert_eq!((ratio_only[0].alt_count, ratio_only[0].depth), (1, 18));

        let depth_only = call_indels(&locus(), &with(IndelNoiseExclusion::DepthOnly));
        assert_eq!(verdicts(&depth_only), vec![Some(IndelVerdict::Pass)]);
        assert_eq!((depth_only[0].alt_count, depth_only[0].depth), (6, 23));
        assert!(depth_only[0].genotype.is_heterozygous());
    }

    /// Noisy fragments stay in the significance test's numerator even though they
    /// are excluded from genotyping. A two-sided but skewed allele with noisy
    /// support is still rejected by the opt-in significance test.
    #[test]
    fn significance_test_still_sees_noisy_support() {
        let mut c = counts(20, vec![insertion("A", 18, 2)]);
        c.alleles[0].noisy = 4;
        assert_eq!(
            verdicts(&call_indels(&c, &with_strand_test())),
            vec![Some(IndelVerdict::FailStrand)]
        );
    }

    /// Known cost of concordance: a homozygous indel at a locus with no OB coverage
    /// at all cannot satisfy the both-strand requirement, so it fails — even though
    /// it is a genuine call the significance test (whose null is the locus' own
    /// one-strand mix) would pass. Matters at capture edges and low coverage.
    #[test]
    fn homozygous_indel_at_a_one_strand_locus_fails_concordance() {
        let mut c = counts(0, vec![insertion("A", 9, 0)]);
        c.ot_depth = 9;
        c.ob_depth = 0;

        let calls = call_indels(&c, &IndelParams::default());
        assert_eq!(verdicts(&calls), vec![Some(IndelVerdict::FailStrand)]);
    }

    /// The counterpart: once there *are* non-supporting fragments on the other
    /// strand, the same one-sided support is judged against them and rejected.
    #[test]
    fn one_sided_support_against_other_strand_coverage_is_still_biased() {
        let c = counts(12, vec![insertion("A", 9, 0)]);
        assert_eq!(
            verdicts(&call_indels(&c, &with_strand_test())),
            vec![Some(IndelVerdict::FailStrand)]
        );
    }

    #[test]
    fn heterozygous_deletion_on_both_strands_passes() {
        // Deletions go through the same chain as insertions (all other cases use insertions).
        let c = counts(5, vec![deletion("AC", 3, 2)]);
        let calls = call_indels(&c, &IndelParams::default());
        assert_eq!(verdicts(&calls), vec![Some(IndelVerdict::Pass)]);
        assert!(calls[0].genotype.is_heterozygous());
        assert!(calls[0].allele.is_deletion());
    }

    #[test]
    fn min_ao_gate_is_inclusive_at_threshold() {
        // The gate is `alt_count < min_ao`, so exactly min_indel_ao (2) is kept and 1 is dropped.
        let at_threshold =
            call_indels(&counts(10, vec![insertion("A", 1, 1)]), &IndelParams::default());
        assert_eq!(at_threshold.len(), 1);
        let below = call_indels(&counts(10, vec![insertion("A", 1, 0)]), &IndelParams::default());
        assert!(below.is_empty());
    }

    /// Strand bias is judged on OT/OB, so it is unaffected by which mate of a
    /// fragment survived deduplication (both mates share an OT/OB assignment but
    /// have opposite reverse flags).
    #[test]
    fn strand_bias_is_on_ot_ob_not_read_direction() {
        let single = counts(20, vec![insertion("A", 7, 0)]);
        assert_eq!(
            verdicts(&call_indels(&single, &with_strand_test())),
            vec![Some(IndelVerdict::FailStrand)],
            "7/0 is p = 2*0.5^7 = 0.016"
        );

        let split = counts(20, vec![insertion("A", 6, 1)]);
        assert_eq!(
            verdicts(&call_indels(&split, &with_strand_test())),
            vec![Some(IndelVerdict::Pass)],
            "the same 7 fragments, one moved to the other strand, is p = 0.125"
        );
    }

    /// Under concordance, unknown-orientation fragments cannot satisfy the
    /// both-strand requirement — they are neither OT nor OB — so an allele that
    /// leans on them fails. (The significance test read them as no-evidence and
    /// could pass such an allele; another cost of the presence rule.)
    #[test]
    fn unknown_strand_cannot_satisfy_concordance() {
        let mut allele = insertion("A", 3, 0);
        allele.unknown_strand = 2;
        assert_eq!(allele.total(), 5, "unknown-strand fragments still count as support");
        let calls = call_indels(&counts(15, vec![allele]), &IndelParams::default());
        assert_eq!(verdicts(&calls), vec![Some(IndelVerdict::FailStrand)]);

        let mut opaque = insertion("A", 0, 0);
        opaque.unknown_strand = 9;
        let calls = call_indels(&counts(15, vec![opaque]), &IndelParams::default());
        assert_eq!(verdicts(&calls), vec![Some(IndelVerdict::FailStrand)]);
    }

    #[test]
    fn noise_counts_cannot_drive_depth_below_zero() {
        // Both noise counts are accumulated over the same fragments the depth comes
        // from, so each can equal but never exceed its side; saturating to zero must
        // still drop the locus rather than wrap.
        let mut c = counts(8, vec![insertion("A", 1, 1)]);
        c.noisy_ref_count = 100;
        c.alleles[0].noisy = 100;
        assert!(call_indels(&c, &IndelParams::default()).is_empty());
    }

    #[test]
    fn multiple_alleles_receive_independent_verdicts() {
        // Same locus: a strand-balanced het (Pass) and a strand-skewed one (FailStrand).
        let c = counts(6, vec![insertion("A", 3, 3), insertion("T", 8, 0)]);
        let calls = call_indels(&c, &with_strand_test());
        assert_eq!(calls.len(), 2);
        for call in &calls {
            match call.allele.bases().first().copied() {
                Some(Base::A) => assert_eq!(call.hard_filter_verdict, Some(IndelVerdict::Pass)),
                Some(Base::T) => {
                    assert_eq!(call.hard_filter_verdict, Some(IndelVerdict::FailStrand))
                }
                other => panic!("unexpected allele lead base {other:?}"),
            }
        }
    }
}
