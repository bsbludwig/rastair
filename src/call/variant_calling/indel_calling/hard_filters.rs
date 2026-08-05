//! Non-ML hard-filter indel calling.
//!
//! The pathway used by `rastair call --experimental-indels --no-ml`: candidate
//! indels are accepted by a fixed chain of hard filters — minimum filtered depth,
//! minimum alternate observations, forward/reverse strand concordance, and a
//! binomial genotype gate (homozygous-reference alleles rejected) — with no ML
//! scoring. Deliberately self-contained.
//!
//! To remove the whole hard-filter pathway, delete this file and the hook points
//! marked with a `hard-filter` comment:
//!   * `pub mod hard_filters;`, the `hard_filter_verdict` field on [`IndelCall`],
//!     and the `!ml_enabled` delegation — all in the parent `indel_calling.rs`
//!   * the verdict-rendering arm in `build_indel_records` (`src/vcf/metrics_to_vcf.rs`)
//!   * the `indel_strand` / `indel_hom_ref` FILTER ids (`src/vcf.rs`)
//!
//! Notes:
//!   * Genotyping depth subtracts a reference-read noise offset
//!     (`IndelCounts::ref_noise_offset`: reference reads with a terminal
//!     homopolymer/dinucleotide repeat or a soft-clip), kept separate from the
//!     ML-facing `depth_offset`.
//!   * The read-level mismatch filter is TAPS-aware (expected C→T on OT / G→A on
//!     OB conversions are not counted as mismatches).
//!   * Failing alleles are emitted with a FILTER tag rather than dropped, matching
//!     rastair's "always emit indels" convention.
//!   * Tie-break: `super::binomial_genotype` resolves exact-probability ties toward
//!     hom-ref (`>=`); such exact ties do not occur with continuous binomial masses.

use super::{IndelCall, IndelParams};
use crate::call::pileup::indels::IndelCounts;
use crate::call::variant_calling::GenotypeTag;

/// Verdict of the non-ML hard-filter indel chain.
///
/// Every allele that reaches genotyping is emitted; failures carry the reason so
/// `build_indel_records` can render the matching VCF FILTER tag.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum IndelVerdict {
    /// Passed every hard filter.
    Pass,
    /// Not supported on both the forward and reverse strands.
    FailStrand,
    /// Genotyped as homozygous reference by the binomial model.
    FailHomRef,
}

impl IndelVerdict {
    pub fn passes(self) -> bool {
        matches!(self, Self::Pass)
    }
}

/// Call indels with the non-ML hard-filter chain.
///
/// Shares the locus depth gate, per-allele min-AO and binomial genotype of
/// [`super::call_indels`], and adds forward/reverse strand concordance. Unlike the
/// ML path, every genotyped allele is emitted with an [`IndelVerdict`] instead of
/// being dropped.
pub fn call_indels(indels: &IndelCounts, params: &IndelParams) -> Vec<IndelCall> {
    let mut calls = Vec::new();

    let total_reads = indels.ref_count + indels.total_indel_reads();
    // Genotyping depth excludes reference reads flagged as noise (terminal
    // repeat / soft-clip); kept separate from the ML-facing `depth_offset`.
    let filtered_depth = total_reads.saturating_sub(indels.ref_noise_offset);

    if filtered_depth < params.min_indel_depth {
        return calls;
    }

    for allele_counts in &indels.alleles {
        let alt_count = allele_counts.total();
        if alt_count < params.min_indel_ao {
            continue;
        }

        let Some(genotype) = super::binomial_genotype(
            alt_count as usize,
            filtered_depth as usize,
            params.indel_error_rate,
        ) else {
            continue;
        };

        let verdict = if !allele_counts.on_both_strands() {
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
    use super::*;
    use crate::call::pileup::indels::{IndelAllele, IndelAlleleCounts};
    use seqair_types::{Base, SmallVec};

    fn insertion(seq: &str, fwd: u32, rev: u32) -> IndelAlleleCounts {
        let bases: SmallVec<Base, 4> = seq.bytes().map(Base::from).collect();
        IndelAlleleCounts { allele: IndelAllele::Insertion(bases), fwd, rev }
    }

    fn deletion(seq: &str, fwd: u32, rev: u32) -> IndelAlleleCounts {
        let bases: SmallVec<Base, 4> = seq.bytes().map(Base::from).collect();
        IndelAlleleCounts { allele: IndelAllele::Deletion(bases), fwd, rev }
    }

    fn counts(ref_count: u32, alleles: Vec<IndelAlleleCounts>) -> IndelCounts {
        IndelCounts { alleles: alleles.into_iter().collect(), ref_count, ..Default::default() }
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
    fn single_stranded_indel_fails_strand_concordance() {
        // Ample support but only on the forward strand → rejected.
        let c = counts(5, vec![insertion("A", 5, 0)]);
        let calls = call_indels(&c, &IndelParams::default());
        assert_eq!(verdicts(&calls), vec![Some(IndelVerdict::FailStrand)]);
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
    fn ref_noise_offset_reduces_genotyping_depth() {
        // 1 ref + 2 alt = depth 3, but 2 reference reads are noise (soft-clip/repeat)
        // → filtered depth 1 < min_indel_depth, locus dropped.
        let mut c = counts(1, vec![insertion("A", 1, 1)]);
        c.ref_noise_offset = 2;
        let calls = call_indels(&c, &IndelParams::default());
        assert!(calls.is_empty());
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
        let at_threshold = call_indels(&counts(10, vec![insertion("A", 1, 1)]), &IndelParams::default());
        assert_eq!(at_threshold.len(), 1);
        let below = call_indels(&counts(10, vec![insertion("A", 1, 0)]), &IndelParams::default());
        assert!(below.is_empty());
    }

    #[test]
    fn multiple_alleles_receive_independent_verdicts() {
        // Same locus: a both-strand het (Pass) and a single-strand allele (FailStrand).
        let c = counts(6, vec![insertion("A", 3, 3), insertion("T", 4, 0)]);
        let calls = call_indels(&c, &IndelParams::default());
        assert_eq!(calls.len(), 2);
        for call in &calls {
            match call.allele.bases().first().copied() {
                Some(Base::A) => assert_eq!(call.hard_filter_verdict, Some(IndelVerdict::Pass)),
                Some(Base::T) => assert_eq!(call.hard_filter_verdict, Some(IndelVerdict::FailStrand)),
                other => panic!("unexpected allele lead base {other:?}"),
            }
        }
    }
}
