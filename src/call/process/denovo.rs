#[cfg(doc)]
use crate::call::pair_hmm::Verdict;
use crate::call::pair_hmm::{Allele, Event, EventIdx, HapIdx};
use crate::{
    call::pileup::indels::IndelAllele,
    metrics::{Alt, PileupMetrics},
    vcf::CpgSide,
};
use enumset::EnumSet;
use seqair_types::{Base, Pos0};
use std::sync::Arc;
use tracing::instrument;

/// The CpGs called variants make with their neighbours: the called indels'
/// once the column loop has called them, so that the forests see them, and the
/// called alts' too once the alt calls are in.
///
/// A called `C` alt before a reference `G`, or a called `G` alt after a reference `C`, makes that
/// reference base the partner of a de-novo CpG ([`PileupMetrics::denovo_partner`]). So does a
/// called insertion ending in `C` before a reference `G`, or one starting with `G` after its
/// reference `C` anchor: the base beside it on the called chromosome decides. A called `C` alt
/// before a called `G` alt forms a CpG of two variants (`TA>CG`), which the reference context
/// cannot see: both alts form it ([`PileupMetrics::cpg_side_of`]), so that genotyping and
/// methylation weigh their reads as at any de-novo CpG. A CpG that is not the reference's needs its
/// `C` and `G` on one chromosome: neither two alts the pair-HMM puts apart, nor an alt and a
/// reference partner whose own called alt sits on every chromosome carrying the alt
/// ([`Verdict::has_haplotype`]).
///
/// Recomputed from the calls, so a test's re-run of the calls can clear it.
// shortcut: a deletion joining a `C` and a `G` across its gap makes no CpG here, as the base
// after the gap is no neighbour of its anchor's column; that needs a pass over the segment.
#[instrument(level = "debug", skip_all)]
pub fn mark_called_cpgs(
    before: Option<&PileupMetrics>,
    current: &mut PileupMetrics,
    after: Option<&PileupMetrics>,
) {
    let called_sides = |column: Option<&PileupMetrics>| -> EnumSet<CpgSide> {
        let alts = column.into_iter().flat_map(|c| c.alts.iter());
        alts.filter(|a| a.called()).filter_map(|a| CpgSide::of_base(a.base)).collect()
    };
    let inserted = |column: Option<&PileupMetrics>, end: fn(&[Base]) -> Option<&Base>| {
        let called = column.into_iter().flat_map(PileupMetrics::called_indels);
        let ends = called.filter_map(|allele| match allele {
            IndelAllele::Insertion(bases) => end(bases).copied(),
            IndelAllele::Deletion(_) => None,
        });
        ends.filter_map(CpgSide::of_base).collect::<EnumSet<CpgSide>>()
    };
    // The C and the G of a CpG that is not the reference's are on one chromosome or make none: a
    // called alt's partner, and a reference partner whose own called alt the other's chromosome
    // carries.
    let split_before = before.is_some_and(|c| !on_one_chromosome(c, current));
    let split_after = after.is_some_and(|g| !on_one_chromosome(current, g));
    let mut before_alts = called_sides(before);
    let mut after_alts = called_sides(after);
    current.lost_to_called_alts = EnumSet::empty();
    if split_before {
        before_alts.remove(CpgSide::C);
        if current.reference_base != Base::G {
            current.lost_to_called_alts.insert(CpgSide::G);
        }
    }
    if split_after {
        after_alts.remove(CpgSide::G);
        if current.reference_base != Base::C {
            current.lost_to_called_alts.insert(CpgSide::C);
        }
    }
    let before_sides = before_alts | inserted(before, <[Base]>::last);
    let after_sides = after_alts | inserted(Some(current), <[Base]>::first);
    current.beside_called_alts = CpgSide::beside(before_sides, after_sides);
}

/// Whether one chromosome reads `C` at `c` and `G` at the next column `g`, as the pair-HMM scored
/// them: by the event of an alt, by none of the column's events for the reference base, the two
/// columns' events present where their alts are called and the cluster's others where the
/// pair-HMM calls them. One does when they are not in one cluster, or a column has neither base.
// shortcut: as beside an indel that takes a cluster's slots; score their phase there too if such
// pairs turn out to matter.
fn on_one_chromosome(c: &PileupMetrics, g: &PileupMetrics) -> bool {
    let (Some(scores), Some(other)) = (c.snv_scores.as_ref(), g.snv_scores.as_ref()) else {
        return true;
    };
    let allele = |column: &PileupMetrics, base: Base| {
        if column.reference_base == base {
            Some(None)
        } else {
            scores.base_event_at(column.pos, base).map(Some)
        }
    };
    let (Some(at_c), Some(at_g)) = (allele(c, Base::C), allele(g, Base::G)) else {
        return true;
    };
    if !Arc::ptr_eq(scores, other) {
        return true;
    }
    let verdict = scores.verdict();
    let called =
        |column: &PileupMetrics, base: Base| column.find_alt(base).is_some_and(Alt::called);
    let held = |event: EventIdx| match scores.event(event) {
        Some(Event { anchor, allele: Allele::Base(base) }) if *anchor == c.pos => called(c, *base),
        Some(Event { anchor, allele: Allele::Base(base) }) if *anchor == g.pos => called(g, *base),
        _ => verdict.is_called(event),
    };
    // One substitution per anchor on a haplotype, so this is the base it reads.
    let reads = |hap: HapIdx, pos: Pos0, allele: Option<EventIdx>| {
        scores.base_events_at(pos).find(|&e| hap.carries(e)) == allele
    };
    let cpg = |hap: HapIdx| reads(hap, c.pos, at_c) && reads(hap, g.pos, at_g);
    verdict.has_haplotype(held, cpg).unwrap_or(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        call::{
            phasing::fixtures::{confirming, insertion},
            tests::fixtures::{column, with_indels},
        },
        metrics::ml::{features::features, types::MlModel},
    };
    use seqair_types::Base::*;

    /// A column over `reference` with `inserted` scored after it, six reads on
    /// each of the haplotypes `calls`: `[0, 0]` leaves it uncalled.
    fn anchor(pos: u32, reference: Base, inserted: &[Base], calls: [usize; 2]) -> PileupMetrics {
        let anchor = column(pos, reference, &[], None, &[]);
        with_indels(anchor, &[insertion(inserted)], &confirming(calls))
    }

    /// A reference `G` with an `A` alt, after `before`.
    fn g_after(before: &PileupMetrics) -> PileupMetrics {
        let pos = u32::try_from(u64::from(before.pos) + 1).expect("a test position");
        let mut g = column(pos, G, &[(A, false)], None, &[]);
        if let Some(alt) = g.find_alt_mut(A) {
            alt.metrics.base = A;
        }
        mark_called_cpgs(Some(before), &mut g, None);
        g
    }

    /// The `C` the insertion puts before the `G` makes it a CpG's `G`, so an
    /// `A` there is what its methylation reads as, for the CpG forest.
    #[test]
    fn a_called_insertion_ending_in_c_makes_the_next_g_a_cpg() {
        let g = g_after(&anchor(1040, A, &[T, C], [0, 1]));
        assert_eq!(g.denovo_partner(), Some(CpgSide::G));
        let alt = g.alt_metrics(A).expect("the A alt");
        assert!(alt.is_conversion());
        let (model, row) = features(&alt, None, None);
        assert_eq!(model, MlModel::Cpg);
        assert!(row.is_ok(), "{row:?}");
    }

    #[test]
    fn a_called_insertion_starting_with_g_makes_its_c_anchor_a_cpg() {
        let mut c = anchor(1040, C, &[G, T], [1, 1]);
        mark_called_cpgs(None, &mut c, None);
        assert_eq!(c.denovo_partner(), Some(CpgSide::C));
    }

    /// A `G` inserted between a reference CpG's `C` and `G` leaves the `C` a
    /// CpG's on both chromosomes, no de-novo one.
    #[test]
    fn an_insertion_inside_a_reference_cpg_makes_no_denovo_one() {
        let mut c = anchor(1040, C, &[G], [0, 1]);
        c.context.after_1 = Some(G);
        mark_called_cpgs(None, &mut c, None);
        assert_eq!(c.cpg_side_of(C), Some(CpgSide::C));
        assert_eq!(c.denovo_partner(), None);
        assert!(!c.in_denovo_cpg());
    }

    #[test]
    fn an_uncalled_insertion_or_one_ending_in_another_base_makes_no_cpg() {
        for before in [anchor(1040, A, &[C], [0, 0]), anchor(1040, A, &[C, A], [0, 1])] {
            let g = g_after(&before);
            assert_eq!(g.denovo_partner(), None);
            assert!(!g.alt_metrics(A).expect("the A alt").is_conversion());
        }
    }
}
