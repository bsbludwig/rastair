#[cfg(doc)]
use crate::call::pair_hmm::Verdict;
use crate::call::{
    pair_hmm::{Allele, Event, EventIdx, HapIdx, HaplotypeScores},
    phasing::PhaseSet,
};
use crate::{
    call::pileup::indels::IndelAllele,
    metrics::{Alt, PileupMetrics},
    vcf::CpgSide,
};
use enumset::EnumSet;
use seqair_types::{Base, Pos0};
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
/// methylation weigh their reads as at any de-novo CpG. Every CpG, the reference's too, needs its
/// `C` and `G` on one chromosome ([`Verdict::has_haplotype`]): not two alts the pair-HMM puts
/// apart, an alt and a reference partner whose own called alt sits on every chromosome carrying
/// the alt, or a reference CpG whose called alts leave no chromosome with both bases. A side whose methylated-read base no chromosome reads beside the
/// partner is [`PileupMetrics::unconfounded`].
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
    // A CpG's C and G are on one chromosome or make none: two called alts, an alt and a reference
    // partner, or the reference's own two bases where called alts replace them.
    let split_before = before.is_some_and(|c| !on_one_chromosome(c, current));
    let split_after = after.is_some_and(|g| !on_one_chromosome(current, g));
    let mut before_alts = called_sides(before);
    let mut after_alts = called_sides(after);
    current.lost_to_called_alts = EnumSet::empty();
    if split_before {
        before_alts.remove(CpgSide::C);
        current.lost_to_called_alts.insert(CpgSide::G);
    }
    if split_after {
        after_alts.remove(CpgSide::G);
        current.lost_to_called_alts.insert(CpgSide::C);
    }
    let before_sides = before_alts | inserted(before, <[Base]>::last);
    let after_sides = after_alts | inserted(Some(current), <[Base]>::first);
    current.beside_called_alts = CpgSide::beside(before_sides, after_sides);
    current.unconfounded = EnumSet::empty();
    if after.is_some_and(|g| !one_chromosome_reads(current, CpgSide::C.mod_base(), g, Base::G)) {
        current.unconfounded.insert(CpgSide::C);
    }
    if before.is_some_and(|c| !one_chromosome_reads(c, Base::C, current, CpgSide::G.mod_base())) {
        current.unconfounded.insert(CpgSide::G);
    }
}

/// Whether one chromosome reads `C` at `c` and `G` at the next column `g`.
fn on_one_chromosome(c: &PileupMetrics, g: &PileupMetrics) -> bool {
    one_chromosome_reads(c, Base::C, g, Base::G)
}

/// Whether one chromosome reads `first_base` at `first` and `next_base` at the next column `next`.
/// A column without a called alt reads its reference on every chromosome; otherwise the pair-HMM
/// decides, in a cluster of either column holding both columns' called alts: by the event of an
/// alt, by none of the column's events for the reference base, the two columns' events present
/// where their alts are called and the cluster's others where the pair-HMM calls them. Without
/// such a cluster the calls decide ([`from_calls`]), and one chromosome does when they cannot.
// shortcut: two hets in no common phase block are taken to be on one chromosome; score the pair
// in a cluster of its own if such pairs turn out to matter.
fn one_chromosome_reads(
    first: &PileupMetrics,
    first_base: Base,
    next: &PileupMetrics,
    next_base: Base,
) -> bool {
    let varies = |column: &PileupMetrics| column.alts.iter().any(Alt::called);
    let columns = [(first, first_base), (next, next_base)];
    if columns.iter().any(|&(column, base)| !varies(column) && column.reference_base != base) {
        return false;
    }
    let holds = |scores: &HaplotypeScores, column: &PileupMetrics| {
        let called = column.alts.iter().filter(|alt| alt.called());
        called.into_iter().all(|alt| scores.base_event_at(column.pos, alt.base).is_some())
    };
    let clusters = [first, next].into_iter().filter_map(|column| column.snv_scores.as_deref());
    let Some(scores) = clusters
        .into_iter()
        .find(|scores| columns.iter().all(|&(column, _)| holds(scores, column)))
    else {
        return from_calls(first, first_base, next, next_base).unwrap_or(true);
    };
    let allele = |column: &PileupMetrics, base: Base| {
        if column.reference_base == base {
            Some(None)
        } else {
            scores.base_event_at(column.pos, base).map(Some)
        }
    };
    let (Some(at_first), Some(at_next)) = (allele(first, first_base), allele(next, next_base))
    else {
        return false;
    };
    let verdict = scores.verdict();
    let called =
        |column: &PileupMetrics, base: Base| column.find_alt(base).is_some_and(Alt::called);
    let held = |event: EventIdx| match scores.event(event) {
        Some(Event { anchor, allele: Allele::Base(base) }) if *anchor == first.pos => {
            called(first, *base)
        }
        Some(Event { anchor, allele: Allele::Base(base) }) if *anchor == next.pos => {
            called(next, *base)
        }
        _ => verdict.is_called(event),
    };
    // One substitution per anchor on a haplotype, so this is the base it reads.
    let reads = |hap: HapIdx, pos: Pos0, allele: Option<EventIdx>| {
        scores.base_events_at(pos).find(|&e| hap.carries(e)) == allele
    };
    let both = |hap: HapIdx| reads(hap, first.pos, at_first) && reads(hap, next.pos, at_next);
    verdict
        .has_haplotype(held, both)
        .or_else(|| from_calls(first, first_base, next, next_base))
        .unwrap_or(true)
}

/// [`one_chromosome_reads`] from the two columns' genotypes, where one is homozygous, and from their
/// phase where both are hets in one phase block; `None` otherwise.
fn from_calls(
    first: &PileupMetrics,
    first_base: Base,
    next: &PileupMetrics,
    next_base: Base,
) -> Option<bool> {
    let (a, a_set) = chromosomes(first)?;
    let (b, b_set) = chromosomes(next)?;
    if a[0] == a[1] || b[0] == b[1] {
        return Some(a.contains(&first_base) && b.contains(&next_base));
    }
    if a_set.is_none() || a_set != b_set {
        return None;
    }
    Some(a.iter().zip(&b).any(|(&x, &y)| x == first_base && y == next_base))
}

/// The bases `column`'s two chromosomes carry by its calls, in haplotype order for a het in a
/// phase block (with the block), the reference on both without a called alt; `None` without a
/// genotype for its called alts.
fn chromosomes(column: &PileupMetrics) -> Option<([Base; 2], Option<PhaseSet>)> {
    if !column.alts.iter().any(Alt::called) {
        return Some(([column.reference_base; 2], None));
    }
    let genotype = column.pos_metrics.genotype.as_ref()?.genotype;
    let [a, b] = genotype.alleles();
    let phase = column.phase.as_ref().and_then(|phase| phase.snv);
    let (one, two) = match phase {
        Some(call) if genotype.is_heterozygous() => {
            (call.first, if call.first == a { b } else { a })
        }
        _ => (a, b),
    };
    let bases = [column.allele_base(one)?, column.allele_base(two)?];
    Some((bases, phase.map(|call| call.set)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        call::{
            phasing::{
                RecordPhases,
                fixtures::{call, confirming, insertion},
            },
            tests::fixtures::{column, with_indels},
            variant_calling::GenotypeTag,
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

    /// A called SNV column with no pair-HMM scores, so only its calls can say
    /// what its chromosomes read: `phase` is the alt's haplotype in block 1001.
    fn called(pos: u32, reference: Base, alt: Base, hom: bool, phase: Option<u8>) -> PileupMetrics {
        let one = std::num::NonZeroU8::MIN;
        let genotype = if hom { GenotypeTag::HomAlt(one) } else { GenotypeTag::RefHet(one) };
        let mut c = column(pos, reference, &[(alt, true)], Some(genotype), &[]);
        c.phase = phase.map(|alt_on| {
            let first = u8::from(alt_on == 1);
            Box::new(RecordPhases { snv: Some(call(1001, first, 30)), indel: None })
        });
        c
    }

    #[test]
    fn a_homozygous_neighbour_decides_without_a_cluster() {
        // `1/1 T` over a reference `C`: no chromosome reads the `C`.
        let lost = called(1000, C, T, true, None);
        assert!(!one_chromosome_reads(&lost, C, &called(1001, A, G, false, None), G));
        // `1/1 C` beside a het `G`: the `G`'s chromosome has the `C` too.
        let both = called(1000, T, C, true, None);
        assert!(one_chromosome_reads(&both, C, &called(1001, A, G, false, None), G));
    }

    #[test]
    fn two_hets_in_one_phase_block_decide_by_their_phase() {
        let c = called(1000, T, C, false, Some(1));
        assert!(one_chromosome_reads(&c, C, &called(1001, A, G, false, Some(1)), G));
        assert!(!one_chromosome_reads(&c, C, &called(1001, A, G, false, Some(2)), G));
    }
}
