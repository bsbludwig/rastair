use crate::{
    call::pileup::indels::IndelAllele,
    metrics::{PairedCounts, PileupMetrics, ReadKey},
    vcf::CpgSide,
};
use seqair_types::{Base, Strand};
use tracing::instrument;

/// What called variants next door put beside a column on each of its two chromosomes
/// ([`PileupMetrics::called_before`], [`PileupMetrics::called_after`]): the called indels' once
/// the column loop has called them, so that the forests see the CpGs they make, and the called
/// alts' too once the alt calls are in.
///
/// A CpG is a `C` and a `G` next to each other on one chromosome, so every CpG question
/// ([`PileupMetrics::cpg_side_of`]) is asked of the two chromosomes' bases at a column and the
/// one beside it. The genotypes give each column's two bases; where both columns are
/// heterozygous the reads covering both say which base goes with which ([`crossed`]). So a
/// called `C` alt before a reference `G` makes a CpG, a `TA>CG` on one chromosome makes one and
/// leaves the other chromosome's reads out of it, a `C>T` and a `G>A` on different chromosomes
/// leave none of a reference CpG, and neither does a homozygous alt.
///
/// A called insertion ending in `C` before a `G`, or one starting with `G` after its `C` anchor,
/// puts that base beside the column on one chromosome.
///
/// Recomputed from the calls, so a test's re-run of the calls can clear it.
// shortcut: an insertion is put on one chromosome whatever its genotype, and the reference base
// stays beside the column on the other; a homozygous insertion inside a reference CpG still
// leaves it one. Read the insertion's genotype if such CpGs turn out to matter.
// shortcut: a deletion joining a `C` and a `G` across its gap makes no CpG here, as the base
// after the gap is no neighbour of its anchor's column; that needs a pass over the segment.
#[instrument(level = "debug", skip_all)]
pub fn mark_called_cpgs(
    before: Option<&PileupMetrics>,
    current: &mut PileupMetrics,
    after: Option<&PileupMetrics>,
) {
    let inserts = |anchor: Option<&PileupMetrics>, end: fn(&[Base]) -> Option<&Base>, base| {
        anchor.into_iter().flat_map(PileupMetrics::called_indels).any(|allele| match allele {
            IndelAllele::Insertion(bases) => end(bases) == Some(&base),
            IndelAllele::Deletion(_) => false,
        })
    };
    let in_order = |[first, second]: [Base; 2], crossed| {
        if crossed { [second, first] } else { [first, second] }
    };
    let own = current.chromosomes();
    current.called_before = before.and_then(|column| {
        let bases = column.chromosomes()?;
        Some(in_order(bases, own.is_some_and(|own| crossed(bases, own, &column.after_counts))))
    });
    current.called_after = after.and_then(|column| {
        let bases = column.chromosomes()?;
        Some(in_order(bases, own.is_some_and(|own| crossed(own, bases, &current.after_counts))))
    });
    let beside = |column: &PileupMetrics, side| column.beside(side).unwrap_or([Base::Unknown; 2]);
    if inserts(before, <[Base]>::last, Base::C) {
        current.called_before = Some(inserted(own, beside(current, CpgSide::G), CpgSide::G));
    }
    if inserts(Some(current), <[Base]>::first, Base::G) {
        current.called_after = Some(inserted(own, beside(current, CpgSide::C), CpgSide::C));
    }
}

/// `beside` with `side`'s partner inserted on one chromosome: the one whose base here takes the
/// side, where the partner is not beside it already.
fn inserted(own: Option<[Base; 2]>, [first, second]: [Base; 2], side: CpgSide) -> [Base; 2] {
    let partner = side.adjacent_base();
    match own {
        Some([base, _]) if base == side.unmod_base() && first != partner => [partner, second],
        _ => [first, partner],
    }
}

/// Whether the chromosome reading `first`'s first base reads `next`'s second: `first` and `next`
/// are two adjacent columns' [`PileupMetrics::chromosomes`], `reads` the first's
/// [`PileupMetrics::after_counts`].
///
/// Only two heterozygous columns leave the question open, and every read covering one covers
/// the other. Each answer says which two-base reads can exist; the one leaving fewer reads
/// unexplained holds, and a tie goes to the genotypes' own order, the two alts together, as
/// neighbouring variants mostly arise.
fn crossed(first: [Base; 2], next: [Base; 2], reads: &PairedCounts) -> bool {
    let ([a, b], [c, d]) = (first, next);
    if a == b || c == d {
        return false;
    }
    unexplained(reads, [[a, d], [b, c]]) < unexplained(reads, [[a, c], [b, d]])
}

/// The reads of two adjacent columns that neither of `chromosomes` (each the bases one reads at
/// the two columns) would show.
fn unexplained(reads: &PairedCounts, chromosomes: [[Base; 2]; 2]) -> u32 {
    let pairs = Base::KNOWN.into_iter().flat_map(|x| Base::KNOWN.map(|y| [x, y]));
    let keys = pairs.flat_map(|read| [Strand::OT, Strand::OB].map(|strand| (strand, read)));
    keys.filter(|&(strand, read)| !chromosomes.iter().any(|&bases| shows(strand, bases, read)))
        .map(|(strand, [current, adj])| reads.get(ReadKey { strand, current, adj }))
        .sum()
}

/// Whether a read on `strand` of a chromosome reading `bases` at two adjacent columns can show
/// `read` there: the bases as they are, or a CpG's converted on the strand that converts it
/// (`TG` on OT, `CA` on OB).
fn shows(strand: Strand, bases: [Base; 2], read: [Base; 2]) -> bool {
    bases == read
        || bases == [Base::C, Base::G]
            && match strand {
                Strand::OT => read == [Base::T, Base::G],
                Strand::OB => read == [Base::C, Base::A],
            }
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
    use hegel::{TestCase, generators as gs};
    use seqair_types::Base::*;

    /// Reads of both chromosomes on both strands, a CpG's converted or not,
    /// say which of two heterozygous columns' bases share a chromosome.
    #[hegel::test]
    fn the_reads_of_both_chromosomes_say_which_bases_share_one(tc: TestCase) {
        let any = || gs::sampled_from(Base::KNOWN.as_slice());
        let other = |base: Base| {
            gs::sampled_from(Base::KNOWN.into_iter().filter(|&b| b != base).collect::<Vec<_>>())
        };
        let (a, c) = (tc.draw(any()), tc.draw(any()));
        let (b, d) = (tc.draw(other(a)), tc.draw(other(c)));
        let is_crossed = tc.draw(gs::booleans());
        let chromosomes = if is_crossed { [[a, d], [b, c]] } else { [[a, c], [b, d]] };

        let mut reads = PairedCounts::default();
        for strand in [Strand::OT, Strand::OB] {
            for bases in chromosomes {
                for _ in 0..tc.draw(gs::integers::<u8>().min_value(1).max_value(5)) {
                    let [current, adj] = match (bases, strand, tc.draw(gs::booleans())) {
                        ([C, G], Strand::OT, true) => [T, G],
                        ([C, G], Strand::OB, true) => [C, A],
                        _ => bases,
                    };
                    reads.increment(ReadKey { strand, current, adj });
                }
            }
        }

        assert_eq!(crossed([a, b], [c, d], &reads), is_crossed);
    }

    /// One heterozygous column leaves nothing to decide, whatever the reads.
    #[test]
    fn a_homozygous_column_is_never_crossed() {
        let mut reads = PairedCounts::default();
        reads.increment(ReadKey { strand: Strand::OT, current: T, adj: A });
        assert!(!crossed([C, C], [G, A], &reads));
        assert!(!crossed([C, T], [G, G], &reads));
    }

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
