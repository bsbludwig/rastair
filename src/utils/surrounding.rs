use crate::metrics::PileupMetrics;
use seqair_types::Offset;
use std::sync::Arc;

/// Are these two pileups genuine neighbours in the genome?
///
/// By the time this runs the sequence has usually been filtered, so two entries
/// sitting next to each other in memory are often nowhere near each other on
/// the chromosome.
///
/// The contig check is `Arc::ptr_eq`, not a name comparison. Every pileup of a
/// region shares one region, so pointer identity settles it — and comparing
/// the `SmolStr` instead was 2.2 % of worker CPU in a `--gpu` chr12 profile,
/// spent establishing at every position that chr12 is chr12. Two pileups from
/// *different* `Arc`s are correctly not adjacent even if the contigs match by
/// name, which is the conservative direction: they cannot be neighbours in one
/// region's output anyway.
pub fn adjacent(left: &PileupMetrics, right: &PileupMetrics) -> bool {
    Arc::ptr_eq(&left.region, &right.region) && right.pos - left.pos == Offset::new(1)
}

/// Apply `f` to every pileup together with its immediate neighbours, in place.
///
/// `before` is the previous pileup *after* `f` has already run on it, `after`
/// the next one *before* it has, so a mutation propagates forwards — which is
/// what the de-novo CpG passes rely on. Either is `None` unless it is genuinely
/// adjacent.
///
/// Deliberately not an iterator adapter: yielding each pileup while also
/// keeping it as the next one's `before` would clone it once per pass.
pub fn map_surrounding<F>(pileups: &mut [PileupMetrics], mut f: F)
where
    F: FnMut(Option<&PileupMetrics>, &mut PileupMetrics, Option<&PileupMetrics>),
{
    for i in 0..pileups.len() {
        let (left, rest) = pileups.split_at_mut(i);
        // `i < len`, so `rest` is never empty; `break` keeps the loop total
        // without an unreachable branch.
        let Some((current, right)) = rest.split_first_mut() else { break };

        let before = left.last().filter(|previous| adjacent(previous, current));
        let after = right.first().filter(|next| adjacent(current, next));

        f(before, current, after);
    }
}

#[cfg(test)]
#[allow(clippy::cast_possible_truncation, reason = "test code")]
mod tests {
    use super::*;
    use crate::{
        call::tests::pileup::{Pileup, SimpleReads},
        sequence::{ChunkRegion, Region, Segment},
        vcf::SequenceContext,
    };
    use seqair_types::Base;
    use seqair_types::SmolStr;
    use std::cell::RefCell;
    use std::collections::HashMap;
    use std::rc::Rc;

    /// One shared region per contig, as production has it: every pileup of a
    /// region is built from the same `Arc`, which is what `adjacent` asks
    /// about. Building a fresh one per pileup would make every pair
    /// non-adjacent and the tests below vacuous.
    fn region_for(contig: &str, start: u64, end: u64) -> Arc<ChunkRegion> {
        thread_local! {
            static REGIONS: RefCell<HashMap<SmolStr, Arc<ChunkRegion>>> =
                RefCell::new(HashMap::new());
        }
        REGIONS.with(|regions| {
            Arc::clone(regions.borrow_mut().entry(contig.into()).or_insert_with(|| {
                Arc::new(ChunkRegion {
                    region: Region { contig: contig.into(), start, end },
                    last_position: end,
                    overlap_start: 0,
                    overlap_end: 0,
                })
            }))
        })
    }

    /// An uncovered `PileupMetrics` for testing: the C of a CpG, since only a
    /// CpG has metrics without reads.
    fn make_pileup(contig: &str, pos: u64) -> PileupMetrics {
        // The region is shared per contig, so the sequence is laid out from
        // its start, not from this position.
        let range = region_for(contig, pos.saturating_sub(10), pos + 20);
        let pos_idx = (pos - range.region.start) as usize;
        let mut sequence = vec![b'C'; pos_idx + 20];
        sequence[pos_idx + 1] = b'G';
        let segment = Rc::new(Segment { range, sequence });

        let context = SequenceContext::new(pos_idx, &segment).expect("valid context");

        let pileup = Pileup {
            region: segment.range.clone(),
            context,
            pos: pos as u32,
            reads: SimpleReads(vec![].into()),
            reference_base: Base::C,
        };

        pileup.metrics(&segment).unwrap()
    }

    #[test]
    fn an_empty_sequence_never_calls_the_mapper() {
        let mut items: Vec<PileupMetrics> = vec![];
        let mut calls = 0;
        map_surrounding(&mut items, |_, _, _| {
            calls += 1;
        });
        assert!(items.is_empty());
        assert_eq!(calls, 0);
    }

    #[test]
    fn a_lone_pileup_has_no_neighbours() {
        let mut items = vec![make_pileup("chr1", 100)];
        let mut calls = 0;
        map_surrounding(&mut items, |before, current, after| {
            calls += 1;
            assert!(before.is_none());
            assert!(after.is_none());
            assert_eq!(current.pos.as_u32(), 100);
        });
        assert_eq!(calls, 1);
        assert_eq!(items.len(), 1);
    }

    #[test]
    fn consecutive_pileups_see_each_other() {
        let mut items: Vec<_> = (0..10).map(|i| make_pileup("chr1", 100 + i)).collect();
        map_surrounding(&mut items, |before, current, after| {
            let pos = current.pos.as_u32();
            assert_eq!(before.map(|p| p.pos.as_u32()), (pos > 100).then(|| pos - 1));
            assert_eq!(after.map(|p| p.pos.as_u32()), (pos < 109).then(|| pos + 1));
        });
        assert_eq!(items.len(), 10);
    }

    /// Neighbouring in the vector is not neighbouring in the genome: the
    /// sequence has usually been filtered before this runs.
    #[test]
    fn gaps_and_contig_changes_break_adjacency() {
        let mut items = vec![
            make_pileup("chr1", 100),
            make_pileup("chr1", 101),
            make_pileup("chr1", 105), // gap
            make_pileup("chr2", 106), // consecutive position, different contig
            make_pileup("chr2", 107),
        ];

        let mut seen = Vec::new();
        map_surrounding(&mut items, |before, current, after| {
            seen.push((
                current.contig_name().to_owned(),
                current.pos.as_u32(),
                before.map(|p| p.pos.as_u32()),
                after.map(|p| p.pos.as_u32()),
            ));
        });

        assert_eq!(
            seen,
            vec![
                ("chr1".to_owned(), 100, None, Some(101)),
                ("chr1".to_owned(), 101, Some(100), None),
                ("chr1".to_owned(), 105, None, None),
                ("chr2".to_owned(), 106, None, Some(107)),
                ("chr2".to_owned(), 107, Some(106), None),
            ]
        );
    }

    #[test]
    fn mutations_are_kept() {
        let mut items =
            vec![make_pileup("chr1", 100), make_pileup("chr1", 101), make_pileup("chr1", 102)];
        map_surrounding(&mut items, |_, current, _| {
            current.pos = seqair_types::Pos0::new(current.pos.as_u32() + 1000).expect("a position");
        });
        assert_eq!(
            items.iter().map(|p| p.pos.as_u32()).collect::<Vec<_>>(),
            vec![1100, 1101, 1102]
        );
    }

    /// The mapper sees the *mutated* previous pileup and the *untouched* next
    /// one. The de-novo CpG passes depend on that direction.
    #[test]
    fn mutation_propagates_forwards_only() {
        let mut items =
            vec![make_pileup("chr1", 100), make_pileup("chr1", 101), make_pileup("chr1", 102)];
        let mut seen = Vec::new();
        map_surrounding(&mut items, |before, current, after| {
            seen.push((before.map(|p| p.pos_metrics.mapq0), after.map(|p| p.pos_metrics.mapq0)));
            current.pos_metrics.mapq0 = 7;
        });
        assert_eq!(seen, vec![(None, Some(0)), (Some(7), Some(0)), (Some(7), None)]);
    }
}
