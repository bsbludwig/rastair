use crate::metrics::PileupMetrics;
use color_eyre::Result;
use tracing::warn;

/// Are these two pileups genuine neighbours in the genome?
///
/// By the time this runs the sequence has usually been filtered, so two entries
/// sitting next to each other in memory are often nowhere near each other on
/// the chromosome.
fn adjacent(left: &PileupMetrics, right: &PileupMetrics) -> bool {
    left.contig_name() == right.contig_name() && left.pos.checked_add(1) == Some(right.pos)
}

/// Apply `f` to every pileup together with its immediate neighbours, in place.
///
/// `before` is the previous pileup *after* `f` has already run on it, `after`
/// the next one *before* it has, so a mutation propagates forwards — which is
/// what the de-novo CpG passes rely on. Either is `None` unless it is genuinely
/// adjacent.
///
/// A pileup `f` fails on is logged with `context` and dropped, but only once
/// the pass is over: until then it stays visible as a neighbour, which is what
/// the sliding window this replaced did. Errors never stop the pass.
///
/// Deliberately not an iterator adapter. Yielding each pileup while also
/// keeping it as the next one's `before` meant cloning all 944 bytes of it once
/// per pass, and the window itself moved each pileup three times more; that
/// `memcpy` was 8 % of worker CPU.
pub fn map_surrounding<F>(pileups: &mut Vec<PileupMetrics>, mut f: F, context: &str)
where
    F: FnMut(Option<&PileupMetrics>, &mut PileupMetrics, Option<&PileupMetrics>) -> Result<()>,
{
    let mut failed: Vec<usize> = Vec::new();

    for i in 0..pileups.len() {
        let (left, rest) = pileups.split_at_mut(i);
        // `i < len`, so `rest` is never empty; `break` keeps the loop total
        // without an unreachable branch.
        let Some((current, right)) = rest.split_first_mut() else { break };

        let before = left.last().filter(|previous| adjacent(previous, current));
        let after = right.first().filter(|next| adjacent(current, next));

        if let Err(error) = f(before, current, after) {
            warn!(error = format!("{error:#}"), "{context}");
            failed.push(i);
        }
    }

    if failed.is_empty() {
        return;
    }
    // `failed` is ascending, so one pass in step with `retain` is enough.
    let mut failed = failed.into_iter().peekable();
    let mut idx = 0usize;
    pileups.retain(|_| {
        let keep = failed.peek() != Some(&idx);
        if !keep {
            failed.next();
        }
        idx += 1;
        keep
    });
}

#[cfg(test)]
#[allow(clippy::cast_possible_truncation, reason = "test code")]
mod tests {
    use super::*;
    use crate::{
        call::pileup::{Pileup, SimpleReads},
        sequence::{ChunkRegion, Region, Segment},
        vcf::SequenceContext,
    };
    use seqair_types::Base;
    use std::rc::Rc;

    /// Helper to create a minimal `PileupMetrics` for testing
    fn make_pileup(contig: &str, pos: u64) -> PileupMetrics {
        // Create a minimal segment with enough context
        let start = pos.saturating_sub(10);
        let end = pos + 20;
        let segment = Rc::new(Segment {
            range: ChunkRegion {
                region: Region { contig: contig.into(), start, end },
                last_position: end,
                overlap_start: 0,
                overlap_end: 0,
            },
            sequence: vec![b'C'; (end - start) as usize],
            overlap_start: 0,
            overlap_end: 0,
        });

        let pos_idx = (pos - start) as usize;
        let context = SequenceContext::new(pos_idx, &segment).expect("valid context");

        let pileup = Pileup {
            region: segment.range.clone(),
            context,
            pos: pos as u32,
            reads: SimpleReads(vec![].into()),
            reference_base: Base::C,
            indel_observations: Default::default(),
            noisy_ref_count: 0,
            homopolymer_run: 0,
            dinucleotide_run: 0,
            soft_clip_count: 0,
            indel_ref_window: Default::default(),
            indel_ref_anchor: 0,
        };

        PileupMetrics::new(pileup).unwrap()
    }

    /// Run `f` over `items` and hand back what survived, so each test reads as
    /// "these went in, these came out".
    fn run<F>(items: &mut Vec<PileupMetrics>, f: F)
    where
        F: FnMut(Option<&PileupMetrics>, &mut PileupMetrics, Option<&PileupMetrics>) -> Result<()>,
    {
        map_surrounding(items, f, "test mapper failed");
    }

    #[test]
    fn an_empty_sequence_never_calls_the_mapper() {
        let mut items: Vec<PileupMetrics> = vec![];
        let mut calls = 0;
        run(&mut items, |_, _, _| {
            calls += 1;
            Ok(())
        });
        assert!(items.is_empty());
        assert_eq!(calls, 0);
    }

    #[test]
    fn a_lone_pileup_has_no_neighbours() {
        let mut items = vec![make_pileup("chr1", 100)];
        let mut calls = 0;
        run(&mut items, |before, current, after| {
            calls += 1;
            assert!(before.is_none());
            assert!(after.is_none());
            assert_eq!(current.pos, 100);
            Ok(())
        });
        assert_eq!(calls, 1);
        assert_eq!(items.len(), 1);
    }

    #[test]
    fn consecutive_pileups_see_each_other() {
        let mut items: Vec<_> = (0..10).map(|i| make_pileup("chr1", 100 + i)).collect();
        run(&mut items, |before, current, after| {
            let pos = current.pos;
            assert_eq!(before.map(|p| p.pos), (pos > 100).then(|| pos - 1));
            assert_eq!(after.map(|p| p.pos), (pos < 109).then(|| pos + 1));
            Ok(())
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
        run(&mut items, |before, current, after| {
            seen.push((
                current.contig_name().to_owned(),
                current.pos,
                before.map(|p| p.pos),
                after.map(|p| p.pos),
            ));
            Ok(())
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
        run(&mut items, |_, current, _| {
            current.pos += 1000;
            Ok(())
        });
        assert_eq!(items.iter().map(|p| p.pos).collect::<Vec<_>>(), vec![1100, 1101, 1102]);
    }

    /// The mapper sees the *mutated* previous pileup and the *untouched* next
    /// one. The de-novo CpG passes depend on that direction.
    #[test]
    fn mutation_propagates_forwards_only() {
        let mut items =
            vec![make_pileup("chr1", 100), make_pileup("chr1", 101), make_pileup("chr1", 102)];
        let mut seen = Vec::new();
        run(&mut items, |before, current, after| {
            seen.push((before.map(|p| p.pos_metrics.mapq0), after.map(|p| p.pos_metrics.mapq0)));
            current.pos_metrics.mapq0 = 7;
            Ok(())
        });
        assert_eq!(seen, vec![(None, Some(0)), (Some(7), Some(0)), (Some(7), None)]);
    }

    /// A failing pileup is dropped, but not before the one after it has had a
    /// chance to see it — and the pass runs to the end regardless.
    #[test]
    fn a_failing_pileup_is_dropped_after_serving_as_a_neighbour() {
        let mut items: Vec<_> = (0..5).map(|i| make_pileup("chr1", 100 + i)).collect();
        let mut saw_101_as_before = false;
        let mut calls = 0;

        run(&mut items, |before, current, _| {
            calls += 1;
            if before.is_some_and(|p| p.pos == 101) {
                saw_101_as_before = true;
            }
            if current.pos == 101 {
                color_eyre::eyre::bail!("simulated failure at 101");
            }
            Ok(())
        });

        assert_eq!(calls, 5, "a failure must not stop the pass");
        assert!(saw_101_as_before, "the failed pileup still served as a neighbour");
        assert_eq!(items.iter().map(|p| p.pos).collect::<Vec<_>>(), vec![100, 102, 103, 104]);
    }

    #[test]
    fn several_failures_are_all_dropped() {
        let mut items: Vec<_> = (0..6).map(|i| make_pileup("chr1", 100 + i)).collect();
        run(&mut items, |_, current, _| {
            if current.pos % 2 == 0 {
                color_eyre::eyre::bail!("simulated failure at {}", current.pos);
            }
            Ok(())
        });
        assert_eq!(items.iter().map(|p| p.pos).collect::<Vec<_>>(), vec![101, 103, 105]);
    }
}
