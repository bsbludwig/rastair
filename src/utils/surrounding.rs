use crate::metrics::PileupMetrics;
use std::fmt;

/// Struct holding references to current and its surrounding items
pub struct Surrounding<'a, T> {
    pub before: Option<&'a T>,
    pub current: &'a mut T,
    pub after: Option<&'a T>,
}

impl<'a, T: fmt::Debug> fmt::Debug for Surrounding<'a, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Surrounding")
            .field("before", &self.before)
            .field("current", &self.current)
            .field("after", &self.after)
            .finish()
    }
}

/// Get the pileup metrics for the record at `index` and its direct neighbors, if they exist
pub fn surrounding_pileups(
    records: &mut [PileupMetrics],
    index: usize,
) -> Surrounding<'_, PileupMetrics> {
    // To appease the borrow checker and get a mutable reference to the current record,
    // we split the records into three parts.
    let (left, right) = records.split_at_mut(index);
    let (current_slice, next_slice) = right.split_at_mut(1);
    let current = &mut current_slice[0];

    let before = left.last();
    let after = next_slice.first();
    // we might not have the direct neighbors
    let before = before.filter(|p| {
        p.pileup.contig() == current.pileup.contig()
            && Some(p.pileup.pos) == current.pileup.pos.checked_sub(1)
    });
    let after = after.filter(|p| {
        p.pileup.contig() == current.pileup.contig()
            && Some(p.pileup.pos) == current.pileup.pos.checked_add(1)
    });

    Surrounding { before, current, after }
}

/// Iterator adapter that processes elements with access to their validated neighbors.
///
/// This adapter maintains a sliding window of 3 elements and validates that neighbors
/// are actually adjacent (same contig, consecutive positions) before passing them to
/// the mapping function.
pub struct SurroundingMap<I, F>
where
    I: Iterator<Item = PileupMetrics>,
    F: FnMut(Option<&PileupMetrics>, &mut PileupMetrics, Option<&PileupMetrics>),
{
    iter: I,
    window: [Option<PileupMetrics>; 3],
    mapper: F,
    started: bool,
}

impl<I, F> SurroundingMap<I, F>
where
    I: Iterator<Item = PileupMetrics>,
    F: FnMut(Option<&PileupMetrics>, &mut PileupMetrics, Option<&PileupMetrics>),
{
    fn new(iter: I, mapper: F) -> Self {
        Self { iter, window: [None, None, None], mapper, started: false }
    }

    /// Validate if the element at window[idx] is a true neighbor of current
    fn is_valid_neighbor(&self, idx: usize, current: &PileupMetrics, is_before: bool) -> bool {
        let Some(neighbor) = self.window[idx].as_ref() else { return false };

        if neighbor.pileup.contig() != current.pileup.contig() {
            return false;
        }

        if is_before {
            // Check if neighbor.pos + 1 == current.pos
            neighbor.pileup.pos.checked_add(1) == Some(current.pileup.pos)
        } else {
            // Check if neighbor.pos - 1 == current.pos (i.e., current.pos + 1 == neighbor.pos)
            current.pileup.pos.checked_add(1) == Some(neighbor.pileup.pos)
        }
    }
}

impl<I, F> Iterator for SurroundingMap<I, F>
where
    I: Iterator<Item = PileupMetrics>,
    F: FnMut(Option<&PileupMetrics>, &mut PileupMetrics, Option<&PileupMetrics>),
{
    type Item = PileupMetrics;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.started {
            // Bootstrap: fill window as [None, first, second]
            // This allows the first element to have before=None
            self.window[1] = self.iter.next();
            self.window[2] = self.iter.next();
            self.started = true;
        } else {
            // Slide the window: pull a new element and shift
            let new_elem = self.iter.next();
            self.window[0] = self.window[1].take();
            self.window[1] = self.window[2].take();
            self.window[2] = new_elem;
        }

        // Load-bearing `?`: If window[1] is None, we're done
        let current = self.window[1].as_ref()?;

        // Phase 1: Validate neighbors while everything is still borrowed
        let before_valid = { self.is_valid_neighbor(0, current, true) };
        let after_valid = { self.is_valid_neighbor(2, current, false) };

        // Phase 2: Take current element to get mutable access
        let Some(mut current) = self.window[1].take() else {
            unreachable!("could not take current element even though we just had it");
        };

        // Phase 3: Call mapper with validated neighbor references
        let before_ref = if before_valid { self.window[0].as_ref() } else { None };
        let after_ref = if after_valid { self.window[2].as_ref() } else { None };

        (self.mapper)(before_ref, &mut current, after_ref);

        // Phase 4: Clone the result to return, then put current back in window for next slide
        let result = current.clone();
        self.window[1] = Some(current);
        Some(result)
    }
}

/// Extension trait to add `map_surrounding` method to iterators of `PileupMetrics`
pub trait PileupMetricsIteratorExt: Iterator<Item = PileupMetrics> + Sized {
    /// Map over elements while providing access to validated neighbors.
    ///
    /// The mapping function receives:
    /// - `before`: Reference to the previous element if it exists and is adjacent
    ///   (same contig, position - 1)
    /// - `current`: Mutable reference to the current element
    /// - `after`: Reference to the next element if it exists and is adjacent
    ///   (same contig, position + 1)
    ///
    /// The function should mutate `current` in place. The modified element is
    /// then yielded by the iterator.
    ///
    /// # Example
    /// ```ignore
    /// use rastair::utils::surrounding::PileupMetricsIteratorExt;
    ///
    /// let processed = pileups
    ///     .map_surrounding(|before, current, after| {
    ///         // Mutate current based on neighbors
    ///         if before.is_some() && after.is_some() {
    ///             current.has_both_neighbors = true;
    ///         }
    ///     })
    ///     .collect::<Vec<_>>();
    /// ```
    fn map_surrounding<F>(self, f: F) -> SurroundingMap<Self, F>
    where
        F: FnMut(Option<&PileupMetrics>, &mut PileupMetrics, Option<&PileupMetrics>),
    {
        SurroundingMap::new(self, f)
    }
}

/// Implement the extension trait for all iterators that yield `PileupMetrics`
impl<I> PileupMetricsIteratorExt for I where I: Iterator<Item = PileupMetrics> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        call::pileup::{Pileup, SimpleReads},
        sequence::{ChunkRegion, Region, Segment},
        vcf::SequenceContext,
    };
    use rastair_types::Base;
    use std::rc::Rc;

    /// Helper to create a minimal PileupMetrics for testing
    fn make_pileup(contig: &str, pos: u64) -> PileupMetrics {
        // Create a minimal segment with enough context
        let start = pos.saturating_sub(10);
        let end = pos + 20;
        let segment = Rc::new(Segment {
            range: ChunkRegion {
                region: Region { contig: contig.into(), start, end },
                last_position: end,
            },
            sequence: vec![b'C'; (end - start) as usize],
        });

        let pos_idx = (pos - start) as usize;
        let context = SequenceContext::new(pos_idx, &segment).expect("valid context");

        let pileup = Pileup {
            region: segment.range.clone(),
            context,
            pos: pos as u32,
            reads: SimpleReads(vec![]),
            reference_base: Base::C,
        };

        PileupMetrics::new(pileup).unwrap()
    }

    #[test]
    fn test_empty_iterator() {
        let items: Vec<PileupMetrics> = vec![];
        let mut calls = 0;

        let result: Vec<_> = items
            .into_iter()
            .map_surrounding(|before, _current, after| {
                calls += 1;
                assert!(before.is_none());
                assert!(after.is_none());
            })
            .collect();

        assert_eq!(result.len(), 0);
        assert_eq!(calls, 0);
    }

    #[test]
    fn test_single_element() {
        let items = vec![make_pileup("chr1", 100)];
        let mut calls = 0;

        let result: Vec<_> = items
            .into_iter()
            .map_surrounding(|before, current, after| {
                calls += 1;
                assert!(before.is_none(), "Single element should have no before");
                assert!(after.is_none(), "Single element should have no after");
                assert_eq!(current.pileup.pos, 100);
            })
            .collect();

        assert_eq!(result.len(), 1);
        assert_eq!(calls, 1);
        assert_eq!(result[0].pileup.pos, 100);
    }

    #[test]
    fn test_two_consecutive_elements() {
        let items = vec![make_pileup("chr1", 100), make_pileup("chr1", 101)];
        let mut first_call = true;

        let result: Vec<_> = items
            .into_iter()
            .map_surrounding(|before, current, after| {
                if first_call {
                    assert!(before.is_none());
                    assert!(after.is_some());
                    assert_eq!(after.unwrap().pileup.pos, 101);
                    assert_eq!(current.pileup.pos, 100);
                    first_call = false;
                } else {
                    assert!(before.is_some());
                    assert_eq!(before.unwrap().pileup.pos, 100);
                    assert!(after.is_none());
                    assert_eq!(current.pileup.pos, 101);
                }
            })
            .collect();

        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_three_consecutive_elements() {
        let items =
            vec![make_pileup("chr1", 100), make_pileup("chr1", 101), make_pileup("chr1", 102)];

        let result: Vec<_> = items
            .into_iter()
            .map_surrounding(|before, current, after| match current.pileup.pos {
                100 => {
                    assert!(before.is_none());
                    assert_eq!(after.unwrap().pileup.pos, 101);
                }
                101 => {
                    assert_eq!(before.unwrap().pileup.pos, 100);
                    assert_eq!(after.unwrap().pileup.pos, 102);
                }
                102 => {
                    assert_eq!(before.unwrap().pileup.pos, 101);
                    assert!(after.is_none());
                }
                _ => panic!("Unexpected position"),
            })
            .collect();

        assert_eq!(result.len(), 3);
    }

    #[test]
    fn test_non_consecutive_positions() {
        // Gap between positions means they shouldn't be treated as neighbors
        let items = vec![
            make_pileup("chr1", 100),
            make_pileup("chr1", 105), // Gap of 5
            make_pileup("chr1", 106),
        ];

        let result: Vec<_> = items
            .into_iter()
            .map_surrounding(|before, current, after| match current.pileup.pos {
                100 => {
                    assert!(before.is_none());
                    assert!(after.is_none(), "105 is not consecutive to 100");
                }
                105 => {
                    assert!(before.is_none(), "100 is not consecutive to 105");
                    assert_eq!(after.unwrap().pileup.pos, 106);
                }
                106 => {
                    assert_eq!(before.unwrap().pileup.pos, 105);
                    assert!(after.is_none());
                }
                _ => panic!("Unexpected position"),
            })
            .collect();

        assert_eq!(result.len(), 3);
    }

    #[test]
    fn test_different_contigs() {
        let items = vec![
            make_pileup("chr1", 100),
            make_pileup("chr2", 101), // Different contig
            make_pileup("chr2", 102),
        ];

        let result: Vec<_> = items
            .into_iter()
            .map_surrounding(|before, current, after| {
                match (current.pileup.contig().as_str(), current.pileup.pos) {
                    ("chr1", 100) => {
                        assert!(before.is_none());
                        assert!(after.is_none(), "chr2 is different contig");
                    }
                    ("chr2", 101) => {
                        assert!(before.is_none(), "chr1 is different contig");
                        assert_eq!(after.unwrap().pileup.pos, 102);
                    }
                    ("chr2", 102) => {
                        assert_eq!(before.unwrap().pileup.pos, 101);
                        assert!(after.is_none());
                    }
                    _ => panic!("Unexpected contig/position"),
                }
            })
            .collect();

        assert_eq!(result.len(), 3);
    }

    #[test]
    fn test_mutation_is_preserved() {
        let items =
            vec![make_pileup("chr1", 100), make_pileup("chr1", 101), make_pileup("chr1", 102)];

        let result: Vec<_> = items
            .into_iter()
            .map_surrounding(|_before, current, _after| {
                // Mutate the position by adding 1000
                current.pileup.pos += 1000;
            })
            .collect();

        assert_eq!(result.len(), 3);
        assert_eq!(result[0].pileup.pos, 1100);
        assert_eq!(result[1].pileup.pos, 1101);
        assert_eq!(result[2].pileup.pos, 1102);
    }

    #[test]
    fn test_long_sequence() {
        // Test with more elements to ensure sliding window works correctly
        let items: Vec<_> = (0..10).map(|i| make_pileup("chr1", 100 + i)).collect();

        let result: Vec<_> = items
            .into_iter()
            .map_surrounding(|before, current, after| {
                let pos = current.pileup.pos;
                if pos == 100 {
                    assert!(before.is_none());
                    assert!(after.is_some());
                } else if pos == 109 {
                    assert!(before.is_some());
                    assert!(after.is_none());
                } else {
                    assert!(before.is_some());
                    assert!(after.is_some());
                    assert_eq!(before.unwrap().pileup.pos, pos - 1);
                    assert_eq!(after.unwrap().pileup.pos, pos + 1);
                }
            })
            .collect();

        assert_eq!(result.len(), 10);
    }

    #[test]
    fn test_mixed_gaps_and_contigs() {
        let items = vec![
            make_pileup("chr1", 100),
            make_pileup("chr1", 101),
            make_pileup("chr1", 105), // Gap
            make_pileup("chr2", 106), // Different contig (but consecutive pos)
            make_pileup("chr2", 107),
        ];

        let result: Vec<_> = items
            .into_iter()
            .map_surrounding(|before, current, after| {
                match (current.pileup.contig().as_str(), current.pileup.pos) {
                    ("chr1", 100) => {
                        assert!(before.is_none());
                        assert_eq!(after.unwrap().pileup.pos, 101);
                    }
                    ("chr1", 101) => {
                        assert_eq!(before.unwrap().pileup.pos, 100);
                        assert!(after.is_none(), "Gap to 105");
                    }
                    ("chr1", 105) => {
                        assert!(before.is_none(), "Gap from 101");
                        assert!(after.is_none(), "Different contig");
                    }
                    ("chr2", 106) => {
                        assert!(before.is_none(), "Different contig from chr1");
                        assert_eq!(after.unwrap().pileup.pos, 107);
                    }
                    ("chr2", 107) => {
                        assert_eq!(before.unwrap().pileup.pos, 106);
                        assert!(after.is_none());
                    }
                    _ => panic!("Unexpected combination"),
                }
            })
            .collect();

        assert_eq!(result.len(), 5);
    }
}
