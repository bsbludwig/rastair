use super::ChunkRegion;
use super::Region;
use super::SelectedRegion;

pub(crate) struct ChunkedRegions {
    pub(crate) full_regions: Vec<SelectedRegion>,
    pub(crate) current_region_idx: usize,
    pub(crate) current_start: u64,
    pub(crate) max_length: u64,
    pub(crate) overlap: u64,
}

impl Iterator for ChunkedRegions {
    type Item = ChunkRegion;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let full_region = self.full_regions.get(self.current_region_idx)?;
            if self.current_start >= full_region.end {
                // Move to next region when we've finished the current one
                self.current_region_idx += 1;
                if let Some(next) = self.full_regions.get(self.current_region_idx) {
                    self.current_start = next.start;
                }
                continue;
            }

            let logical_start = self.current_start;
            let logical_end =
                self.current_start.saturating_add(self.max_length).min(full_region.end);

            let last_position = full_region.last_position;
            // The overlap reaches past the region's own ends, as far as the
            // contig goes, so the bases on a requested region's edges have the
            // same context as those on a chunk boundary inside it.
            let chunk_start = self.current_start.saturating_sub(self.overlap);
            let chunk_end = logical_end.saturating_add(self.overlap).min(last_position);

            let overlap_start = logical_start - chunk_start;
            let overlap_end = chunk_end - logical_end;

            let chunk = ChunkRegion {
                region: Region {
                    contig: full_region.contig.clone(),
                    start: chunk_start,
                    end: chunk_end,
                },
                last_position,
                overlap_start,
                overlap_end,
            };

            self.current_start = self.current_start.saturating_add(self.max_length);

            return Some(chunk);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hegel::{TestCase, generators as gs};
    use std::collections::HashSet;

    #[hegel::test]
    fn chunks_cover_full_region(tc: TestCase) {
        let start = tc.draw(gs::integers::<u64>().max_value(999));
        let length = tc.draw(gs::integers::<u64>().min_value(100).max_value(9_999));
        let tail = tc.draw(gs::integers::<u64>().max_value(499));
        let max_length = tc.draw(gs::integers::<u64>().min_value(300).max_value(999));
        let overlap = tc.draw(gs::integers::<u64>().min_value(10).max_value(99));
        let end = start + length;
        let last_position = end + tail;
        let full_regions = vec![SelectedRegion {
            region: Region { contig: "chr1".into(), start, end },
            last_position,
        }];

        let chunked = ChunkedRegions {
            full_regions,
            current_region_idx: 0,
            current_start: start,
            max_length,
            overlap,
        };

        let chunks: Vec<_> = chunked.collect();
        assert!(!chunks.is_empty());

        // The cores tile the region exactly: each position is one
        // chunk's, in order, with no gap.
        let mut next_core = start;
        for chunk in &chunks {
            let core_start = chunk.region.start + chunk.overlap_start;
            let core_end = chunk.region.end - chunk.overlap_end;
            assert_eq!(core_start, next_core);
            assert!(core_end > core_start);
            next_core = core_end;
        }
        assert_eq!(next_core, end);

        for chunk in &chunks {
            // The overlap reaches past the region on both sides, as far as
            // the contig and `overlap` allow.
            assert_eq!(chunk.overlap_start, overlap.min(chunk.region.start + chunk.overlap_start));
            assert!(chunk.overlap_end <= overlap);
            assert!(chunk.region.end <= last_position);
            assert!(chunk.region.end - chunk.region.start <= max_length + 2 * overlap);
            assert_eq!(chunk.last_position, last_position);
        }
        let last = chunks.last().unwrap();
        assert_eq!(last.overlap_end, overlap.min(last_position - end));
    }

    #[hegel::test]
    fn handles_multiple_regions(tc: TestCase) {
        let region1_length = tc.draw(gs::integers::<u64>().min_value(100).max_value(999));
        let region2_length = tc.draw(gs::integers::<u64>().min_value(100).max_value(999));
        let max_length = tc.draw(gs::integers::<u64>().min_value(10).max_value(99));
        let overlap = tc.draw(gs::integers::<u64>().min_value(1).max_value(9));
        let full_regions = vec![
            SelectedRegion::entire_contig(Region {
                contig: "chr1".into(),
                start: 1,
                end: region1_length,
            }),
            SelectedRegion::entire_contig(Region {
                contig: "chr2".into(),
                start: 1,
                end: region2_length,
            }),
        ];

        let chunked = ChunkedRegions {
            full_regions,
            current_region_idx: 0,
            current_start: 1,
            max_length,
            overlap,
        };

        let chunks: Vec<_> = chunked.collect();

        // Property 1: Should have chunks from both regions
        let chr1_chunks: Vec<_> = chunks.iter().filter(|c| c.region.contig == "chr1").collect();
        let chr2_chunks: Vec<_> = chunks.iter().filter(|c| c.region.contig == "chr2").collect();
        assert!(!chr1_chunks.is_empty());
        assert!(!chr2_chunks.is_empty());

        // Property 2: Chunks should be in order by chromosome
        for pair in chunks.windows(2) {
            if pair[0].region.contig == pair[1].region.contig {
                assert!(pair[0].region.start <= pair[1].region.start);
            }
        }

        // Property 3: Each region's chunks should cover their full length
        let mut chr1_covered = HashSet::new();
        let mut chr2_covered = HashSet::new();

        for chunk in &chunks {
            let covered =
                if chunk.region.contig == "chr1" { &mut chr1_covered } else { &mut chr2_covered };
            for pos in chunk.region.start..chunk.region.end {
                covered.insert(pos);
            }
        }

        for pos in 1..region1_length {
            assert!(chr1_covered.contains(&pos));
        }
        for pos in 1..region2_length {
            assert!(chr2_covered.contains(&pos));
        }
    }
}
