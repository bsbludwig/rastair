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
            if self.current_region_idx >= self.full_regions.len() {
                return None;
            }

            let full_region = &self.full_regions[self.current_region_idx];
            if self.current_start >= full_region.end {
                // Move to next region when we've finished the current one
                self.current_region_idx += 1;
                if self.current_region_idx < self.full_regions.len() {
                    self.current_start = self.full_regions[self.current_region_idx].start;
                }
                continue;
            }

            let logical_start = self.current_start;
            let logical_end =
                self.current_start.saturating_add(self.max_length).min(full_region.end);

            let last_position = match full_region {
                SelectedRegion::EntireContig(region) => region.end,
                SelectedRegion::UserDefinedSubset { last_position, .. } => *last_position,
            };
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
    use proptest::prelude::*;
    use std::collections::HashSet;

    proptest! {
        #[test]
        fn chunks_cover_full_region(
            // Generate reasonable chromosome lengths and chunk parameters
            start in 0u64..1000u64,
            length in 100u64..10000u64,
            tail in 0u64..500u64,
            max_length in 300u64..1000u64,
            overlap in 10u64..100u64
        ) {
            let end = start + length;
            let last_position = end + tail;
            let full_regions = vec![SelectedRegion::UserDefinedSubset {
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
            prop_assert!(!chunks.is_empty());

            // The cores tile the region exactly: each position is one
            // chunk's, in order, with no gap.
            let mut next_core = start;
            for chunk in &chunks {
                let core_start = chunk.region.start + chunk.overlap_start;
                let core_end = chunk.region.end - chunk.overlap_end;
                prop_assert_eq!(core_start, next_core);
                prop_assert!(core_end > core_start);
                next_core = core_end;
            }
            prop_assert_eq!(next_core, end);

            for chunk in &chunks {
                // The overlap reaches past the region on both sides, as far as
                // the contig and `overlap` allow.
                prop_assert_eq!(chunk.overlap_start, overlap.min(chunk.region.start + chunk.overlap_start));
                prop_assert!(chunk.overlap_end <= overlap);
                prop_assert!(chunk.region.end <= last_position);
                prop_assert!(chunk.region.end - chunk.region.start <= max_length + 2 * overlap);
                prop_assert_eq!(chunk.last_position, last_position);
            }
            let last = chunks.last().unwrap();
            prop_assert_eq!(last.overlap_end, overlap.min(last_position - end));
        }

        #[test]
        fn handles_multiple_regions(
            region1_length in 100u64..1000u64,
            region2_length in 100u64..1000u64,
            max_length in 10u64..100u64,
            overlap in 1u64..10u64
        ) {
            let full_regions = vec![
                SelectedRegion::EntireContig(Region {
                    contig: "chr1".into(),
                    start: 1,
                    end: region1_length,
                }),
                SelectedRegion::EntireContig(Region {
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
            prop_assert!(!chr1_chunks.is_empty());
            prop_assert!(!chr2_chunks.is_empty());

            // Property 2: Chunks should be in order by chromosome
            for pair in chunks.windows(2) {
                if pair[0].region.contig == pair[1].region.contig {
                    prop_assert!(pair[0].region.start <= pair[1].region.start);
                }
            }

            // Property 3: Each region's chunks should cover their full length
            let mut chr1_covered = HashSet::new();
            let mut chr2_covered = HashSet::new();

            for chunk in &chunks {
                let covered = if chunk.region.contig == "chr1" {
                    &mut chr1_covered
                } else {
                    &mut chr2_covered
                };
                for pos in chunk.region.start..chunk.region.end {
                    covered.insert(pos);
                }
            }

            for pos in 1..region1_length {
                prop_assert!(chr1_covered.contains(&pos));
            }
            for pos in 1..region2_length {
                prop_assert!(chr2_covered.contains(&pos));
            }
        }
    }
}
