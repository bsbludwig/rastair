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

            let end = self.current_start.saturating_add(self.max_length).min(full_region.end);
            let chunk = ChunkRegion {
                region: Region {
                    contig: full_region.contig.clone(),
                    start: self.current_start,
                    end,
                },
                last_position: full_region.end,
            };

            self.current_start = end;
            if self.current_start < full_region.end {
                self.current_start = self.current_start.saturating_sub(self.overlap);
            }

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
            start in 1u64..1000u64,
            length in 100u64..10000u64,
            max_length in 300u64..1000u64,
            overlap in 10u64..100u64
        ) {
            let end = start + length;
            let full_regions = vec![SelectedRegion::EntireContig(Region {
                contig: "chr1".into(),
                start,
                end,
            })];

            let chunked = ChunkedRegions {
                full_regions,
                current_region_idx: 0,
                current_start: start,
                max_length,
                overlap,
            };

            // Collect all chunks
            let chunks: Vec<_> = chunked.collect();

            // Property 1: There should be at least one chunk
            prop_assert!(!chunks.is_empty());

            // Property 2: First chunk should start at the beginning
            prop_assert_eq!(chunks[0].region.start, start);

            // Property 3: Last chunk should end at the end
            prop_assert_eq!(chunks.last().unwrap().region.end, end);

            // Property 4: No chunk should exceed max_length
            for chunk in &chunks {
                prop_assert!(chunk.region.end - chunk.region.start <= max_length);
            }

            // Property 5: Adjacent chunks should overlap by the specified amount (except possibly the last one)
            for pair in chunks.windows(2) {
                let current = &pair[0];
                let next = &pair[1];

                // The overlap between chunks
                let actual_overlap = current.region.end - next.region.start;

                // If this isn't the last chunk leading to the end, it should have the specified overlap
                if next.region.end < end {
                    prop_assert_eq!(actual_overlap, overlap);
                }
            }

            // Property 6: No gaps between chunks
            let mut covered = HashSet::new();
            for chunk in &chunks {
                for pos in chunk.region.start..chunk.region.end {
                    covered.insert(pos);
                }
            }
            for pos in start..end {
                prop_assert!(covered.contains(&pos), "Position {} not covered by any chunk", pos);
            }

            // Property 7: All chunks should have the correct last_position
            for chunk in &chunks {
                prop_assert_eq!(chunk.last_position, end);
            }
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
