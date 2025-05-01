use super::ChunkRegion;
use super::FullRegion;
use super::Region;

pub(crate) struct ChunkedRegions {
    pub(crate) full_regions: Vec<FullRegion>,
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
            if self.current_start >= full_region.0.end {
                // Move to next region when we've finished the current one
                self.current_region_idx += 1;
                if self.current_region_idx < self.full_regions.len() {
                    self.current_start = self.full_regions[self.current_region_idx].0.start;
                }
                continue;
            }

            let end = self.current_start.saturating_add(self.max_length).min(full_region.0.end);
            let chunk = ChunkRegion {
                region: Region {
                    chromosome: full_region.0.chromosome.clone(),
                    start: self.current_start,
                    end,
                },
                last_position: full_region.0.end,
            };

            self.current_start = end;
            if self.current_start < full_region.0.end {
                self.current_start = self.current_start.saturating_sub(self.overlap);
            }

            return Some(chunk);
        }
    }
}
