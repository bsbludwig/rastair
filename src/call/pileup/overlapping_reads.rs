// TODO: Replace this with a more efficient dedupe algorithm
//  - track qnames separately and not in SimpleRead
//  - use a FxHashMap for quick lookup by qname suffix
// TODO: Add tests for this module

use super::SimpleReads;
use smallvec::SmallVec;

impl SimpleReads {
    /// Remove overlapping reads from the same fragment.
    pub fn remove_overlapping_pairs(&mut self) {
        // For each read, check if we already saw one with the same name.
        //
        // If the bases agree, keep only the first one. If they disagree, keep none.
        //
        // But this is rust -- so we can't just remove elements while iterating.
        // Instead, we keep a little list of indices to remove, and then remove them afterwards.
        // This should be fine since the amount of items to remove is typically small.
        let mut to_remove = SmallVec::<usize, 16>::new();
        for i in 0..self.0.len() {
            let base_i = &self.0[i];
            for j in (i + 1)..self.0.len() {
                let base_j = &self.0[j];
                if base_i.qname == base_j.qname {
                    // Same read name
                    if base_i.base == base_j.base {
                        // Same base, keep only the first one
                        to_remove.push(j);
                    } else {
                        // Different bases, ignore the second in pair
                        // NOTE: This is different from rastair1
                        if base_i.second {
                            to_remove.push(i);
                        } else {
                            to_remove.push(j);
                        }
                    }
                    // No need to check further
                    break;
                }
            }
        }
        // Remove duplicates
        to_remove.sort_unstable();
        for &idx in to_remove.iter().rev() {
            self.0.swap_remove(idx);
        }
    }
}
