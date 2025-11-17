use rustc_hash::{FxBuildHasher, FxHashSet};
use smallvec::SmallVec;

/// Efficient read name deduplicator optimized for genomic sequencing data.
///
/// Uses a two-stage approach:
/// 1. Fast suffix-based filtering using the last 8 bytes as u64
/// 2. Full string comparison only when suffix collisions occur (rare)
///
/// This design exploits the fact that sequencing read names often have common
/// prefixes but unique suffixes, making suffix collisions rare.
pub struct ReadDeduplicator {
    /// Set of suffixes (last 8 bytes) for fast initial filtering
    suffixes: FxHashSet<u64>,
    /// Full read names, compared only when suffixes collide
    full_strings: Vec<SmallVec<u8, 48>>,
}

impl ReadDeduplicator {
    /// Create a new read deduplicator with estimated capacity
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            suffixes: FxHashSet::with_capacity_and_hasher(capacity, FxBuildHasher),
            full_strings: Vec::with_capacity(capacity),
        }
    }

    /// Extract the last 8 bytes of a read name as a u64 for fast comparison
    #[inline]
    fn extract_suffix(bytes: &[u8]) -> u64 {
        let len = bytes.len();

        if len >= 8 {
            // Take last 8 bytes and convert to u64 (little-endian)
            let start = len - 8;
            u64::from_le_bytes([
                bytes[start],
                bytes[start + 1],
                bytes[start + 2],
                bytes[start + 3],
                bytes[start + 4],
                bytes[start + 5],
                bytes[start + 6],
                bytes[start + 7],
            ])
        } else {
            // For short strings, pad with zeros
            let mut padded = [0u8; 8];
            padded[..len].copy_from_slice(bytes);
            u64::from_le_bytes(padded)
        }
    }

    /// Check if a read name is a duplicate. Returns `true` if duplicate, `false` if new.
    ///
    /// Inserts the read name into the deduplicator if it's new.
    ///
    /// # Performance
    /// - O(1) average case when no suffix collision
    /// - O(n) worst case where n is number of previous suffix collisions (typically very small)
    pub fn is_duplicate(&mut self, read_name: &[u8]) -> bool {
        let suffix = Self::extract_suffix(read_name);

        // Try to insert the suffix - if successful, it's definitely new
        if self.suffixes.insert(suffix) {
            // New suffix, definitely not a duplicate, add the full string
            self.full_strings.push(SmallVec::from(read_name));
            false
        } else {
            // Already saw this suffix -- make sure it's a real duplicate
            //
            // Assumption: This list is small (<100 items) so linear search is
            // acceptable and faster than inserting into a `HashSet`
            if self.full_strings.iter().any(|stored_name| stored_name.as_slice() == read_name) {
                // Found exact duplicate
                true
            } else {
                // No exact match found, but suffix collision occurred
                // Store the full read name for future comparisons
                self.full_strings.push(SmallVec::from(read_name));
                false
            }
        }
    }

    /// Total number of unique read names processed
    pub fn len(&self) -> usize {
        self.suffixes.len()
    }

    /// Check if no reads have been processed
    pub fn is_empty(&self) -> bool {
        self.suffixes.is_empty()
    }

    /// Clear all stored data
    pub fn clear(&mut self) {
        self.suffixes.clear();
        self.full_strings.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_deduplication() {
        let mut dedup = ReadDeduplicator::with_capacity(5);

        // First occurrence should not be duplicate
        assert!(!dedup.is_duplicate(b"LH00407:74:22NV2TLT3:1:1116:38058:14010"));
        assert_eq!(dedup.len(), 1);

        // Different read should not be duplicate
        assert!(!dedup.is_duplicate(b"LH00407:74:22NV2TLT3:3:2186:24726:25065"));
        assert_eq!(dedup.len(), 2);

        // Same read should be duplicate
        assert!(dedup.is_duplicate(b"LH00407:74:22NV2TLT3:1:1116:38058:14010"));
        assert_eq!(dedup.len(), 2); // Length shouldn't change
    }

    #[test]
    fn test_suffix_extraction() {
        let read1 = b"LH00407:74:22NV2TLT3:1:1116:38058:14010";
        let read2 = b"LH00407:74:22NV2TLT3:3:2186:24726:25065";

        let suffix1 = ReadDeduplicator::extract_suffix(read1);
        let suffix2 = ReadDeduplicator::extract_suffix(read2);

        // Should be different suffixes
        assert_ne!(suffix1, suffix2);
    }

    #[test]
    fn test_short_strings() {
        let mut dedup = ReadDeduplicator::with_capacity(5);

        assert!(!dedup.is_duplicate(b"short"));
        assert!(!dedup.is_duplicate(b"abc"));
        assert!(dedup.is_duplicate(b"short"));
        assert!(dedup.is_duplicate(b"abc"));
    }

    #[test]
    fn test_collision_handling() {
        let mut dedup = ReadDeduplicator::with_capacity(5);

        // Create strings that might have same suffix but different content
        let read1 = b"prefix1:12345678";
        let read2 = b"prefix2:12345678"; // Same suffix, different prefix

        assert!(!dedup.is_duplicate(read1));
        assert!(!dedup.is_duplicate(read2));

        // Now test actual duplicates
        assert!(dedup.is_duplicate(read1));
        assert!(dedup.is_duplicate(read2));
    }
}
