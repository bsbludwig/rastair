use crate::sequence::Segment;
use color_eyre::Result;

const HALF_WINDOW: usize = 50;

pub struct SlidingEntropy<'a> {
    sequence: &'a [u8],
    // Order: [C, T, A, G] — matches the original Counter::entries() iteration order to preserve
    // floating-point summation order and keep entropy values numerically identical to the old code.
    counts: [usize; 4],
    unknown_count: usize,
    window_start: usize,
    window_end: usize,
    initialized: bool,
}

impl<'a> SlidingEntropy<'a> {
    pub fn new(segment: &'a Segment) -> Self {
        Self {
            sequence: &segment.sequence,
            counts: [0; 4],
            unknown_count: 0,
            window_start: 0,
            window_end: 0,
            initialized: false,
        }
    }

    pub fn entropy_at(&mut self, idx: usize) -> f64 {
        let new_start = idx.saturating_sub(HALF_WINDOW).min(self.sequence.len());
        let new_end = idx.saturating_add(HALF_WINDOW + 1).min(self.sequence.len());

        if !self.initialized || new_start >= self.window_end || new_end <= self.window_start {
            self.initialize(new_start, new_end);
        } else {
            self.slide_to(new_start, new_end);
        }

        self.compute_entropy()
    }

    fn initialize(&mut self, start: usize, end: usize) {
        self.counts = [0; 4];
        self.unknown_count = 0;
        for &byte in &self.sequence[start..end] {
            self.add(byte);
        }
        self.window_start = start;
        self.window_end = end;
        self.initialized = true;
    }

    fn slide_to(&mut self, new_start: usize, new_end: usize) {
        // Remove bases that left the window on the left (window moved right)
        if new_start > self.window_start {
            for &byte in &self.sequence[self.window_start..new_start] {
                self.remove(byte);
            }
        }
        // Remove bases that left the window on the right
        if new_end < self.window_end {
            for &byte in &self.sequence[new_end..self.window_end] {
                self.remove(byte);
            }
        }
        // Add bases that entered the window on the right
        if new_end > self.window_end {
            for &byte in &self.sequence[self.window_end..new_end] {
                self.add(byte);
            }
        }
        // Add bases that entered the window on the left
        if new_start < self.window_start {
            for &byte in &self.sequence[new_start..self.window_start] {
                self.add(byte);
            }
        }
        self.window_start = new_start;
        self.window_end = new_end;
    }

    fn add(&mut self, byte: u8) {
        match byte {
            b'C' | b'c' => self.counts[0] += 1,
            b'T' | b't' => self.counts[1] += 1,
            b'A' | b'a' => self.counts[2] += 1,
            b'G' | b'g' => self.counts[3] += 1,
            _ => self.unknown_count += 1,
        }
    }

    fn remove(&mut self, byte: u8) {
        match byte {
            b'C' | b'c' => {
                debug_assert!(self.counts[0] > 0, "underflow removing C");
                self.counts[0] = self.counts[0].saturating_sub(1);
            }
            b'T' | b't' => {
                debug_assert!(self.counts[1] > 0, "underflow removing T");
                self.counts[1] = self.counts[1].saturating_sub(1);
            }
            b'A' | b'a' => {
                debug_assert!(self.counts[2] > 0, "underflow removing A");
                self.counts[2] = self.counts[2].saturating_sub(1);
            }
            b'G' | b'g' => {
                debug_assert!(self.counts[3] > 0, "underflow removing G");
                self.counts[3] = self.counts[3].saturating_sub(1);
            }
            _ => {
                debug_assert!(self.unknown_count > 0, "underflow removing unknown");
                self.unknown_count = self.unknown_count.saturating_sub(1);
            }
        }
    }

    fn compute_entropy(&self) -> f64 {
        let total = (self.counts.iter().sum::<usize>() + self.unknown_count) as f64;
        if total == 0.0 {
            return 0.0;
        }

        self.counts
            .iter()
            .filter(|&&count| count > 0)
            .map(|&count| {
                let p = (count as f64) / total;
                -p * p.log2()
            })
            .sum::<f64>()
    }
}

impl Segment {
    pub fn entropy_around<const N: usize>(&self, idx: usize) -> Result<f64> {
        let seq_context = self.get(idx.saturating_sub(N / 2), idx.saturating_add(N / 2 + 1))?;

        Ok(entropy(seq_context))
    }
}

/// Calculate Shannon entropy for sequence context
pub fn entropy(sequence: &[u8]) -> f64 {
    // counts order: [C, T, A, G] — matches the original Counter::entries() order to preserve
    // floating-point summation stability across refactors.
    let mut counts = [0usize; 4];
    let mut unknown = 0usize;
    for &byte in sequence {
        match byte {
            b'C' | b'c' => counts[0] += 1,
            b'T' | b't' => counts[1] += 1,
            b'A' | b'a' => counts[2] += 1,
            b'G' | b'g' => counts[3] += 1,
            _ => unknown += 1,
        }
    }

    let total = (counts.iter().sum::<usize>() + unknown) as f64;
    if total == 0.0 {
        return 0.0;
    }

    counts
        .iter()
        .filter(|&&count| count > 0)
        .map(|&count| {
            let p = (count as f64) / total;
            -p * p.log2()
        })
        .sum::<f64>()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::{ChunkRegion, Region, Segment};
    use proptest::{prop_assert, proptest};
    use seqair_types::SmolStr;
    use std::iter::repeat_n;

    fn test_segment(sequence: Vec<u8>) -> Segment {
        let len = sequence.len();
        Segment {
            range: ChunkRegion {
                region: Region { contig: SmolStr::new("chr13"), start: 1, end: len as u64 },
                last_position: len as u64,
                overlap_start: 0,
                overlap_end: 0,
            },
            sequence,
            overlap_start: 0,
            overlap_end: 0,
        }
    }

    fn assert_sliding_matches(segment: &Segment, positions: &[usize]) {
        let mut sliding = SlidingEntropy::new(segment);
        for &idx in positions {
            let standalone = segment.entropy_around::<100>(idx).unwrap();
            let incremental = sliding.entropy_at(idx);
            assert!(
                (standalone - incremental).abs() < 1e-12,
                "Mismatch at idx {idx}: standalone={standalone}, incremental={incremental}"
            );
        }
    }

    #[test]
    fn low_entropy() {
        let segment = test_segment(repeat_n(b'A', 100).collect());

        let entropy = segment.entropy_around::<100>(50).unwrap();

        assert_eq!(entropy, 0.0);
    }

    #[test]
    fn high_entropy() {
        let sequence = repeat_n(b"ACTG", 25).flat_map(|x| *x).collect();
        let segment = test_segment(sequence);

        let entropy = segment.entropy_around::<100>(50).unwrap();

        assert!(entropy > 0.0);
    }

    proptest! {
        #[test]
        fn entropy_is_non_negative(sequence in "[ACTG]{100}") {
            let segment = test_segment(sequence.into_bytes());

            let entropy = segment.entropy_around::<100>(50).unwrap();

            // Shannon entropy over an alphabet of 4 bases is bounded to [0, log2(4)] = [0, 2]
            prop_assert!(entropy >= 0.0);
            prop_assert!(entropy <= 2.0 + f64::EPSILON);
        }
    }

    #[test]
    fn sliding_matches_standalone_full_forward_scan() {
        let sequence: Vec<u8> = (0..500).map(|i| b"ACGTAACCGGTT"[i % 12]).collect();
        let segment = test_segment(sequence);
        assert_sliding_matches(&segment, &(0..500).collect::<Vec<_>>());
    }

    #[test]
    fn sliding_with_sparse_forward_positions() {
        let sequence: Vec<u8> = (0..1000).map(|i| b"ACGTAACCGGTT"[i % 12]).collect();
        let segment = test_segment(sequence);
        assert_sliding_matches(&segment, &[10, 50, 200, 201, 202, 500, 999]);
    }

    #[test]
    fn sliding_edge_position_zero() {
        let sequence: Vec<u8> = (0..200).map(|i| b"ACGTAACCGGTT"[i % 12]).collect();
        let segment = test_segment(sequence);
        assert_sliding_matches(&segment, &[0, 1, 2, 50, 100]);
    }

    #[test]
    fn sliding_edge_position_end() {
        let sequence: Vec<u8> = (0..200).map(|i| b"ACGTAACCGGTT"[i % 12]).collect();
        let len = sequence.len();
        let segment = test_segment(sequence);
        assert_sliding_matches(&segment, &[len - 1, len - 2, len - 51, len - 52]);
    }

    #[test]
    fn sliding_short_segment_smaller_than_window() {
        // Sequence shorter than the full 101-base window — clamping must work correctly.
        let sequence: Vec<u8> = b"ACGTACGT".to_vec();
        let segment = test_segment(sequence);
        assert_sliding_matches(&segment, &[0, 1, 4, 7]);
    }

    #[test]
    fn sliding_with_n_bases() {
        // N bases should be counted as unknown (in denominator, not entropy sum).
        let sequence: Vec<u8> = b"ACGTNNNACGT".to_vec();
        let segment = test_segment(sequence);
        assert_sliding_matches(&segment, &[0, 5, 10]);
    }

    #[test]
    fn sliding_backward_jump() {
        // Backward access triggers reinit; result must still match standalone.
        let sequence: Vec<u8> = (0..500).map(|i| b"ACGTAACCGGTT"[i % 12]).collect();
        let segment = test_segment(sequence);
        assert_sliding_matches(&segment, &[300, 400, 100, 200, 50]);
    }

    #[test]
    fn sliding_same_position_twice() {
        let sequence: Vec<u8> = (0..200).map(|i| b"ACGTAACCGGTT"[i % 12]).collect();
        let segment = test_segment(sequence);
        assert_sliding_matches(&segment, &[100, 100, 101, 100]);
    }
}
