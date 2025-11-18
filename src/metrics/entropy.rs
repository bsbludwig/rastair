use crate::{
    sequence::Segment,
    utils::{Base, Counter},
};
use color_eyre::Result;

impl Segment {
    pub fn entropy_around<const N: usize>(&self, idx: usize) -> Result<f64> {
        let seq_context = self.get(idx.saturating_sub(N / 2), idx.saturating_add(N / 2 + 1))?;

        Ok(entropy(seq_context))
    }
}

///  Calculate Shannon entropy for sequence context
pub fn entropy(sequence: &[u8]) -> f64 {
    let counts: Counter = sequence.iter().map(|&b| Base::from(b)).collect();
    let total = sequence.len() as f64;

    counts
        .entries()
        .iter()
        .filter(|(_base, count)| *count > 0)
        .map(|(_base, count)| {
            let p = (*count as f64) / total;
            -p * p.log2()
        })
        .sum::<f64>()
}

#[cfg(test)]
mod tests {
    use crate::sequence::{ChunkRegion, Region, Segment};
    use proptest::proptest;
    use rastair_types::SmolStr;
    use std::iter::repeat_n;

    #[test]
    fn low_entropy() {
        let region = ChunkRegion {
            region: Region { contig: SmolStr::new("chr13"), start: 1, end: 100 },
            last_position: 100,
        };
        let segment = Segment { range: region, sequence: repeat_n(b'A', 100).collect() };
        let pos = 50;

        let entropy = segment.entropy_around::<100>(pos).unwrap();

        assert_eq!(entropy, 0.0);
    }

    #[test]
    fn high_entropy() {
        let region = ChunkRegion {
            region: Region { contig: SmolStr::new("chr13"), start: 1, end: 100 },
            last_position: 100,
        };
        let sequence = repeat_n(b"ACTG", 25).flat_map(|x| *x).collect();
        let segment = Segment { range: region, sequence };
        let pos = 50;

        let entropy = segment.entropy_around::<100>(pos).unwrap();

        assert!(entropy > 0.0);
    }

    proptest! {
        #[test]
        fn entropy_is_non_negative(sequence in "[ACTG]{100}") {
            let region = ChunkRegion {
                region: Region { contig: SmolStr::new("chr13"), start: 1, end: 100 },
                last_position: 100,
            };
            let segment = Segment { range: region, sequence: sequence.into_bytes(),
            };
            let pos = 50;

            let entropy = segment.entropy_around::<100>(pos).unwrap();

            assert!(entropy >= 0.0);
        }
    }
}
