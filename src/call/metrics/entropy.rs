use crate::{
    call::{variants::VariantCandidatePileup, vcf::Entropy},
    utils::{Base, Counter},
};
use smallvec::smallvec;

impl VariantCandidatePileup {
    ///  Calculate Shannon entropy for 100bp context around variant position
    pub(crate) fn entropy(&self) -> Entropy {
        let idx = self.idx();
        let seq_context = self
            .segment
            .get(idx.saturating_sub(50), idx.saturating_add(51))
            .expect("sequence context indices are valid");

        let counts: Counter = seq_context.iter().filter_map(|&b| Base::try_from(b).ok()).collect();
        let total = seq_context.len() as f64;
        let entropy = counts
            .entries()
            .iter()
            .filter(|(_base, count)| *count > 0)
            .map(|(_base, count)| {
                let p = (*count as f64) / total;
                -p * p.log2()
            })
            .sum::<f64>();

        Entropy(smallvec![entropy])
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        call::variants::SeenBases,
        sequence::{ChunkRegion, Region, Segment},
    };
    use proptest::proptest;
    use smol_str::SmolStr;
    use std::{iter::repeat_n, rc::Rc};

    use super::*;

    #[test]
    fn low_entropy() {
        let segment = Segment {
            range: ChunkRegion {
                region: Region { chromosome: SmolStr::new("chr13"), start: 1, end: 100 },
                last_position: 100,
            },
            sequence: repeat_n(b'A', 100).collect(),
        };
        let pileup = VariantCandidatePileup {
            segment: Rc::new(segment),
            pos: 50,
            bases: SeenBases(smallvec![]),
            reference_base: Base::A,
        };
        let entropy = pileup.entropy();
        assert_eq!(entropy.0[0], 0.0);
    }

    #[test]
    fn high_entropy() {
        let segment = Segment {
            range: ChunkRegion {
                region: Region { chromosome: SmolStr::new("chr13"), start: 1, end: 100 },
                last_position: 100,
            },
            sequence: repeat_n(b"ACTG", 25).flat_map(|x| *x).collect(),
        };
        let pileup = VariantCandidatePileup {
            segment: Rc::new(segment),
            pos: 50,
            bases: SeenBases(smallvec![]),
            reference_base: Base::A,
        };
        let entropy = pileup.entropy();
        assert!(entropy.0[0] > 0.0);
    }

    proptest! {
        #[test]
        fn entropy_is_non_negative(sequence in "[ACTG]{100}") {
            let segment = Segment {
                range: ChunkRegion {
                    region: Region { chromosome: SmolStr::new("chr13"), start: 1, end: 100 },
                    last_position: 100,
                },
                sequence: sequence.into_bytes(),
            };
            let pileup = VariantCandidatePileup {
                segment: Rc::new(segment),
                pos: 50,
                bases: SeenBases(smallvec![]),
                reference_base: Base::A,
            };
            let entropy = pileup.entropy();
            assert!(entropy.0[0] >= 0.0);
        }
    }
}
