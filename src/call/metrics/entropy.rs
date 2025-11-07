use crate::{
    call::pileup::Pileup,
    utils::{Base, Counter},
};

impl Pileup {
    ///  Calculate Shannon entropy for 100bp context around variant position
    pub(crate) fn entropy(&self) -> f64 {
        let idx = self.idx();
        let seq_context = self
            .segment
            .get(idx.saturating_sub(50), idx.saturating_add(51))
            .expect("sequence context indices are valid");

        let counts: Counter = seq_context.iter().map(|&b| Base::from(b)).collect();
        let total = seq_context.len() as f64;

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
}

#[cfg(test)]
mod tests {
    use crate::{
        call::pileup::SimpleReads,
        sequence::{ChunkRegion, Region, Segment},
        vcf::SequenceContext,
    };
    use proptest::proptest;
    use smallvec::smallvec;
    use smol_str::SmolStr;
    use std::{iter::repeat_n, sync::Arc};

    use super::*;

    #[test]
    fn low_entropy() {
        let region = ChunkRegion {
            region: Region { contig: SmolStr::new("chr13"), start: 1, end: 100 },
            last_position: 100,
        };
        let pileup = Pileup {
            region: region.clone(),
            context: SequenceContext::default(),
            segment: Arc::new(Segment { range: region, sequence: repeat_n(b'A', 100).collect() }),
            pos: 50,
            reads: SimpleReads(smallvec![]),
            reference_base: Base::A,
            is_cpg: false,
        };
        let entropy = pileup.entropy();
        assert_eq!(entropy, 0.0);
    }

    #[test]
    fn high_entropy() {
        let region = ChunkRegion {
            region: Region { contig: SmolStr::new("chr13"), start: 1, end: 100 },
            last_position: 100,
        };
        let pileup = Pileup {
            region: region.clone(),
            context: SequenceContext::default(),
            segment: Arc::new(Segment {
                range: region,
                sequence: repeat_n(b"ACTG", 25).flat_map(|x| *x).collect(),
            }),
            pos: 50,
            reads: SimpleReads(smallvec![]),
            reference_base: Base::A,
            is_cpg: false,
        };
        let entropy = pileup.entropy();
        assert!(entropy > 0.0);
    }

    proptest! {
        #[test]
        fn entropy_is_non_negative(sequence in "[ACTG]{100}") {
            let region = ChunkRegion {
                region: Region { contig: SmolStr::new("chr13"), start: 1, end: 100 },
                last_position: 100,
            };
            let pileup = Pileup {
                region: region.clone(),
                context: SequenceContext::default(),
                segment: Arc::new(Segment {
                    range: region,
                    sequence: sequence.into_bytes(),
                }),
                pos: 50,
                reads: SimpleReads(smallvec![]),
                reference_base: Base::A,
                is_cpg: false,
            };
            let entropy = pileup.entropy();
            assert!(entropy >= 0.0);
        }
    }
}
