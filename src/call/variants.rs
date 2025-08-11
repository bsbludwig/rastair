use crate::{
    sequence::Segment,
    utils::{Base, Counter, Strand},
};
use color_eyre::eyre::ContextCompat as _;
use smallvec::SmallVec;
use smol_str::SmolStr;
use std::{fmt, ops::Deref, rc::Rc};

#[derive(Debug, Clone)]
pub struct VariantCandidatePileup {
    pub segment: Rc<Segment>,
    /// Position in the sequence, 0-based
    pub pos: u32,
    pub bases: SeenBases,
    pub reference_base: Base,
    pub is_cpg: bool,
}

impl VariantCandidatePileup {
    /// Chromosome name of the segment
    pub fn chrom(&self) -> SmolStr {
        self.segment.range.contig.clone()
    }

    /// Position in the segment sequence, 0-based
    pub fn idx(&self) -> usize {
        let pos = usize::try_from(self.pos).expect("pos fits usize");
        usize::try_from(self.pos)
            .expect("position fits in usize")
            .checked_sub(usize::try_from(self.segment.range.start).expect("index fits in usize"))
            .wrap_err_with(|| {
                format!(
                    "pile position {} is not in segment range {}..{}",
                    pos, self.segment.range.start, self.segment.range.end
                )
            })
            .expect("valid index")
    }

    /// Sequence slice before the variant position
    pub fn sequence_before<const N: usize>(&self) -> SmallVec<Base, N> {
        let idx = self.idx();

        self.segment.sequence_slice::<N>(idx.saturating_sub(N), idx).unwrap_or_default()
    }

    /// Reference base right before the variant position
    pub fn ref_before(&self) -> Option<Base> {
        self.sequence_before::<1>().first().copied()
    }

    /// Sequence slice after the variant position
    pub fn sequence_after<const N: usize>(&self) -> SmallVec<Base, N> {
        let idx = self.idx();

        self.segment.sequence_slice::<N>(idx + 1, idx + N + 1).unwrap_or_default()
    }

    /// Reference base right after the variant position
    pub fn ref_after(&self) -> Option<Base> {
        self.sequence_after::<1>().first().copied()
    }

    pub fn alleles(&self) -> SmallVec<Base, 4> {
        let mut res = SmallVec::new();
        res.push(self.reference_base);
        self.bases.iter().map(|b| b.base).fold(res, |mut acc, base| {
            if !acc.contains(&base) {
                acc.push(base);
            }
            acc
        })
    }

    /// Get tuples of alleles (in order) and their corresponding evidence
    pub fn by_allele(&self) -> SmallVec<(Base, SmallVec<&SeenBase, 20>), 4> {
        self.alleles()
            .iter()
            .map(|base| {
                let matching_bases = self.bases.iter().filter(|b| b.base == *base).collect();
                (*base, matching_bases)
            })
            .collect()
    }

    pub fn alts(&self) -> SmallVec<Base, 4> {
        self.alleles()
            .into_iter()
            .filter(|base| self.reference_base != *base)
            .collect::<SmallVec<Base, 4>>()
    }
}

/// A collection of bases seen in a pileup
#[derive(Clone)]
pub struct SeenBases(pub(crate) SmallVec<SeenBase, 20>);

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for SeenBases {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.0.iter()).finish()
    }
}

impl Deref for SeenBases {
    type Target = [SeenBase];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// A base seen in a pileup
#[derive(Clone)]
pub struct SeenBase {
    pub base: Base,
    pub qual: u8,
    pub mapq: u8,
    pub strand: Strand,
    pub reverse: bool,
    pub position: PositionInRead,
    pub matching_bases: u32,
    pub indels: u32,
    // /// Query Name of the read this base belongs to
    // ///
    // /// <!-- LLM explanation -->
    // /// The `qname` field is the first mandatory field in each alignment record in a BAM/SAM file. It contains:
    // ///
    // /// - **Read identifier**: A unique string that identifies the sequencing read
    // /// - **Paired-end reads**: For paired-end sequencing, both reads in a pair typically share the same `qname` base
    // /// - **Format**: Usually follows the format from the sequencing instrument
    // ///
    // /// For example, in a SAM file (the text version of BAM), you might see:
    // /// ```text
    // /// SRR123456.1     99      chr1    1000    60      50M     =       1200    250     AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG    *
    // /// SRR123456.2     147     chr1    1200    60      50M     =       1000    -250    TGCAGGCCTATGCAGCTGACTGCATAGCGTCAGCT    *
    // /// ```
    // ///
    // /// In this example:
    // /// - `SRR123456.1` and `SRR123456.2` are the `qname` values
    // /// - These represent a paired-end read pair (same base name, different suffixes)
    // ///
    // /// The `qname` is essential for:
    // /// - Tracking reads through analysis pipelines
    // /// - Identifying paired reads
    // /// - Debugging alignment issues
    // /// - Linking back to original FASTQ files
    // pub qname: SmallVec<u8, 42>,
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for SeenBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.base)?;
        write!(f, " {} Q{} MQ{}", self.strand, self.qual, self.mapq,)
    }
}

impl SeenBases {
    pub fn matches(&self, base: Base) -> bool {
        self.0.iter().all(|b| b.base == base)
    }

    pub fn is_variant_candidate(&self) -> bool {
        let counter: Counter = self.0.iter().map(|x| x.base).collect();
        counter.multiple_bases()
    }
}

#[derive(Clone, Copy)]
pub struct PositionInRead {
    /// Position in the read, 0-based
    pub pos: u32,
    /// Length of the read
    pub read_length: u32,
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Display for PositionInRead {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} / {}", self.pos, self.read_length)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::{ChunkRegion, Region};
    use insta::assert_debug_snapshot;

    fn fake_segment() -> Rc<Segment> {
        Rc::new(Segment {
            range: ChunkRegion {
                region: Region { contig: "chr19".into(), start: 1000, end: 1100 },
                last_position: 2000,
            },
            sequence: vec![],
        })
    }

    #[test]
    fn test_alleles_in_order() {
        let bases = SeenBases(SmallVec::from_vec(vec![
            SeenBase {
                base: Base::A,
                qual: 30,
                mapq: 20,
                strand: Strand::OT,
                reverse: false,
                position: PositionInRead { pos: 0, read_length: 100 },
                matching_bases: 10,
                indels: 0,
                // qname: SmallVec::from_vec(b"read1".to_vec()),
            },
            SeenBase {
                base: Base::C,
                qual: 30,
                mapq: 20,
                strand: Strand::OB,
                reverse: true,
                position: PositionInRead { pos: 1, read_length: 100 },
                matching_bases: 5,
                indels: 0,
                // qname: SmallVec::from_vec(b"read2".to_vec()),
            },
            SeenBase {
                base: Base::A,
                qual: 30,
                mapq: 20,
                strand: Strand::OT,
                reverse: false,
                position: PositionInRead { pos: 2, read_length: 100 },
                matching_bases: 15,
                indels: 0,
                // qname: SmallVec::from_vec(b"read3".to_vec()),
            },
        ]));

        let variant_candidate = VariantCandidatePileup {
            segment: fake_segment(),
            pos: 1002, // Corresponds to index in the segment
            bases,
            reference_base: Base::T, // Assume T is the reference base at this position
            is_cpg: false,
        };

        let alleles = variant_candidate.alleles();
        assert_eq!(alleles[0], Base::T, "Reference base should be first");

        assert_debug_snapshot!(alleles, @r"
        [
            T,
            A,
            C,
        ]
        ");
    }
}
