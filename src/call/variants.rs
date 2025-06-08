use crate::{
    sequence::Segment,
    utils::{Base, Counter},
};
use color_eyre::eyre::ContextCompat as _;
use smallvec::SmallVec;
use smol_str::SmolStr;
use std::{fmt, ops::Deref, rc::Rc};

#[derive(Debug, Clone)]
pub struct VariantCandidatePileup {
    pub segment: Rc<Segment>,
    pub pos: u32,
    pub bases: SeenBases,
    pub reference_base: Base,
}

impl VariantCandidatePileup {
    /// Chromosome name of the segment
    pub fn chrom(&self) -> SmolStr {
        self.segment.range.chromosome.clone()
    }

    /// Position in the segment, 1-based
    fn idx(&self) -> usize {
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
    pub fn sequence_before(&self) -> SmallVec<Base, 2> {
        let idx = self.idx();

        self.segment.sequence_slice::<2>(idx.saturating_sub(2), idx).unwrap_or_default()
    }

    /// Sequence slice after the variant position
    pub fn sequence_after(&self) -> SmallVec<Base, 2> {
        let idx = self.idx();

        self.segment.sequence_slice::<2>(idx + 1, idx + 3).unwrap_or_default()
    }
}

/// A collection of bases seen in a pileup
#[derive(Clone)]
pub struct SeenBases(pub(crate) SmallVec<SeenBase, 20>);

impl SeenBases {
    pub fn alleles(&self) -> SmallVec<Base, 4> {
        self.iter().map(|b| b.base).fold(smallvec::SmallVec::new(), |mut acc, base| {
            if !acc.contains(&base) {
                acc.push(base);
            }
            acc
        })
    }

    pub fn alts(&self, reference: Base) -> SmallVec<Base, 4> {
        self.alleles().into_iter().filter(|base| reference != *base).collect::<SmallVec<Base, 4>>()
    }
}

#[test]
fn test_seen_bases_alts() {
    let bases = SeenBases(smallvec::smallvec![
        SeenBase {
            base: Base::A,
            qual: 30,
            mapq: 20,
            reverse: false,
            position: PositionInRead { pos: 0, read_length: 100 },
            matching_bases: 90,
            indels: 0,
            qname: SmallVec::from(&b"read1"[..]),
        },
        SeenBase {
            base: Base::C,
            qual: 30,
            mapq: 20,
            reverse: false,
            position: PositionInRead { pos: 1, read_length: 100 },
            matching_bases: 90,
            indels: 0,
            qname: SmallVec::from(&b"read1"[..]),
        },
    ]);
    assert_eq!(bases.alts(Base::A).as_slice(), &[Base::C]);
}

#[cfg(not(tarpaulin_include))]
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
    pub reverse: bool,
    pub position: PositionInRead,
    pub matching_bases: u32,
    pub indels: u32,
    /// Query Name of the read this base belongs to
    ///
    /// <!-- LLM explanation -->
    /// The `qname` field is the first mandatory field in each alignment record in a BAM/SAM file. It contains:
    ///
    /// - **Read identifier**: A unique string that identifies the sequencing read
    /// - **Paired-end reads**: For paired-end sequencing, both reads in a pair typically share the same `qname` base
    /// - **Format**: Usually follows the format from the sequencing instrument
    ///
    /// For example, in a SAM file (the text version of BAM), you might see:
    /// ```text
    /// SRR123456.1     99      chr1    1000    60      50M     =       1200    250     AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG    *
    /// SRR123456.2     147     chr1    1200    60      50M     =       1000    -250    TGCAGGCCTATGCAGCTGACTGCATAGCGTCAGCT    *
    /// ```
    ///
    /// In this example:
    /// - `SRR123456.1` and `SRR123456.2` are the `qname` values
    /// - These represent a paired-end read pair (same base name, different suffixes)
    ///
    /// The `qname` is essential for:
    /// - Tracking reads through analysis pipelines
    /// - Identifying paired reads
    /// - Debugging alignment issues
    /// - Linking back to original FASTQ files
    pub qname: SmallVec<u8, 16>,
}

#[cfg(not(tarpaulin_include))]
impl fmt::Debug for SeenBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.base)?;
        write!(
            f,
            " {}{} {}/{}",
            if self.reverse { "rev" } else { "fwd" },
            self.position,
            self.qual,
            self.mapq,
        )
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

impl fmt::Display for PositionInRead {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} / {}", self.pos, self.read_length)
    }
}
