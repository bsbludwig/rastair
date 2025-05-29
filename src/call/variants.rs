use crate::utils::Base;
use smallvec::SmallVec;
use smol_str::SmolStr;
use std::{fmt, ops::Deref};

#[derive(Debug, Clone)]
pub struct VariantCandidatePileup {
    pub chrom: SmolStr,
    pub pos: u32,
    pub bases: SeenBases,
    pub reference_base: Base,
    pub next_base: Option<Base>,
}

impl VariantCandidatePileup {
    /// Is this a C->G variant candidate?
    pub fn is_cpg(&self) -> bool {
        self.reference_base == Base::C && self.next_base == Some(Base::G)
    }
}

/// A collection of bases seen in a pileup
#[derive(Clone)]
pub struct SeenBases(pub(crate) SmallVec<SeenBase, 20>);

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
    /// ```
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
        counter.interesting()
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

#[derive(Debug, Default)]
pub struct Counter {
    pub c: usize,
    pub t: usize,
    pub a: usize,
    pub g: usize,
}

impl Counter {
    /// Interesting if there are multiple different bases seen
    fn interesting(&self) -> bool {
        let mut count = 0;
        if self.c > 0 {
            count += 1;
        }
        if self.t > 0 {
            count += 1;
        }
        if self.a > 0 {
            count += 1;
        }
        if self.g > 0 {
            count += 1;
        }
        count >= 1
    }
}

impl FromIterator<Base> for Counter {
    fn from_iter<I: IntoIterator<Item = Base>>(iter: I) -> Self {
        let mut counter = Counter { c: 0, t: 0, a: 0, g: 0 };
        for c in iter {
            match c {
                Base::C => counter.c += 1,
                Base::T => counter.t += 1,
                Base::A => counter.a += 1,
                Base::G => counter.g += 1,
            }
        }
        counter
    }
}
