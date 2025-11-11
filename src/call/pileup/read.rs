use crate::utils::{Base, Counter, Strand};
use smallvec::SmallVec;
use std::{fmt, ops::Deref};

/// A collection of bases seen in a pileup
#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[serde(transparent)]
pub struct SimpleReads(pub(crate) SmallVec<SimpleRead, 20>);

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for SimpleReads {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.0.iter()).finish()
    }
}

impl Deref for SimpleReads {
    type Target = [SimpleRead];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// A base seen in a pileup
#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[cfg_attr(test, derive(better_default::Default))] // for easier test construction
pub struct SimpleRead {
    /// The base seen
    #[cfg_attr(test, default(Base::Unknown))]
    pub base: Base,
    /// Base quality
    #[cfg_attr(test, default(30))]
    pub qual: u8,
    /// Mapping quality of the read this base belongs to
    #[cfg_attr(test, default(20))]
    pub mapq: u8,
    /// Strand the read was mapped to
    #[cfg_attr(test, default(Strand::OT))]
    pub strand: Strand,
    /// Whether the read was mapped to the reverse strand
    pub reverse: bool,
    /// Whether this base was seen in the first or second read of a pair
    pub second: bool,
    /// Position of the base in the read
    #[cfg_attr(test, default(PositionInRead { pos: 50, read_length: 100 }))]
    pub position: PositionInRead,
    /// Number of matching bases in the read this base belongs to
    #[cfg_attr(test, default(90))]
    pub matching_bases: u32,
    /// Number of indels in the read this base belongs to
    #[cfg_attr(test, default(2))]
    pub indels: u32,
    /// Query Name of the read this base belongs to
    // FIXME: Move out if only used for overlapping reads detection
    pub qname: SmallVec<u8, 42>,
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for SimpleRead {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.base)?;
        write!(f, " {} Q{} MQ{}", self.strand, self.qual, self.mapq,)
    }
}

impl SimpleReads {
    pub fn matches(&self, base: Base) -> bool {
        self.0.iter().all(|b| b.base == base)
    }

    pub fn is_variant_candidate(&self) -> bool {
        let counter: Counter = self.0.iter().map(|x| x.base).collect();
        counter.multiple_bases()
    }
}

#[derive(Clone, Copy, serde::Serialize, serde::Deserialize)]
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
