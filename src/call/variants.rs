use crate::{
    sequence::Segment,
    utils::{Base, Counter, Strand},
};
use better_default::Default;
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
#[cfg_attr(test, derive(Default))] // for easier test construction
pub struct SeenBase {
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
    pub qname: SmallVec<u8, 42>,
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

    fn default<T: Default>() -> T {
        T::default()
    }

    #[test]
    fn test_alleles_in_order() {
        let bases = SeenBases(SmallVec::from_vec(vec![
            SeenBase {
                base: Base::A,
                strand: Strand::OT,
                qname: SmallVec::from_vec(b"read1".to_vec()),
                ..default()
            },
            SeenBase {
                base: Base::C,
                strand: Strand::OB,
                reverse: true,
                second: false,
                qname: SmallVec::from_vec(b"read2".to_vec()),
                ..default()
            },
            SeenBase {
                base: Base::A,
                strand: Strand::OT,
                qname: SmallVec::from_vec(b"read3".to_vec()),
                ..default()
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
