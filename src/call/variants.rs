use crate::utils::{Base, TryAsBase as _};
use rust_htslib::bam::pileup::Alignment;
use smallvec::SmallVec;
use std::ops::Deref;

pub fn pileup_mapper(a: Alignment<'_>) -> Option<SeenBase> {
    let pos = a.qpos()?;
    let record = a.record();
    if !record.is_proper_pair() {
        // fixme: maybe be more lenient here
        return None;
    }
    if record.is_quality_check_failed() {
        return None;
    }
    // fixme: understand this better:
    // if record.cigar().iter().any(|c| matches!(c, Cigar::SoftClip(_))) {
    //     return None;
    // }

    Some(SeenBase {
        qname: SmallVec::from(record.qname()),
        base: record.seq()[pos].as_base().ok()?, // fixme: handle error or at least check usual error modes
        qual: record.qual()[pos],
        mapq: record.mapq(),
        reverse: record.is_reverse(),
        at_fringe: pos == 0 || pos == record.seq().len() - 1,
    })
}

#[derive(Debug)]
pub(crate) struct VariantCandidatePileup {
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
pub struct SeenBases(pub(crate) SmallVec<SeenBase, 20>);

#[cfg(not(tarpaulin_include))]
impl std::fmt::Debug for SeenBases {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
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
pub struct SeenBase {
    pub(crate) base: Base,
    pub(crate) qual: u8,
    pub(crate) mapq: u8,
    pub(crate) reverse: bool,
    pub(crate) at_fringe: bool,
    pub(crate) qname: SmallVec<u8, 16>,
}

#[cfg(not(tarpaulin_include))]
impl std::fmt::Debug for SeenBase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.base)?;
        write!(
            f,
            " {}{} {}/{}",
            if self.reverse { "rev" } else { "fwd" },
            if self.at_fringe { " fr" } else { "" },
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
        count > 1
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
