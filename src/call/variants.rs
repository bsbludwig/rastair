use super::scores::{self, Calc as _, StrandBias, VariantCandidatePileupMetrics};
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
        base: record.seq()[pos].as_base().unwrap(),
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
    pub fn metrics(&self) -> VariantCandidatePileupMetrics {
        let reference_bases = self.bases.iter().filter(|b| b.base == self.reference_base);
        let reference_count = reference_bases.clone().count();
        let alt_bases = self.bases.iter().filter(|b| b.base != self.reference_base);
        let alt_count = alt_bases.clone().count();

        let vaf = scores::VariantAlleleFrequency {
            reference_count: reference_count as u64,
            alt_count: alt_count as u64,
        };

        let binomial = scores::BinomialTest {
            reference_count: reference_count as u64,
            alt_count: alt_count as u64,
            error_rate: 0.01,
        };

        let mapq = scores::MappingQuality {
            reference_mapq: SmallVec::from_iter(reference_bases.clone().map(|b| b.mapq)),
            alt_mapq: SmallVec::from_iter(alt_bases.clone().map(|b| b.mapq)),
        };

        let baseq = scores::BaseQuality {
            reference_baseq: SmallVec::from_iter(reference_bases.clone().map(|b| b.qual)),
            alt_baseq: SmallVec::from_iter(alt_bases.clone().map(|b| b.qual)),
        };

        VariantCandidatePileupMetrics {
            reference_count,
            alt_count,
            vaf: vaf.calculate(),
            binomial: binomial.calculate(),
            mapq: mapq.calculate(),
            baseq: baseq.calculate(),
            strand_bias: StrandBias {
                reference_ot: reference_bases.clone().filter(|b| !b.reverse).count() as u64,
                reference_ob: reference_bases.clone().filter(|b| b.reverse).count() as u64,
                alt_ot: alt_bases.clone().filter(|b| !b.reverse).count() as u64,
                alt_ob: alt_bases.clone().filter(|b| b.reverse).count() as u64,
            }
            .calculate(),
        }
    }
}

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

pub struct SeenBase {
    base: Base,
    qual: u8,
    mapq: u8,
    reverse: bool,
    at_fringe: bool,
    qname: SmallVec<u8, 48>,
}

#[cfg(not(tarpaulin_include))]
impl std::fmt::Debug for SeenBase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if f.alternate() {
            // use the alternate format with `{:#?}` for no color
            write!(f, "{}", self.base)?;
        } else {
            write!(f, "{}", self.base.display_colored())?;
        }
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
struct Counter {
    c: usize,
    t: usize,
    a: usize,
    g: usize,
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
