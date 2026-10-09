//! Read-backed phasing.
//!
//! Two heterozygous sites carried by one DNA fragment are in known relative
//! phase, so the read pairs spanning them say which alleles share a
//! chromosome. The per-read observations have to be captured while a pileup
//! column is built, because `PileupMetrics` keeps no reads; only the seqair
//! backend does that.

use crate::{
    metrics::{PileupMetrics, methylation::CpgSide},
    utils::logging::ThisIsABug as _,
};
use color_eyre::eyre::OptionExt as _;
use seqair_types::{Base, BaseQuality, Pos0, Pos1, Strand};
use std::num::{NonZeroU8, NonZeroU64};
use tracing::{debug, warn};

mod problem;
mod solver;

use problem::{SiteAllele, build_problem};
use solver::solve;

/// Observations on one alternative base a column needs before its reads are
/// kept for phasing.
const MIN_ALT_READS: u32 = 2;
/// Mapping quality a read needs to contribute a phase observation.
const MIN_MAPQ: u8 = 20;
/// Base quality a read needs to contribute a phase observation.
const MIN_BASEQ: u8 = 20;

/// Is this read good enough to contribute a phase observation?
///
/// Stricter than the calling filters, so a read can count towards depth and
/// genotype without carrying phase.
#[must_use]
pub fn accepts_read(mapq: u8, baseq: BaseQuality) -> bool {
    mapq >= MIN_MAPQ && baseq.get().is_some_and(|baseq| baseq >= MIN_BASEQ)
}

/// What a haplotype is charged for disagreeing with one observation: the
/// observation's base quality, capped.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Cost(u8);

impl Cost {
    /// Most a single observation can be charged.
    const CAP: u8 = 40;

    const fn of(qual: BaseQuality) -> Self {
        match qual.get() {
            Some(qual) if qual < Self::CAP => Self(qual),
            Some(_) => Self(Self::CAP),
            None => Self(0),
        }
    }
}

impl From<Cost> for i32 {
    fn from(cost: Cost) -> Self {
        Self::from(cost.0)
    }
}

/// The template a read belongs to: seqair's seed-fixed hash of its qname.
///
/// Both mates share it, and it is stable across segments and runs. Two
/// templates whose qnames collide merge into one fragment, so everything built
/// on it tolerates a merge rather than assuming a bijection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct FragmentId(NonZeroU64);

impl FragmentId {
    /// The fragment a read with this qname hash belongs to.
    ///
    /// A read without a qname has no hash at all (`qname_hash()` is `None`) and
    /// never gets here; lumping such reads together would link every
    /// heterozygote in the region, so they are dropped instead. A zero hash is
    /// refused for the same reason.
    #[must_use]
    pub fn new(qname_hash: u64) -> Option<Self> {
        NonZeroU64::new(qname_hash).map(Self)
    }
}

/// The 1-based position a set of phased records is filed under, written out as
/// VCF `PS`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PhaseSet(Pos1);

impl PhaseSet {
    /// The set a column at this 0-based position anchors.
    #[must_use]
    pub fn at(pos: Pos0) -> Option<Self> {
        pos.to_one_based().ok().map(Self)
    }

    #[must_use]
    pub const fn as_i32(self) -> i32 {
        self.0.as_i32()
    }
}

/// An allele numbered the way [`GenotypeTag`](crate::call::variant_calling::GenotypeTag)
/// numbers it: the reference, or the `n`-th entry of `PileupMetrics::alts`.
///
/// The VCF writer remaps this onto `ALT` column order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum GtAllele {
    Ref,
    Alt(NonZeroU8),
}

/// One read's base at a column that might take part in phasing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PhaseObservation {
    pub fragment: FragmentId,
    pub base: Base,
    /// Needed to tell a converted base from a variant.
    pub strand: Strand,
    /// Already past [`accepts_read`].
    pub qual: BaseQuality,
}

/// Where a heterozygous column sits in its phase block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PhaseCall {
    /// VCF `PS`: the block's first emitted heterozygous site.
    pub set: PhaseSet,
    /// The allele haplotype 1 carries.
    pub first: GtAllele,
}

/// What phasing keeps on one column: its reads' observations until the segment
/// is phased, then its place in a block, if it has one.
#[derive(Debug, Clone)]
pub enum ColumnPhase {
    Observed(Box<[PhaseObservation]>),
    Phased(PhaseCall),
}

/// The phase call a column carries, if any.
pub(crate) fn call_of(pileup: &PileupMetrics) -> Option<PhaseCall> {
    match pileup.phase.as_deref() {
        Some(ColumnPhase::Phased(call)) => Some(*call),
        _ => None,
    }
}

fn observations_of(pileup: &PileupMetrics) -> Option<&[PhaseObservation]> {
    match pileup.phase.as_deref() {
        Some(ColumnPhase::Observed(observations)) => Some(observations),
        _ => None,
    }
}

/// Could TAPS conversion alone explain this mismatch against the reference?
///
/// Deliberately blind to sequence context: it decides which columns are
/// collected, before anything knows where the de-novo CpGs are.
fn taps_explains_mismatch(reference_base: Base, base: Base, strand: Strand) -> bool {
    [CpgSide::C, CpgSide::G].into_iter().any(|side| {
        reference_base == side.unmod_base() && base == side.mod_base() && strand == side.strand()
    })
}

/// Keep a column's observations only where it may be heterozygous: at least
/// [`MIN_ALT_READS`] mismatches on one base that TAPS conversion cannot
/// explain. The whole column is kept, TAPS-ambiguous reads included.
///
/// `Some` always carries at least one observation.
#[must_use]
pub fn candidate_observations(
    observations: &[PhaseObservation],
    reference_base: Base,
) -> Option<Box<[PhaseObservation]>> {
    // An `N` in the reference can never be genotyped.
    reference_base.known_index()?;

    let mut informative = [0u32; Base::KNOWN.len()];
    for obs in observations.iter().filter(|obs| {
        obs.base != reference_base && !taps_explains_mismatch(reference_base, obs.base, obs.strand)
    }) {
        if let Some(count) = obs.base.known_index().and_then(|idx| informative.get_mut(idx)) {
            *count += 1;
        }
    }

    informative.iter().any(|&count| count >= MIN_ALT_READS).then(|| observations.into())
}

/// Phase a segment in place, and drop the observations, which nothing reads
/// after this.
///
/// Runs before the overlap is trimmed, since a site in the overlap links reads
/// reaching into the core, but only the sites the segment emits (`is_core`)
/// get a call: `PS` names each block's first emitted site, so it names a
/// record of this file and is unique between segments, and a block left with
/// one emitted site loses its phase. Each block is oriented so that haplotype
/// 1 carries the lower allele at that site, which makes its first record
/// `0|1` rather than `1|0`.
pub fn phase_segment(pileups: &mut [PileupMetrics], is_core: impl Fn(u64) -> bool) {
    let problem = build_problem(pileups);
    let solution = solve(&problem);

    let mut calls: Vec<(usize, PhaseCall)> = Vec::new();
    for block in &solution.blocks {
        let emitted: Vec<_> = block
            .iter()
            .filter(|&&(index, _)| is_core(problem.sites[index].pos.as_u64()))
            .collect();
        let (Some(&&(anchor, anchor_allele)), true) = (emitted.first(), emitted.len() >= 2) else {
            continue;
        };
        let pos = problem.sites[anchor].pos;
        // Fails only at i32::MAX, which no BAM contig reaches.
        let set = match PhaseSet::at(pos).ok_or_eyre("No 1-based position").this_is_a_bug() {
            Ok(set) => set,
            Err(error) => {
                warn!(?error, %pos, "Skipping a phase block whose anchor has no PS");
                continue;
            }
        };
        let flip = anchor_allele == SiteAllele::Second;
        for &&(index, allele) in &emitted {
            let site = &problem.sites[index];
            let allele = if flip { allele.flipped() } else { allele };
            calls.push((site.pileup, PhaseCall { set, first: site.allele(allele).allele }));
        }
    }

    for pileup in pileups.iter_mut() {
        pileup.phase = None;
    }
    for (index, call) in calls {
        // Sites index the slice they were built from.
        match pileups.get_mut(index).ok_or_eyre("Phase site outside its segment").this_is_a_bug() {
            Ok(pileup) => pileup.phase = Some(Box::new(ColumnPhase::Phased(call))),
            Err(error) => warn!(?error, index, "Dropping a phase call"),
        }
    }

    let stats = solution.stats;
    debug!(
        sites = stats.sites,
        phased_sites = stats.phased_sites,
        blocks = stats.blocks,
        largest_block = stats.largest_block,
        mec_cost = stats.mec_cost,
        "Phased segment"
    );
}

#[cfg(test)]
pub(crate) mod fixtures;

#[cfg(test)]
mod tests;
