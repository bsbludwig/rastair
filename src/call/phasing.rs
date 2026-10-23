//! Read-backed phasing.
//!
//! Two heterozygous sites carried by one DNA fragment are in known relative
//! phase, so the read pairs spanning them say which alleles share a
//! chromosome. The per-read observations have to be captured while a pileup
//! column is built, because `PileupMetrics` keeps no reads.
//!
//! Per-column data sits in `PileupMetrics.phase`; thresholds are constants
//! except the block cut, `--phase-min-switch-quality`. Phasing runs in every
//! output but CpG-only (`--cpgs-only`, or `--bed` alone), where the pre-filter
//! drops the variant sites before it could (`CallParams::calls_variants`). `problem` builds the sites
//! and fragments, `solver` splits them.
//!
//! **Observations** are captured inside the column loop. A column keeps its
//! whole list when it has `MIN_ALT_READS` reads of one non-reference base on
//! any strand (`candidate_observations`): a converted `T` over a `C` counts
//! too, since only the phasing gate knows whether the site is in a CpG
//! (de-novo CpGs are not known yet). An `N` reference or an empty column is
//! excluded, so `phase.is_some()` means there is something to use.
//! `accepts_read` (MAPQ and base quality at least 20) applies at collection,
//! so a read can count for depth and genotype without carrying phase.
//! `FragmentId` is a `NonZeroU64`: seqair reports "no qname" as a missing
//! hash, and bucketing those reads together would link every het in the
//! region. They are dropped, and `call` warns once per run, since a missing
//! qname is a whole-file property (CRAM with `RN=false`).
//!
//! **Blocks.** [`phase_segment`] runs before the overlap trim, because a site
//! in the overlap links reads reaching into the core, and takes the segment's
//! `is_core` predicate: each block is anchored on its first emitted site, so
//! `PS` names a record of this file and is unique across segments, and a block
//! with fewer than two emitted sites is not written. The solver returns raw
//! orientations, defined only up to a flip per block; `phase_segment` is the
//! one place that normalises them, onto the anchor, so a block opens on `0|1`.
//! `PS` and `PQ` are written only where the genotype is phased; homozygotes
//! never are.
//!
//! **Testing.** `tests/phasing_cli.rs` runs `call --phase` on a `Scenario`
//! (`tests/utils/scenario.rs`): a pseudo-random reference with no CpG unless
//! the test writes one, het SNVs and indels per haplotype, per-haplotype
//! methylation converted on the read's own strand, and read pairs per
//! haplotype and strand. `check_phasing` reads the VCF back and fails on a
//! switch error, a phased non-het, a `PS` not naming its block's first record,
//! or a block not opening on `0|1`. Read observations through this real path
//! rather than attaching them to a `pileups!` fixture, and make a scenario
//! decisive: the broken code should give a different block structure, not a
//! weaker edge (check by breaking the code it guards). A scenario linking two
//! blocks through one read pair gets two blocks: use two or more copies to
//! mean "linked"; a scenario whose reads all sit in one half of its contig
//! gives the other half a `DepthScale` below 1. `tests/data/test.bam` is a
//! poor phasing fixture (its few het calls are C>T on a fully methylated
//! control). `phasing_accuracy_on_an_na12878_slice` scores 200 kb of a chr12
//! NA12878 TAPS BAM against Platinum Genomes and is the only test that sees
//! phasing get worse on real reads; its reads stay out of git, so build the
//! fixture with `scripts/make_phasing_slice.py tmp/taps tmp/na12878_phasing`
//! and point `RASTAIR_PHASING_SLICE` at it, or it skips. Run it before merging
//! phasing changes.

use crate::{metrics::PileupMetrics, utils::logging::ThisIsABug as _};
use color_eyre::eyre::OptionExt as _;
use seqair_types::{Pos0, Pos1};
use std::num::{NonZeroU8, NonZeroU64};
use tracing::{debug, warn};

mod problem;
mod solver;

pub use problem::DepthScale;
use problem::{SiteAllele, build_problem};
pub use solver::DEFAULT_MIN_SWITCH_QUALITY;
use solver::{PhasedSite, solve};

/// Votes each allele of an SNV site needs.
const MIN_ALT_READS: u32 = 2;
/// Mapping quality a read needs to vote.
const MIN_MAPQ: u8 = 20;

/// What a haplotype is charged for disagreeing with one vote: the read's
/// preference between its two haplotypes, Phred-scaled and capped.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Cost(u8);

impl Cost {
    /// Most a single vote can be charged: one fragment's link is wrong
    /// 0.26 % of the time (Phred 26) whatever its base qualities say
    /// (chimeras, mismapped mates), so no single read may claim more.
    const CAP: u8 = 26;

    /// A read's log10 likelihood ratio between two haplotypes, as a Phred
    /// cost: ten times its size, rounded; an infinite ratio is the cap.
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "non-negative and below the cap; NaN becomes 0"
    )]
    fn of_log10_ratio(ratio: f32) -> Self {
        let phred = (10.0 * ratio.abs()).round();
        if phred >= f32::from(Self::CAP) { Self(Self::CAP) } else { Self(phred as u8) }
    }

    /// The measured Phred error of a raw switch quality, by bucket. Fragment
    /// errors share causes, so however many fragments add up, the error
    /// stops falling at about Q37.
    #[must_use]
    pub const fn calibrated(raw: i32) -> Self {
        Self(match raw {
            ..10 => 3,
            10..20 => 8,
            20..30 => 22,
            30..120 => 31,
            _ => 37,
        })
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

    /// The fragment of a read seqair reports a qname hash for.
    #[must_use]
    pub const fn of(qname_hash: NonZeroU64) -> Self {
        Self(qname_hash)
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

/// Which of a column's records a phase site is: the SNV record or the indel
/// record. One column can hold one of each, the SNV record written first.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiteKind {
    Snv,
    Indel,
}

/// Where a heterozygous site sits in its phase block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PhaseCall {
    /// VCF `PS`: the block's first emitted heterozygous site.
    pub set: PhaseSet,
    /// The allele haplotype 1 carries.
    pub first: GtAllele,
    /// VCF `PQ`: how sure the block is of this site's orientation.
    pub quality: Cost,
}

/// The phase of each of a column's records, once the segment is phased.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RecordPhases {
    pub snv: Option<PhaseCall>,
    pub indel: Option<PhaseCall>,
}

impl RecordPhases {
    const fn get(&self, kind: SiteKind) -> Option<PhaseCall> {
        match kind {
            SiteKind::Snv => self.snv,
            SiteKind::Indel => self.indel,
        }
    }

    const fn slot(&mut self, kind: SiteKind) -> &mut Option<PhaseCall> {
        match kind {
            SiteKind::Snv => &mut self.snv,
            SiteKind::Indel => &mut self.indel,
        }
    }
}

/// The phase call a column's `kind` record carries, if any.
pub(crate) fn call_of(pileup: &PileupMetrics, kind: SiteKind) -> Option<PhaseCall> {
    pileup.phase.as_deref()?.get(kind)
}

/// Phase a segment in place.
///
/// Runs before the overlap is trimmed, since a site in the overlap links reads
/// reaching into the core, but only the sites the segment emits (`is_core`)
/// get a call: `PS` names each block's first emitted site, so it names a
/// record of this file and is unique between segments, and a block left with
/// one emitted site loses its phase. Each block is oriented so that haplotype
/// 1 carries the lower allele at that site, which makes its first record
/// `0|1` rather than `1|0`.
pub fn phase_segment(
    pileups: &mut [PileupMetrics],
    is_core: impl Fn(u64) -> bool,
    depth_scale: Option<DepthScale>,
    min_switch_quality: u8,
) {
    let problem = build_problem(pileups, depth_scale);
    let solution = solve(&problem, min_switch_quality.into());

    let mut calls: Vec<(usize, SiteKind, PhaseCall)> = Vec::new();
    for block in &solution.blocks {
        let emitted: Vec<_> = block
            .iter()
            .filter(|phased| is_core(problem.sites[phased.site].pos.as_u64()))
            .collect();
        let (Some(&&PhasedSite { site: anchor, allele: anchor_allele, .. }), true) =
            (emitted.first(), emitted.len() >= 2)
        else {
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
        for &&PhasedSite { site: index, allele, quality } in &emitted {
            let site = &problem.sites[index];
            let allele = if flip { allele.flipped() } else { allele };
            let call = PhaseCall { set, first: site.allele(allele), quality };
            calls.push((site.pileup, site.kind, call));
        }
    }

    for pileup in pileups.iter_mut() {
        pileup.phase = None;
    }
    for (index, kind, call) in calls {
        // Sites index the slice they were built from.
        match pileups.get_mut(index).ok_or_eyre("Phase site outside its segment").this_is_a_bug() {
            Ok(pileup) => *pileup.phase.get_or_insert_default().slot(kind) = Some(call),
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
