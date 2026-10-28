//! Pair-HMM scores of a cluster of variant columns' reads against its local
//! haplotypes.
//!
//! Once a sub-segment's columns are drafted, the reads spanning a cluster of
//! indel and SNV columns are aligned by a TAPS-aware pair-HMM against the
//! reference window and against the window with every compatible subset of
//! the cluster's [`Event`]s applied. `PileupMetrics` keeps no reads, so what
//! is stored is one row of log10 likelihoods per read, with the [`Verdict`]
//! they add up to: the diploid pair of haplotypes and each event's [`Margin`]
//! of presence. The scoring is compair's; the stored scores are plain data.
//!
//! **Haplotypes are subsets of events.** [`HapIdx`] is a bitmask over
//! [`EventIdx`]; [`Cluster::new`] keeps the subsets one chromosome can carry
//! (`Event::compatible`: not two substitutions or two indels at one anchor,
//! no deletion over another event's anchor), so at most eight haplotypes and
//! alignments per read. A substitution and an indel at one anchor share one:
//! the indel keeps the anchor base and the substitution changes it. The
//! scoring applies a subset's events later anchor first, so an earlier
//! anchor's window index still holds, and keeps one haplotype per sequence
//! (`drop_respellings`, the fewest events): two subsets can spell one sequence
//! (an inserted `T` before a run plus the deletion of the `C` after it is the
//! `C>T`; an insertion and a deletion in one run are the reference), and
//! identical haplotypes tie every pair over them. A dropped event is never
//! called. Subsets matter: two homozygous deletions 13 bp apart need the
//! haplotype carrying both, or one loses its reads to the other's mismatches.
//!
//! **The verdict** ([`Verdict`], once per cluster) is over pairs of
//! haplotypes. An `Account` is a pair, or a pair with the reads of one
//! `Artefact` mixed in at `ARTEFACT_FRACTION` less
//! `ARTEFACT_PRIOR_LOG10`: an event the pair lacks, or the reference for a
//! pair without it. That keeps 7 alt reads in 73 from being a genotype while 4
//! in 8 stay a het, and one stray reference read from making a homozygote a
//! het. Margins are per event, the best account whose pair carries it over the
//! best whose pair does not; an event is called at [`Margin::PRESENT`] and
//! QUAL is ten times the margin. [`Verdict::called_pair`] is the likeliest
//! pair carrying exactly the called events; it decides which indel events are
//! called, not the `GT` a record writes (see **Genotypes**).
//! `genotype_quality` is the called pair's flat-prior quality over the
//! pure pairs, low when it is not the best pair; the phasing gate reads it.
//! No record writes it: a record's `GQ` is over its own alleles.
//!
//! **A read's likelihood under a pair** ([`ReadScores`]):
//!
//! - Each read is mixed with its mismapping chance `10^(-MAPQ/10)`, explained
//!   by its best haplotype (`Mismap`), so a stray low-MAPQ read costs any pair
//!   at most `MAPQ / 10` log10. An unavailable MAPQ (255) counts as 60
//!   (`MappingQualityExt::or_unique`).
//! - Every CpG is scored at `UNKNOWN_BETA` (½), which cancels wherever both
//!   haplotypes of a pair carry it. In an all-substitution cluster a CpG one
//!   haplotype has and another lacks (a de-novo `T>C`, a `C>T` at a reference
//!   CpG) is profiled instead: its carrier is scored at the converting strand's
//!   level as the pair's only carrier ([`ReadScores::log10_alone`]) and as both
//!   chromosomes (`log10_doubled`), and a pair reads the doubled score when the
//!   partner carries all the haplotype's divergent CpGs. At ½ every converted
//!   OT read of a methylated hom-alt de-novo CpG costs the homozygote.
//!   Shortcuts: clusters with an indel stay at ½; `log10_alone` is what
//!   `reads_over` and phasing votes read, having no pair; a haplotype whose
//!   partner carries only some of its CpGs is scored alone at all of them;
//!   `Partner::at` judges the partner from the non-carriers' consensus, else
//!   the reference.
//!
//! **`AD` and `DP`** are per record ([`HaplotypeScores::reads_over`]): a read
//! counts for the allele whose likeliest haplotype beats every other allele's
//! by `INFORMATIVE_LOG10` (GATK's informative-read threshold), the cluster's
//! other events free, so a read over one event of a cluster says as much
//! about it as a read over all. It is far below `Margin::PRESENT` on
//! purpose: the pair-HMM sums a tract gap over every placement, so one read of
//! a 1 bp deletion in a 13-base homopolymer is only ~0.9 log10 likelier under
//! the deletion. `AD` sums to `DP` but for the reads of another allele at the
//! record's site.
//!
//! **Genotypes.** A record's `GT`, `GQ` and `PL` all come from
//! [`Verdict::genotype_likelihoods`], its genotypes over its own alleles with
//! the cluster's other events held at their calls, each the sum of its pairs
//! at their best accounts. Which alleles are present is decided elsewhere (an
//! indel's by the called pair, an SNV's by the forests), so `GT` is the
//! likeliest genotype carrying an alt
//! ([`GenotypeLikelihoods::likeliest_with_an_alt`]), for indel and SNV
//! records alike, and `GQ` and `PL` are over every genotype, `0/0` included.
//! The called pair's own genotype can differ where cis and trans pairs of a
//! het tie and sum past the homozygote.
//! No genotype prior: a het prior of 0.184 log10 (GIAB's het:hom ratio) is
//! finer than the likelihood's tie margins are calibrated, so it does not
//! compose with them; it cost both SNV and indel records on chr12. It is not a
//! per-site marginal either: that one's "without" is a free best-of-N over the
//! cluster's other alleles, so at a hom indel beside a stutter or
//! twin-spelling candidate every read is neutral and one copy wins. Where no
//! read without a linked het deletion reaches an SNV, `(ref, del+SNV)` and
//! `(SNV, del+SNV)` tie and the het or the homozygote is written at a low GQ:
//! the model's honest tie, as is such a call on MAPQ 0-19 reads.
//!
//! **TAPS twin spellings.** A deletion next to a methylated CpG has two
//! spellings whose haplotypes differ by a `C`/`T` at the CpG; converted OT
//! reads match the wrong one for free, so at ½ it sometimes collects them. A
//! cluster calls one of the two; which one is still open.

use crate::call::{
    phasing::FragmentId,
    pileup::indels::IndelAllele,
    variant_calling::{GenotypeTag, GtAllele},
};
use crate::utils::MappingQualityExt;
use seqair_types::{Base, MappingQuality, Phred, Pos0, SmallVec};
use std::num::NonZeroU8;

mod scoring;
#[cfg(test)]
pub use scoring::score;
pub use scoring::{HaplotypeRead, HaplotypeScratch, Prepared, prepare};

/// Events, indels and substitutions, a cluster scores at most: every subset of them is a
/// haplotype, so this bounds the alignments per read at `2^MAX_EVENTS`.
pub const MAX_EVENTS: usize = 3;
/// One haplotype per subset of the events.
pub const MAX_HAPLOTYPES: usize = 1 << MAX_EVENTS;
/// Events one anchor contributes at most: two that exclude each other (two
/// substitutions, or two indels) fill a pair there, and a third of their kind
/// could never be called. The limit counts every event, so a substitution
/// beside two indels at one anchor is left out too.
const MAX_PER_ANCHOR: usize = 2;

/// A three-wide score row, reference and two alternatives at one anchor,
/// widened to [`MAX_HAPLOTYPES`]: with two indels (or two substitutions) at
/// one anchor the haplotypes are `{}`, `{0}` and `{1}`, at indices 0, 1 and 2.
#[cfg(test)]
pub(crate) const fn row(three: [f32; 3]) -> [f32; MAX_HAPLOTYPES] {
    let mut log10 = [f32::NEG_INFINITY; MAX_HAPLOTYPES];
    log10[0] = three[0];
    log10[1] = three[1];
    log10[2] = three[2];
    log10
}

/// What an event changes: the base at its anchor, or an indel after it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Allele {
    Base(Base),
    Indel(IndelAllele),
}

impl Allele {
    #[must_use]
    pub fn as_indel(&self) -> Option<&IndelAllele> {
        match self {
            Self::Indel(indel) => Some(indel),
            Self::Base(_) => None,
        }
    }

    /// Diagonals a read carrying this allele sits off the reference.
    #[must_use]
    fn shift(&self) -> usize {
        self.as_indel().map_or(0, IndelAllele::len)
    }
}

impl From<IndelAllele> for Allele {
    fn from(indel: IndelAllele) -> Self {
        Self::Indel(indel)
    }
}

/// An allele at the column it is anchored on: a substitution's own column,
/// or the base before an indel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Event {
    pub anchor: Pos0,
    pub allele: Allele,
}

impl Event {
    /// Whether one chromosome can carry both: not two substitutions or two
    /// indels at one anchor, and neither a deletion over the other's anchor,
    /// which for a substitution is the base it changes. A substitution and an
    /// indel at one anchor can: the indel keeps the anchor base, and the
    /// substitution changes it.
    fn compatible(&self, other: &Self) -> bool {
        let deletes_over = |a: &Self, b: &Self| match &a.allele {
            Allele::Indel(IndelAllele::Deletion(deleted)) => {
                let end = a.anchor.as_u64().saturating_add(deleted.len() as u64);
                b.anchor > a.anchor && b.anchor.as_u64() <= end
            }
            Allele::Indel(IndelAllele::Insertion(_)) | Allele::Base(_) => false,
        };
        let substitution = |event: &Self| matches!(event.allele, Allele::Base(_));
        let excludes = self.anchor == other.anchor && substitution(self) == substitution(other);
        !excludes && !deletes_over(self, other) && !deletes_over(other, self)
    }
}

/// The events scored together, and the haplotypes they make.
///
/// At most [`MAX_EVENTS`] events and `MAX_PER_ANCHOR` per anchor. [`Verdict`]
/// always finds a pair carrying exactly the events it calls present: they sit
/// in the best pair, and any subset of a haplotype's events is a haplotype.
#[derive(Debug, Clone)]
pub struct Cluster {
    /// An indel cluster's in record order: most reads first, then anchor,
    /// then the allele, as `cluster_candidates` ranks them. An SNV cluster
    /// has its indel chain's in that order, then each SNV column's
    /// substitutions by column, the first-stage forest's favourite first. A column's record covers
    /// the events anchored at it.
    events: SmallVec<Event, MAX_EVENTS>,
    /// The subsets of `events` one chromosome can carry, ascending.
    haplotypes: SmallVec<HapIdx, MAX_HAPLOTYPES>,
    /// By [`HapIdx::index`], a bitmask over the divergent `CpG` positions
    /// (a role at a window index some haplotypes have and others not) the
    /// haplotype carries; all zero unless the scoring profiled their levels.
    divergent_cpgs: [u32; MAX_HAPLOTYPES],
}

impl Cluster {
    /// `None` when `events` breaks the limits [`Self::admits`] checks.
    #[must_use]
    pub fn new(events: SmallVec<Event, MAX_EVENTS>) -> Option<Self> {
        let admitted = events
            .iter()
            .enumerate()
            .all(|(i, event)| events.get(..i).is_some_and(|before| Self::admits(before, event)));
        admitted.then(|| Self {
            haplotypes: HapIdx::all_valid(&events),
            events,
            divergent_cpgs: [0; MAX_HAPLOTYPES],
        })
    }

    /// Whether `partner` carries every divergent `CpG` `hap` does, so a pair
    /// of the two scores `hap` with the levels profiled over both chromosomes.
    fn shares_cpgs(&self, hap: HapIdx, partner: HapIdx) -> bool {
        let mask = |h: HapIdx| self.divergent_cpgs.get(h.index()).copied().unwrap_or(0);
        mask(hap) & !mask(partner) == 0
    }

    /// Whether `candidate` can join a cluster of `events`.
    #[must_use]
    pub fn admits(events: &[Event], candidate: &Event) -> bool {
        events.len() < MAX_EVENTS
            && events.iter().filter(|event| event.anchor == candidate.anchor).count()
                < MAX_PER_ANCHOR
    }

    #[must_use]
    pub fn events(&self) -> &[Event] {
        &self.events
    }

    #[must_use]
    pub fn haplotypes(&self) -> &[HapIdx] {
        &self.haplotypes
    }
}

/// Position of an event in [`Cluster::events`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EventIdx(u8);

impl EventIdx {
    #[must_use]
    pub fn new(position: usize) -> Option<Self> {
        (position < MAX_EVENTS).then(|| u8::try_from(position).ok()).flatten().map(Self)
    }

    #[must_use]
    pub fn position(self) -> usize {
        usize::from(self.0)
    }
}

/// One of a cluster's haplotypes: the reference window with a subset of the
/// events applied, the subset as a bitmask over [`EventIdx`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HapIdx(u8);

impl HapIdx {
    pub const REFERENCE: Self = Self(0);

    /// The haplotype carrying `event` alone.
    #[must_use]
    pub fn of_event(event: EventIdx) -> Self {
        Self(1 << event.0)
    }

    /// Position in [`ReadScores::log10_alone`].
    #[must_use]
    pub fn index(self) -> usize {
        usize::from(self.0)
    }

    #[must_use]
    pub fn carries(self, event: EventIdx) -> bool {
        self.0 & (1 << event.0) != 0
    }

    /// The events this haplotype carries, ascending.
    pub fn events(self) -> impl Iterator<Item = EventIdx> {
        (0..MAX_EVENTS).filter_map(EventIdx::new).filter(move |event| self.carries(*event))
    }

    /// Every subset of `events` one chromosome can carry, ascending; the
    /// reference first and every singleton among them.
    fn all_valid(events: &[Event]) -> SmallVec<Self, MAX_HAPLOTYPES> {
        let n = events.len().min(MAX_EVENTS);
        (0..1u8 << n)
            .map(Self)
            .filter(|hap| {
                let carried: SmallVec<&Event, MAX_EVENTS> =
                    hap.events().filter_map(|e| events.get(e.position())).collect();
                carried
                    .iter()
                    .enumerate()
                    .all(|(i, a)| carried.iter().skip(i + 1).all(|b| a.compatible(b)))
            })
            .collect()
    }
}

/// An unordered pair of a cluster's haplotypes: a diploid genotype over them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HapPair {
    low: HapIdx,
    high: HapIdx,
}

impl HapPair {
    #[must_use]
    pub fn new(a: HapIdx, b: HapIdx) -> Self {
        Self { low: a.min(b), high: a.max(b) }
    }

    /// Every pair over `haplotypes`, in VCF `GL` order over that list:
    /// `(0,0), (0,1), (1,1), (0,2), (1,2), (2,2), ...`.
    fn all(haplotypes: &[HapIdx]) -> impl Iterator<Item = Self> + '_ {
        haplotypes.iter().enumerate().flat_map(|(j, &high)| {
            haplotypes.iter().take(j + 1).map(move |&low| Self::new(low, high))
        })
    }

    /// The two haplotypes, lower first.
    #[must_use]
    pub fn haplotypes(self) -> (HapIdx, HapIdx) {
        (self.low, self.high)
    }

    /// How many of the pair's chromosomes carry `event`.
    #[must_use]
    pub fn copies_of(self, event: EventIdx) -> u8 {
        u8::from(self.low.carries(event)) + u8::from(self.high.carries(event))
    }

    #[must_use]
    pub fn carries(self, event: EventIdx) -> bool {
        self.copies_of(event) > 0
    }

    #[must_use]
    pub fn is_heterozygous(self) -> bool {
        self.low != self.high
    }

    /// For an event on exactly one chromosome: the haplotype carrying it and
    /// the one not.
    #[must_use]
    pub fn split(self, event: EventIdx) -> Option<(HapIdx, HapIdx)> {
        match (self.low.carries(event), self.high.carries(event)) {
            (true, false) => Some((self.low, self.high)),
            (false, true) => Some((self.high, self.low)),
            _ => None,
        }
    }
}

/// One read's likelihood under each of a cluster's haplotypes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReadScores {
    /// `None` for a read without a qname.
    pub fragment: Option<FragmentId>,
    pub mapq: MappingQuality,
    /// `log10 P(read | haplotype)` by [`HapIdx::index`]; `-inf` for a subset
    /// that is not a haplotype, or that the read's band misses. At a
    /// haplotype's divergent `CpG`s, the level profiled for it as the pair's
    /// only carrier.
    log10_alone: [f32; MAX_HAPLOTYPES],
    /// `log10_alone` with the divergent `CpG`s at the level profiled for both
    /// chromosomes carrying them; equal to it for a haplotype without any.
    log10_doubled: [f32; MAX_HAPLOTYPES],
}

impl ReadScores {
    /// A row whose haplotypes carry no divergent `CpG`.
    #[must_use]
    pub fn new(
        fragment: Option<FragmentId>,
        mapq: MappingQuality,
        log10: [f32; MAX_HAPLOTYPES],
    ) -> Self {
        Self { fragment, mapq, log10_alone: log10, log10_doubled: log10 }
    }

    /// The likelihood under `hap`, its divergent `CpG`s at the level profiled
    /// as the pair's only carrier. shortcut: whatever the partner, so
    /// [`HaplotypeScores::reads_over`] (AD/DP) and phasing votes, which have
    /// no pair, read this.
    #[must_use]
    pub fn log10_alone(&self, hap: HapIdx) -> f32 {
        self.log10_alone.get(hap.index()).copied().unwrap_or(f32::NEG_INFINITY)
    }

    /// The likelihood under `hap`, its divergent `CpG`s at the level profiled
    /// for both chromosomes.
    fn log10_doubled(&self, hap: HapIdx) -> f32 {
        self.log10_doubled.get(hap.index()).copied().unwrap_or(f32::NEG_INFINITY)
    }

    /// The likelihood under `hap` as one chromosome of a pair with `partner`.
    fn log10_beside(&self, cluster: &Cluster, hap: HapIdx, partner: HapIdx) -> f32 {
        if cluster.shares_cpgs(hap, partner) {
            self.log10_doubled(hap)
        } else {
            self.log10_alone(hap)
        }
    }

    /// `log10 P(read | a) - log10 P(read | b)`, each at its alone level
    /// ([`Self::log10_alone`]): positive when the read prefers
    /// `a`, infinite when one of the two cannot have produced it. `None` when
    /// neither can, since the read then says nothing between them.
    #[cfg(test)]
    #[must_use]
    pub fn log10_ratio(&self, a: HapIdx, b: HapIdx) -> Option<f32> {
        let (a, b) = (self.log10_alone(a), self.log10_alone(b));
        (a.is_finite() || b.is_finite()).then_some(a - b)
    }

    /// Whether the read is at least [`Margin::PRESENT`] likelier under `hap`
    /// than under `over`: the same "ten times likelier" an event is called on.
    #[cfg(test)]
    #[must_use]
    pub fn prefers(&self, hap: HapIdx, over: HapIdx) -> bool {
        self.log10_ratio(hap, over).is_some_and(|ratio| f64::from(ratio) >= Margin::PRESENT.0)
    }

    /// `log10 P(read | pair)`, the read drawn from either haplotype of the pair
    /// with equal chance.
    #[cfg(test)]
    fn log10_under(&self, cluster: &Cluster, pair: HapPair) -> f64 {
        self.mismap(cluster).mix(self.log10_unmixed(cluster, pair))
    }

    fn log10_unmixed(&self, cluster: &Cluster, pair: HapPair) -> f64 {
        let (a, b) = pair.haplotypes();
        log10_mean(
            f64::from(self.log10_beside(cluster, a, b)),
            f64::from(self.log10_beside(cluster, b, a)),
        )
    }

    /// The read's mismapping chance `10^(-MAPQ/10)` over `cluster`.
    fn mismap(&self, cluster: &Cluster) -> Mismap {
        let mismapped = 10f64.powf(-f64::from(self.mapq.or_unique()) / 10.0);
        let best = cluster
            .haplotypes()
            .iter()
            .map(|&hap| f64::from(self.log10_alone(hap).max(self.log10_doubled(hap))))
            .fold(f64::NEG_INFINITY, f64::max);
        Mismap { mapped: (1.0 - mismapped).log10(), stray: mismapped.log10() + best }
    }

    /// `log10 P(read | pair contaminated by hap)`, before the [`Mismap`]: the
    /// read drawn from `hap` with chance `fraction`, else from the pair.
    fn log10_contaminated(
        &self,
        cluster: &Cluster,
        pair: HapPair,
        hap: HapIdx,
        fraction: f64,
    ) -> f64 {
        log10_sum(
            self.log10_unmixed(cluster, pair) + (1.0 - fraction).log10(),
            f64::from(self.log10_alone(hap)) + fraction.log10(),
        )
    }
}

/// A read's `log10` under a genotype is mixed with its mismapping chance: a
/// mismapped read is explained by whichever haplotype fits it best, so it
/// costs no genotype more than `MAPQ / 10` log10. Both terms depend on the
/// read alone, so they are worked out once per read and cluster.
#[derive(Debug, Clone, Copy)]
struct Mismap {
    /// `log10` of the chance the read is mapped right.
    mapped: f64,
    /// `log10` of the chance it is mismapped, plus its best haplotype's.
    stray: f64,
}

impl Mismap {
    fn mix(self, log10: f64) -> f64 {
        log10_sum(self.mapped + log10, self.stray)
    }
}

/// The fraction at which an event is tried as an artefact: reads drawn
/// from its haplotype this often on top of a genotype without it, as
/// stutter, a mis-mapped paralog or a sample contaminant would supply them. A
/// true event has to beat the likeliest of these by [`Margin::PRESENT`] on
/// top of [`ARTEFACT_PRIOR_LOG10`], which is what keeps seven alt reads in
/// seventy-three from being a genotype while four in eight still are. On
/// chr12 (26x and 73x) adding 0.02 or 0.05 beside it changed no call, and 0.2
/// cost insertions at 26x.
const ARTEFACT_FRACTION: f64 = 0.1;

/// How much likelier a read must be under its best allele than under the
/// runner-up to count for it in `AD` and `DP` (GATK's informative-read
/// threshold). Far below [`Margin::PRESENT`] on purpose: the pair-HMM sums a
/// tract indel over every placement of its gap, so one read of a 1 bp
/// deletion in a 13-base homopolymer is only ~0.9 log10 likelier under the
/// deletion than under the reference, and a read-level "ten times" would
/// count no read at all at exactly the indels that are called most.
const INFORMATIVE_LOG10: f32 = 0.2;

/// Prior odds against an artefact explanation, as `log10`: an artefact
/// account is taken at face value only when it beats the genotype a hundred
/// times over. Chosen on chr12: at 1 it cost 0.006 insertion F1 at 26x, at 0
/// ten times that; at 2 it is neutral at 26x and +0.005 / +0.011 at 73x.
const ARTEFACT_PRIOR_LOG10: f64 = 2.0;

/// A cluster's reads scored against its haplotypes, and what they add up to.
#[derive(Debug, Clone)]
pub struct HaplotypeScores {
    cluster: Cluster,
    reads: Box<[ReadScores]>,
    /// Derived from `reads`.
    verdict: Verdict,
}

/// A record's reads: its `AD` and `DP`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecordReads {
    /// The reference first, then each of the record's alleles.
    pub ad: SmallVec<u32, 3>,
    /// The informative reads, those of another allele at the site included.
    pub depth: u32,
}

impl HaplotypeScores {
    /// `None` without reads, or when no pair can explain them all: there is
    /// then nothing to call.
    #[must_use]
    pub fn new(cluster: Cluster, reads: Box<[ReadScores]>) -> Option<Self> {
        let verdict = Verdict::of(&cluster, &reads)?;
        Some(Self { cluster, reads, verdict })
    }

    /// Only reads that at least one haplotype can explain.
    #[must_use]
    pub fn reads(&self) -> &[ReadScores] {
        &self.reads
    }

    #[must_use]
    pub fn verdict(&self) -> &Verdict {
        &self.verdict
    }

    /// [`Self::new`] with every allele at column 0, for fixtures.
    #[cfg(test)]
    #[must_use]
    pub fn at_column(
        alleles: SmallVec<IndelAllele, MAX_EVENTS>,
        reads: Box<[ReadScores]>,
    ) -> Option<Self> {
        let events = alleles
            .into_iter()
            .map(|allele| Event { anchor: Pos0::ZERO, allele: allele.into() })
            .collect();
        Self::new(Cluster::new(events)?, reads)
    }

    #[must_use]
    pub fn cluster(&self) -> &Cluster {
        &self.cluster
    }

    pub fn event(&self, event: EventIdx) -> Option<&Event> {
        self.cluster.events.get(event.position())
    }

    /// The events anchored at `pos`, in record order.
    pub fn events_at(&self, pos: Pos0) -> impl Iterator<Item = EventIdx> + '_ {
        self.cluster
            .events
            .iter()
            .enumerate()
            .filter(move |(_, event)| event.anchor == pos)
            .filter_map(|(position, _)| EventIdx::new(position))
    }

    /// The indel events anchored at `pos`, in record order: the ones an
    /// indel record writes.
    pub fn indel_events_at(&self, pos: Pos0) -> impl Iterator<Item = EventIdx> + '_ {
        self.events_at(pos)
            .filter(|&e| self.event(e).is_some_and(|e| e.allele.as_indel().is_some()))
    }

    /// The substitution events anchored at `pos`, in record order.
    pub fn base_events_at(&self, pos: Pos0) -> impl Iterator<Item = EventIdx> + '_ {
        self.events_at(pos)
            .filter(|&e| self.event(e).is_some_and(|e| e.allele.as_indel().is_none()))
    }

    /// The substitution of `base` anchored at `pos`, if the cluster scored it.
    pub fn base_event_at(&self, pos: Pos0, base: Base) -> Option<EventIdx> {
        self.base_events_at(pos)
            .find(|&e| self.event(e).is_some_and(|e| e.allele == Allele::Base(base)))
    }

    /// `P(0, 1 or 2 copies of event | reads)` under a flat prior over the
    /// pairs on their own, artefact accounts aside.
    #[must_use]
    pub fn copy_posteriors(&self, event: EventIdx) -> [f64; 3] {
        let pairs = &self.verdict.pure_scores;
        let best = pairs.iter().map(|&(_, l)| l).fold(f64::NEG_INFINITY, f64::max);
        let mut by_copies = [0.0; 3];
        if !best.is_finite() {
            return by_copies;
        }
        for &(pair, log10) in pairs {
            if let Some(slot) = by_copies.get_mut(usize::from(pair.copies_of(event))) {
                *slot += 10f64.powf(log10 - best);
            }
        }
        let total: f64 = by_copies.iter().sum();
        by_copies.map(|p| p / total)
    }

    /// `read`'s best likelihood over the haplotypes without `event` and over
    /// those with it, at the alone level: the read's own say at the event's
    /// site, its other events marginalised out.
    ///
    /// For phasing votes only; a genotype read from it is blind where every
    /// read carrying the event also carries a linked one.
    #[must_use]
    pub fn marginal(&self, read: &ReadScores, event: EventIdx) -> [f32; 2] {
        [false, true].map(|carries| {
            self.cluster
                .haplotypes
                .iter()
                .filter(|hap| hap.carries(event) == carries)
                .map(|&hap| read.log10_alone(hap))
                .fold(f32::NEG_INFINITY, f32::max)
        })
    }

    /// A record's `AD` and `DP` over `events`, its alleles at one site. Each
    /// read counts for the allele its likeliest haplotype carries: none of
    /// `events` (the reference), one of them, or another allele at the site,
    /// one a chromosome cannot carry with them, which counts in `DP` only.
    /// A read counts when that allele beats every other by
    /// [`INFORMATIVE_LOG10`]; one that cannot tell them apart, as one ending
    /// before the site cannot, counts nowhere.
    #[must_use]
    pub fn reads_over(&self, events: &[EventIdx]) -> RecordReads {
        let other = events.len() + 1;
        let rival = |carried: EventIdx| {
            events.iter().any(|&event| {
                self.event(carried)
                    .zip(self.event(event))
                    .is_some_and(|(carried, event)| !carried.compatible(event))
            })
        };
        let allele = |hap: HapIdx| match events.iter().position(|&event| hap.carries(event)) {
            Some(i) => i + 1,
            None if hap.events().any(rival) => other,
            None => 0,
        };
        let alleles: SmallVec<(HapIdx, usize), MAX_HAPLOTYPES> =
            self.cluster.haplotypes.iter().map(|&hap| (hap, allele(hap))).collect();
        let mut reads = RecordReads { ad: std::iter::repeat_n(0, other).collect(), depth: 0 };
        for read in &self.reads {
            let mut best: SmallVec<f32, { MAX_EVENTS + 2 }> =
                std::iter::repeat_n(f32::NEG_INFINITY, other + 1).collect();
            for &(hap, allele) in &alleles {
                if let Some(slot) = best.get_mut(allele) {
                    *slot = slot.max(read.log10_alone(hap));
                }
            }
            let Some(allele) = informative(&best) else { continue };
            reads.depth = reads.depth.saturating_add(1);
            if let Some(count) = reads.ad.get_mut(allele) {
                *count = count.saturating_add(1);
            }
        }
        reads
    }
}

/// The position of the likeliest of `log10`, when it beats the runner-up by
/// [`INFORMATIVE_LOG10`].
fn informative(log10: &[f32]) -> Option<usize> {
    let mut best: Option<(usize, f32)> = None;
    let mut runner_up = f32::NEG_INFINITY;
    for (i, &score) in log10.iter().enumerate() {
        match best {
            Some((_, top)) if score <= top => runner_up = runner_up.max(score),
            Some((_, top)) => {
                runner_up = top;
                best = Some((i, score));
            }
            None => best = Some((i, score)),
        }
    }
    let (i, top) = best?;
    (top.is_finite() && top - runner_up >= INFORMATIVE_LOG10).then_some(i)
}

impl Verdict {
    /// The diploid pair over the cluster's haplotypes that explains `reads`
    /// best, and how sure the reads are of each of its events.
    fn of(cluster: &Cluster, reads: &[ReadScores]) -> Option<Self> {
        if reads.is_empty() {
            return None;
        }
        let haplotypes = cluster.haplotypes();
        let event_indices = || (0..cluster.events.len()).filter_map(EventIdx::new);
        let reads: Vec<(&ReadScores, Mismap)> =
            reads.iter().map(|read| (read, read.mismap(cluster))).collect();
        let reads = &reads;
        // Pairs come in GL order, each pure before its contaminated forms, so
        // ties go to the earliest and reads that say nothing make the
        // reference pair.
        let accounts: SmallVec<Account, MAX_ACCOUNTS> = HapPair::all(haplotypes)
            .flat_map(|pair| {
                let pure = Account {
                    pair,
                    artefact: None,
                    log10: reads
                        .iter()
                        .map(|(read, mismap)| mismap.mix(read.log10_unmixed(cluster, pair)))
                        .sum(),
                };
                // The reference is a contaminant too, of a pair without it: a
                // stray reference read at a homozygous alt is the same kind of
                // artefact as a stray alt read at a homozygous reference.
                let artefacts = event_indices()
                    .filter(move |event| !pair.carries(*event))
                    .map(Artefact::Event)
                    .chain((pair.low != HapIdx::REFERENCE).then_some(Artefact::Reference));
                let contaminated = artefacts.map(move |artefact| {
                    let hap = artefact.haplotype();
                    let log10: f64 = reads
                        .iter()
                        .map(|(read, mismap)| {
                            mismap.mix(read.log10_contaminated(
                                cluster,
                                pair,
                                hap,
                                ARTEFACT_FRACTION,
                            ))
                        })
                        .sum();
                    Account { pair, artefact: Some(artefact), log10: log10 - ARTEFACT_PRIOR_LOG10 }
                });
                std::iter::once(pure).chain(contaminated)
            })
            .collect();
        let pair_scores: SmallVec<(HapPair, f64), MAX_PAIRS> = HapPair::all(haplotypes)
            .map(|pair| {
                let best = accounts
                    .iter()
                    .filter(|account| account.pair == pair)
                    .map(|account| account.log10)
                    .fold(f64::NEG_INFINITY, f64::max);
                (pair, best)
            })
            .collect();
        let best = likeliest(&pair_scores, |_| true).filter(|(_, log10)| log10.is_finite())?.0;

        // Flat prior over the pairs on their own, artefact accounts aside:
        // P(not this pair) = others / (1 + others). The phasing gate on
        // `genotype_quality` was calibrated without them.
        let pure = || accounts.iter().filter(|account| account.artefact.is_none());
        let quality_of = |called: &Account| {
            let others: f64 = pure()
                .filter(|other| other.pair != called.pair)
                .map(|other| 10f64.powf(other.log10 - called.log10))
                .sum();
            phred_of_others(others)
        };

        let best_where = |carries: bool, event: EventIdx| {
            accounts
                .iter()
                .filter(|account| account.pair.carries(event) == carries)
                .map(|account| account.log10)
                .reduce(f64::max)
                .unwrap_or(f64::NEG_INFINITY)
        };
        let margins: SmallVec<Margin, MAX_EVENTS> = event_indices()
            .map(|event| Margin(best_where(true, event) - best_where(false, event)))
            .collect();

        // The called pair carries exactly the events called present; the
        // cluster's limits guarantee such a pair exists.
        let (called, _) = likeliest(&pair_scores, |pair| calls(&margins, pair))?;
        // The called pair's best account may be contaminated; its quality is
        // the pair's on its own.
        let called_pure = pure().find(|account| account.pair == called)?;
        let pure_scores = pure().map(|account| (account.pair, account.log10)).collect();
        Some(Self {
            best,
            called,
            genotype_quality: quality_of(called_pure),
            margins,
            pair_scores,
            pure_scores,
        })
    }
}

/// Pairs over [`MAX_HAPLOTYPES`] haplotypes, each on its own and with each of
/// up to [`MAX_EVENTS`] events it lacks, and the reference, as an artefact.
const MAX_ACCOUNTS: usize = MAX_PAIRS * (2 + MAX_EVENTS);

/// Unordered pairs over [`MAX_HAPLOTYPES`] haplotypes.
const MAX_PAIRS: usize = MAX_HAPLOTYPES * (MAX_HAPLOTYPES + 1) / 2;

/// The haplotype an artefact account draws its stray reads from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Artefact {
    /// The haplotype carrying this event alone, for a pair without it.
    Event(EventIdx),
    /// The reference, for a pair without it.
    Reference,
}

impl Artefact {
    fn haplotype(self) -> HapIdx {
        match self {
            Self::Event(event) => HapIdx::of_event(event),
            Self::Reference => HapIdx::REFERENCE,
        }
    }
}

/// One explanation of a cluster's reads: `pair` on its own, or with the reads
/// of an `artefact` haplotype mixed in at [`ARTEFACT_FRACTION`], its
/// prior already paid in `log10`.
struct Account {
    pair: HapPair,
    artefact: Option<Artefact>,
    log10: f64,
}

/// The likeliest of the pairs `admits` takes, the earliest of equals.
fn likeliest(
    pair_scores: &[(HapPair, f64)],
    admits: impl Fn(HapPair) -> bool,
) -> Option<(HapPair, f64)> {
    pair_scores.iter().filter(|(pair, _)| admits(*pair)).fold(None, |best, &(pair, log10)| {
        match best {
            Some((_, top)) if top >= log10 => best,
            _ => Some((pair, log10)),
        }
    })
}

/// Whether `pair` carries exactly the events `margins` call present.
fn calls(margins: &[Margin], pair: HapPair) -> bool {
    margins.iter().enumerate().all(|(i, margin)| {
        EventIdx::new(i).is_some_and(|event| margin.is_present() == pair.carries(event))
    })
}

/// How much likelier a cluster's reads are with an event present than without
/// it: `log10` of the best account whose pair carries it over the best account
/// without it, a pair without the event or such a pair with the event's reads
/// as an artefact at [`ARTEFACT_FRACTION`], whatever its zygosity. Always
/// finite: the mismapping term (`Mismap::stray`) explains every read under
/// every account.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Margin(f64);

impl Margin {
    /// The least margin an event is called present on: the reads ten times
    /// likelier with it, QUAL 10.
    pub const PRESENT: Self = Self(1.0);

    /// `log10` of how much likelier the reads are with the event.
    #[must_use]
    pub const fn log10(self) -> f64 {
        self.0
    }

    #[must_use]
    pub fn is_present(self) -> bool {
        self >= Self::PRESENT
    }

    #[must_use]
    pub fn min(self, other: Self) -> Self {
        Self(self.0.min(other.0))
    }

    /// The margin as VCF `QUAL`: ten times the margin, between 0 and 255; a
    /// likelihood ratio between two best accounts, not a posterior. Written
    /// out it is capped at 99, like every `QUAL`.
    pub fn qual(self) -> Phred {
        phred_of_log10(-self.0)
    }
}

/// A cluster's diploid genotype from the pair-HMM, and each event's
/// [`Margin`].
#[derive(Debug, Clone)]
pub struct Verdict {
    /// The pair of the account explaining the reads best: each read drawn
    /// from either of its haplotypes with equal chance, or from an artefact
    /// event's haplotype at [`ARTEFACT_FRACTION`] instead.
    pub best: HapPair,
    /// The likeliest pair carrying exactly the events called present, which
    /// is the best pair unless an event in it misses [`Margin::PRESENT`].
    called: HapPair,
    /// Phred of the chance the called pair is wrong, under a flat prior over
    /// the pairs on their own (artefact accounts aside); low whenever it is
    /// not the best pair.
    pub genotype_quality: Phred,
    /// By [`EventIdx`].
    margins: SmallVec<Margin, MAX_EVENTS>,
    /// Each pair's best account, pure or with an artefact, in `GL` order.
    pair_scores: SmallVec<(HapPair, f64), MAX_PAIRS>,
    /// Each pair's account on its own, in `GL` order.
    pure_scores: SmallVec<(HapPair, f64), MAX_PAIRS>,
}

impl Verdict {
    /// `event`'s margin; `None` for an event not scored.
    #[must_use]
    pub fn margin(&self, event: EventIdx) -> Option<Margin> {
        self.margins.get(event.position()).copied()
    }

    /// The likeliest pair over the events called present, which decides the
    /// indel calls; a record's `GT` comes from its own
    /// [`Self::genotype_likelihoods`].
    ///
    /// An event is called present only from inside the best pair (a margin
    /// above zero already implies membership), so the called pair is the best
    /// pair whenever every event in it is called. When one is not, the reads
    /// decide what remains: a pair of two haplotypes whose second event is
    /// not called becomes the homozygote of the first when no read carries
    /// the reference.
    #[must_use]
    pub fn called_pair(&self) -> HapPair {
        self.called
    }

    /// Chromosomes of the called pair carrying `event`.
    #[must_use]
    pub fn copies_of(&self, event: EventIdx) -> u8 {
        self.called.copies_of(event)
    }

    /// A record's genotype likelihoods over `events`, its alternate alleles in
    /// order: each genotype the sum over the pairs whose two haplotypes carry
    /// its alleles, every other event of the cluster carried exactly where
    /// `held` says, each pair at its best account. Summing is the flat-prior
    /// marginal over what the record does not name (the other events'
    /// zygosity and phase). A haplotype carries the record's reference
    /// allele when it carries none of `events`, and makes no genotype of the
    /// record when it carries two. `None` when no pair qualifies.
    #[must_use]
    pub fn genotype_likelihoods(
        &self,
        events: &[EventIdx],
        held: impl Fn(EventIdx) -> bool,
    ) -> Option<GenotypeLikelihoods> {
        let alleles = events.len() + 1;
        let mut log10: SmallVec<f64, 6> =
            std::iter::repeat_n(f64::NEG_INFINITY, alleles * (alleles + 1) / 2).collect();
        let allele = |hap: HapIdx| {
            let mut carried = events.iter().enumerate().filter(|&(_, &event)| hap.carries(event));
            match (carried.next(), carried.next()) {
                (None, _) => Some(0),
                (Some((i, _)), None) => Some(i + 1),
                (Some(_), Some(_)) => None,
            }
        };
        let others: SmallVec<EventIdx, MAX_EVENTS> = (0..self.margins.len())
            .filter_map(EventIdx::new)
            .filter(|event| !events.contains(event))
            .collect();
        for &(pair, score) in &self.pair_scores {
            if others.iter().any(|&event| pair.carries(event) != held(event)) {
                continue;
            }
            let (low, high) = pair.haplotypes();
            let (Some(a), Some(b)) = (allele(low), allele(high)) else { continue };
            if let Some(slot) = log10.get_mut(GenotypeLikelihoods::slot(a, b)) {
                *slot = log10_sum(*slot, score);
            }
        }
        log10.iter().any(|score| score.is_finite()).then_some(GenotypeLikelihoods(log10))
    }

    /// Whether `event` is called present: what [`Self::genotype_likelihoods`]
    /// holds an indel record's other events at.
    #[must_use]
    pub fn is_called(&self, event: EventIdx) -> bool {
        self.margin(event).is_some_and(Margin::is_present)
    }
}

/// A record's diploid genotype likelihoods over its alleles, the reference
/// first: `log10`, in VCF `GL` order (`0/0, 0/1, 1/1, 0/2, 1/2, 2/2`). A
/// genotype no pair makes is `-inf`.
#[derive(Debug, Clone, PartialEq)]
pub struct GenotypeLikelihoods(SmallVec<f64, 6>);

impl GenotypeLikelihoods {
    /// The position of the genotype `a/b` in `GL` order.
    #[must_use]
    pub fn index(a: GtAllele, b: GtAllele) -> usize {
        Self::slot(a.number().into(), b.number().into())
    }

    fn slot(a: usize, b: usize) -> usize {
        let (low, high) = (a.min(b), a.max(b));
        high * (high + 1) / 2 + low
    }

    /// `log10` of the genotype `a/b`.
    #[must_use]
    pub fn get(&self, a: GtAllele, b: GtAllele) -> Option<f64> {
        self.0.get(Self::index(a, b)).copied()
    }

    /// The record's alleles, the reference included.
    #[must_use]
    pub fn alleles(&self) -> usize {
        (1..).find(|n| n * (n + 1) / 2 >= self.0.len()).unwrap_or(1)
    }

    /// The likeliest genotype carrying an alt, the earliest in `GL` order of
    /// equals; `None` when no pair makes one. Where the alleles' presence is
    /// decided elsewhere, by the forests or the margins, this is the record's
    /// `GT`, so that `GT`, `GQ` and `PL` agree.
    #[must_use]
    pub fn likeliest_with_an_alt(&self) -> Option<GenotypeTag> {
        let allele = |i: usize| {
            Some(NonZeroU8::new(u8::try_from(i).ok()?).map_or(GtAllele::Ref, GtAllele::Alt))
        };
        (1..self.alleles())
            .flat_map(|high| (0..=high).map(move |low| (low, high)))
            .filter_map(|(low, high)| {
                let log10 = self.0.get(Self::slot(low, high)).copied()?;
                let genotype = GenotypeTag::from_alleles(allele(low)?, allele(high)?);
                log10.is_finite().then_some((genotype, log10))
            })
            .fold(None, |best: Option<(GenotypeTag, f64)>, next| match best {
                Some(top) if top.1 >= next.1 => Some(top),
                _ => Some(next),
            })
            .map(|(genotype, _)| genotype)
    }

    /// `GQ` of the genotype `a/b`: Phred of the chance it is wrong under a
    /// flat prior over the record's genotypes, `0/0` included; 0 for a
    /// genotype no pair makes.
    pub fn quality(&self, a: GtAllele, b: GtAllele) -> Phred {
        let index = Self::index(a, b);
        let Some(&called) = self.0.get(index).filter(|score| score.is_finite()) else {
            return Phred::from_phred(0);
        };
        let others: f64 = self
            .0
            .iter()
            .enumerate()
            .filter(|&(i, _)| i != index)
            .map(|(_, &score)| 10f64.powf(score - called))
            .sum();
        phred_of_others(others)
    }

    /// `PL`: each genotype's Phred below the likeliest, rounded, held at
    /// 255 (where a genotype no pair makes stands).
    #[must_use]
    pub fn phred(&self) -> SmallVec<i32, 6> {
        let best = self.0.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        self.0
            .iter()
            .map(|&score| {
                let phred = (-10.0 * (score - best)).round();
                #[expect(clippy::cast_possible_truncation, reason = "held within 0..=255 first")]
                let phred = if phred.is_nan() { 255 } else { phred.clamp(0.0, 255.0) as i32 };
                phred
            })
            .collect()
    }

    /// The same likelihoods with the alleles renumbered: allele `i` here is
    /// allele `order[i]` of `self`. `None` when `order` does not name each of
    /// `self`'s alleles once.
    #[must_use]
    pub fn reordered(&self, order: &[usize]) -> Option<Self> {
        let mut seen: SmallVec<usize, 3> = order.into();
        seen.sort_unstable();
        if !seen.iter().copied().eq(0..self.alleles()) {
            return None;
        }
        let mut log10 = SmallVec::new();
        for high in 0..order.len() {
            for low in 0..=high {
                let (&a, &b) = (order.get(low)?, order.get(high)?);
                log10.push(*self.0.get(Self::slot(a, b))?);
            }
        }
        Some(Self(log10))
    }
}

/// Phred of the chance a call is wrong under a flat prior, from `others`, the
/// other options' likelihoods summed relative to the call's own:
/// `others / (1 + others)`, and 0 where `others` is not finite.
fn phred_of_others(others: f64) -> Phred {
    if !others.is_finite() {
        return Phred::from_phred(0);
    }
    phred_of_log10(others.log10() - others.ln_1p() / std::f64::consts::LN_10)
}

/// Phred of a probability given as its `log10`, rounded and held between 0
/// and 255 like the other genotype qualities; NaN is 0.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "clamped into u8 first"
)]
fn phred_of_log10(log10_probability: f64) -> Phred {
    Phred::from_phred((-10.0 * log10_probability).round().clamp(0.0, 255.0) as u8)
}

/// `log10(10^a + 10^b)`, stable for large negative values, and `-inf` only
/// when both are.
fn log10_sum(a: f64, b: f64) -> f64 {
    let (hi, lo) = if a >= b { (a, b) } else { (b, a) };
    if hi == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    hi + (1.0 + 10f64.powf(lo - hi)).log10()
}

/// `log10((10^a + 10^b) / 2)`.
fn log10_mean(a: f64, b: f64) -> f64 {
    log10_sum(a, b) - std::f64::consts::LOG10_2
}

/// Score rows the tests share: a read each haplotype explains alone.
#[cfg(test)]
pub(crate) mod fixtures {
    use super::{ReadScores, row};
    use crate::{call::variant_calling::GtAllele, utils::MappingQuality};

    pub(crate) const NO: f32 = f32::NEG_INFINITY;
    pub(crate) const REF: GtAllele = GtAllele::Ref;
    pub(crate) const ALT: GtAllele = GtAllele::Alt(std::num::NonZeroU8::MIN);
    pub(crate) const REF_READ: [f32; 3] = [-1.0, -4.0, -4.0];
    pub(crate) const ALT1_READ: [f32; 3] = [-4.0, -1.0, -4.0];
    pub(crate) const ALT2_READ: [f32; 3] = [-4.0, -4.0, -1.0];

    /// `n` anonymous reads of each `log10` row, at MAPQ 60.
    pub(crate) fn rows(rows: &[([f32; 3], usize)]) -> Vec<ReadScores> {
        rows.iter()
            .flat_map(|&(log10, n)| {
                std::iter::repeat_n(
                    ReadScores::new(None, MappingQuality::from_byte(60), row(log10)),
                    n,
                )
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::{fixtures::*, *};
    use crate::call::test_helpers::pos0;
    use crate::utils::mapping_quality::MAPQ_UNIQUE;
    use hegel::{TestCase, generators as gs};
    use seqair_types::Base;

    fn deletion(n: usize) -> IndelAllele {
        IndelAllele::Deletion(std::iter::repeat_n(Base::T, n).collect())
    }

    fn read(log10: [f32; 3]) -> ReadScores {
        ReadScores::new(None, MappingQuality::from_byte(60), row(log10))
    }

    /// `alts` scored alleles at one anchor, over `reads` of each kind.
    fn scores(alts: usize, reads: &[([f32; 3], usize)]) -> Option<HaplotypeScores> {
        HaplotypeScores::at_column((1..=alts).map(deletion).collect(), rows(reads).into())
    }

    fn pair(a: u8, b: u8) -> HapPair {
        HapPair::new(HapIdx(a), HapIdx(b))
    }

    fn ev(position: usize) -> EventIdx {
        EventIdx::new(position).expect("an event")
    }

    /// The haplotype carrying the alternative at `position` alone.
    fn alt(position: usize) -> HapIdx {
        HapIdx::of_event(ev(position))
    }

    fn verdict(scores: Option<HaplotypeScores>) -> Verdict {
        scores.expect("a verdict").verdict
    }

    fn likelihoods(log10: &[f64]) -> GenotypeLikelihoods {
        GenotypeLikelihoods(log10.into())
    }

    #[test]
    fn genotypes_are_in_gl_order() {
        let order = [(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2)];
        for (index, (a, b)) in order.into_iter().enumerate() {
            assert_eq!(GenotypeLikelihoods::slot(a, b), index);
            assert_eq!(GenotypeLikelihoods::slot(b, a), index);
        }
        assert_eq!(likelihoods(&[0.0; 3]).alleles(), 2);
        assert_eq!(likelihoods(&[0.0; 6]).alleles(), 3);
    }

    /// `GQ` weighs every other genotype, `0/0` included; `PL` is each one's
    /// distance from the likeliest, a genotype no pair makes at 255.
    #[test]
    fn gq_and_pl_come_from_one_vector() {
        let het = likelihoods(&[-3.0, 0.0, f64::NEG_INFINITY]);
        assert_eq!(*het.quality(REF, ALT), 30.0);
        assert_eq!(*het.quality(REF, REF), 0.0);
        assert_eq!(*het.quality(ALT, ALT), 0.0, "a genotype no pair makes");
        assert_eq!(&het.phred()[..], [30, 0, 255]);
        // Two runners-up close to the call: GQ below either's PL.
        let close = likelihoods(&[-0.5, 0.0, -0.5]);
        assert!(*close.quality(REF, ALT) < 5.0, "{:?}", close.quality(REF, ALT));
        assert_eq!(&close.phred()[..], [5, 0, 5]);
        // A written genotype the reads rule out by more than `f64` spans.
        let forced = likelihoods(&[-400.0, 0.0, -400.0]);
        assert_eq!(*forced.quality(REF, REF), 0.0);
        assert_eq!(&forced.phred()[..], [255, 0, 255]);
    }

    #[test]
    fn reordering_renumbers_the_alts() {
        let gl = likelihoods(&[0.0, -1.0, -2.0, -3.0, -4.0, -5.0]);
        let swapped = gl.reordered(&[0, 2, 1]).expect("a permutation");
        // 0/0, 0/2, 2/2, 0/1, 1/2, 1/1 of the original.
        assert_eq!(swapped, likelihoods(&[0.0, -3.0, -5.0, -1.0, -4.0, -2.0]));
        assert_eq!(gl.reordered(&[0, 1]), None, "an allele left out");
        assert_eq!(gl.reordered(&[0, 1, 1]), None);
    }

    /// One event: each genotype is its one pair. Two at one anchor, the
    /// second held absent: the record over the first never sees a pair
    /// carrying the second.
    #[test]
    fn a_record_holds_the_other_events_of_its_cluster() {
        let single = verdict(scores(1, &[(REF_READ, 4), (ALT1_READ, 4)]));
        let gl = single.genotype_likelihoods(&[ev(0)], |_| false).expect("likelihoods");
        let pairs: SmallVec<f64, 6> = single.pair_scores.iter().map(|&(_, l)| l).collect();
        assert_eq!(gl, GenotypeLikelihoods(pairs));
        assert_eq!(gl.phred().get(1), Some(&0), "the het");

        let two = verdict(scores(2, &[(REF_READ, 4), (ALT1_READ, 4)]));
        let held = two.genotype_likelihoods(&[ev(0)], |_| false).expect("likelihoods");
        let carrying_neither = |pair: HapPair| !pair.carries(ev(1));
        let expected: SmallVec<f64, 6> = two
            .pair_scores
            .iter()
            .filter(|(pair, _)| carrying_neither(*pair))
            .map(|&(_, l)| l)
            .collect();
        assert_eq!(held, GenotypeLikelihoods(expected));
    }

    /// A strong reference read at MAPQ 10 costs the homozygote of the alt at
    /// most its MAPQ once mixed; unmixed it costs the full three log10.
    #[test]
    fn a_low_mapq_read_costs_a_pair_at_most_its_mapq() {
        let cluster = Cluster::new(
            [Event { anchor: pos0(0), allele: deletion(1).into() }].into_iter().collect(),
        )
        .expect("a cluster");
        let read = ReadScores::new(None, MappingQuality::from_byte(10), row([-1.0, -4.0, NO]));
        let hom = pair(1, 1);
        let unmixed = read.log10_unmixed(&cluster, hom);
        assert!((unmixed + 4.0).abs() < 1e-9);
        let mixed = read.log10_under(&cluster, hom);
        assert!((-2.0..-1.0).contains(&mixed), "{mixed}");
        let reference = read.log10_under(&cluster, pair(0, 0));
        assert!((reference + 1.0).abs() < 0.05, "{reference}");
        let unmapped = ReadScores::new(None, MappingQuality::from_byte(0), row([-1.0, -4.0, NO]));
        assert!((unmapped.log10_under(&cluster, hom) + 1.0).abs() < 1e-9);
    }

    #[test]
    fn mapq_255_is_unavailable_not_certain() {
        let cluster = Cluster::new(
            [Event { anchor: pos0(0), allele: deletion(1).into() }].into_iter().collect(),
        )
        .expect("a cluster");
        let hom = pair(1, 1);
        let unavailable = ReadScores::new(None, MappingQuality::UNAVAILABLE, row([-1.0, -4.0, NO]));
        let unique =
            ReadScores::new(None, MappingQuality::from_byte(MAPQ_UNIQUE), row([-1.0, -4.0, NO]));
        let log10 = unavailable.log10_under(&cluster, hom);
        assert!((log10 - unique.log10_under(&cluster, hom)).abs() < 1e-12);
        assert!(log10 > unavailable.log10_unmixed(&cluster, hom));
    }

    /// One stray reference read among twelve alt reads, every read Q40-like:
    /// at face value the stray read costs the homozygote more than the het's
    /// halving costs the twelve, but the homozygote with the reference as a
    /// 10 % artefact explains both.
    #[test]
    fn a_stray_reference_read_is_an_artefact_of_a_homozygote() {
        const STRAY: [f32; 3] = [0.0, -4.5, NO];
        const ALT: [f32; 3] = [-4.5, 0.0, NO];
        let scores = scores(1, &[(STRAY, 1), (ALT, 12)]).expect("scores");
        let under = |p: HapPair| -> f64 {
            scores.reads.iter().map(|read| read.log10_under(scores.cluster(), p)).sum()
        };
        assert!(under(pair(0, 1)) > under(pair(1, 1)), "pure: the het wins");
        assert_eq!(scores.verdict.best, pair(1, 1));
        assert_eq!(scores.verdict.called_pair(), pair(1, 1));
    }

    /// One copy or two by the best pair carrying the event, each pair at its
    /// best account, so the stray reference read leaves the homozygote; GQ
    /// grows with the reads and is low on a near tie.
    #[test]
    fn copies_are_the_best_carrying_pair_at_its_best_account() {
        let gq = |reads: &[([f32; 3], usize)]| {
            let v = verdict(scores(1, reads));
            let copies = v.copies_of(ev(0));
            let gl = v.genotype_likelihoods(&[ev(0)], |_| false).expect("likelihoods");
            let quality = if copies == 2 { gl.quality(ALT, ALT) } else { gl.quality(REF, ALT) };
            (copies, *quality)
        };
        assert_eq!(gq(&[([0.0, -4.5, NO], 1), ([-4.5, 0.0, NO], 12)]).0, 2);
        let (copies, het) = gq(&[([0.0, -4.5, NO], 6), ([-4.5, 0.0, NO], 6)]);
        assert_eq!(copies, 1);
        assert!(het > 30.0, "{het}");
        let (hom, deeper) = (gq(&[([-4.5, 0.0, NO], 3)]), gq(&[([-4.5, 0.0, NO], 12)]));
        assert_eq!((hom.0, deeper.0), (2, 2));
        assert!(deeper.1 > hom.1, "{hom:?} {deeper:?}");
    }

    #[test]
    fn a_clean_het_stays_het_beside_the_reference_artefact() {
        let verdict = verdict(scores(1, &[([0.0, -4.5, NO], 6), ([-4.5, 0.0, NO], 6)]));
        assert_eq!(verdict.called_pair(), pair(0, 1));
    }

    /// Two substitutions, each making its own divergent `CpG`: a haplotype is
    /// scored at its doubled level beside a partner carrying all its `CpG`s,
    /// and alone beside one carrying only some, the haplotype with both beside
    /// one with the first among them.
    #[test]
    fn a_pair_reads_the_doubled_level_only_beside_every_cpg() {
        let snv = |anchor, base| Event { anchor: pos0(anchor), allele: Allele::Base(base) };
        let mut cluster = Cluster::new([snv(10, Base::C), snv(30, Base::C)].into_iter().collect())
            .expect("a cluster");
        assert_eq!(cluster.haplotypes().len(), 4);
        cluster.divergent_cpgs = [0, 0b01, 0b10, 0b11, 0, 0, 0, 0];
        let (reference, first, second, both) = (HapIdx(0), HapIdx(1), HapIdx(2), HapIdx(3));
        assert!(cluster.shares_cpgs(first, first));
        assert!(cluster.shares_cpgs(first, both));
        assert!(!cluster.shares_cpgs(both, first), "mixed: both carries a CpG first lacks");
        assert!(!cluster.shares_cpgs(first, second));
        assert!(!cluster.shares_cpgs(first, reference));
        assert!(cluster.shares_cpgs(reference, both));

        let mut read = ReadScores::new(
            None,
            MappingQuality::from_byte(60),
            [-1.0, -2.0, -3.0, -4.0, NO, NO, NO, NO],
        );
        read.log10_doubled = [-1.0, -0.5, -1.5, -2.5, NO, NO, NO, NO];
        let beside = |h, p| read.log10_beside(&cluster, h, p);
        assert_eq!(beside(both, first), -4.0);
        assert_eq!(beside(first, both), -0.5);
        assert_eq!(beside(both, both), -2.5);
        assert_eq!(beside(second, first), -3.0);
        assert_eq!(read.log10_alone(both), -4.0);
    }

    #[test]
    fn pairs_come_in_gl_order() {
        let haplotypes = [HapIdx(0), HapIdx(1), HapIdx(2)];
        let pairs: Vec<_> = HapPair::all(&haplotypes)
            .map(|p| (p.haplotypes().0.index(), p.haplotypes().1.index()))
            .collect();
        assert_eq!(pairs, [(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2)]);
    }

    /// Two events at one anchor cannot share a chromosome; at different
    /// anchors they can, unless one deletes over the other's anchor.
    #[test]
    fn haplotypes_are_the_compatible_subsets() {
        let at =
            |anchor, allele: IndelAllele| Event { anchor: pos0(anchor), allele: allele.into() };
        let same_anchor = [at(10, deletion(1)), at(10, deletion(2))];
        assert_eq!(HapIdx::all_valid(&same_anchor).as_ref(), [HapIdx(0), HapIdx(1), HapIdx(2)]);
        let apart = [at(10, deletion(1)), at(20, deletion(2))];
        assert_eq!(
            HapIdx::all_valid(&apart).as_ref(),
            [HapIdx(0), HapIdx(1), HapIdx(2), HapIdx(3)]
        );
        let overlapping = [at(10, deletion(5)), at(12, deletion(1))];
        assert_eq!(HapIdx::all_valid(&overlapping).as_ref(), [HapIdx(0), HapIdx(1), HapIdx(2)]);
        let three = [at(10, deletion(1)), at(20, deletion(1)), at(20, deletion(2))];
        let valid: Vec<u8> = HapIdx::all_valid(&three).iter().map(|h| h.0).collect();
        assert_eq!(valid, [0, 1, 2, 3, 4, 5]);
    }

    /// A substitution shares a chromosome with any event not at its own base,
    /// and with no deletion of it: the base a deletion removes cannot also be
    /// changed. At the deletion's anchor (the base it keeps) the two are at one
    /// anchor and so apart as any two events there are.
    #[test]
    fn a_substitution_is_incompatible_inside_a_deleted_span() {
        let snv = |anchor| Event { anchor: pos0(anchor), allele: Allele::Base(Base::A) };
        let del = Event { anchor: pos0(10), allele: deletion(3).into() };
        let valid = |events: &[Event]| -> Vec<u8> {
            HapIdx::all_valid(events).iter().map(|h| h.0).collect()
        };
        for inside in 11..=13 {
            assert_eq!(valid(&[del.clone(), snv(inside)]), [0, 1, 2], "SNV at {inside}");
        }
        for outside in [9, 14, 30] {
            assert_eq!(valid(&[del.clone(), snv(outside)]), [0, 1, 2, 3], "SNV at {outside}");
        }
        assert_eq!(valid(&[del.clone(), snv(10)]), [0, 1, 2, 3], "the base it keeps");
        let ins = Event {
            anchor: pos0(12),
            allele: IndelAllele::Insertion([Base::G].into_iter().collect()).into(),
        };
        assert_eq!(valid(&[snv(12), ins.clone()]), [0, 1, 2, 3], "the base before it");
        assert_eq!(valid(&[snv(13), ins]), [0, 1, 2, 3]);
        assert_eq!(valid(&[snv(5), snv(6)]), [0, 1, 2, 3]);
    }

    /// Two substitutions at one base are a multi-allelic site: never one
    /// chromosome, and at most two of them.
    #[test]
    fn substitutions_at_one_base_follow_the_anchor_limits() {
        let snv = |base| Event { anchor: pos0(7), allele: Allele::Base(base) };
        let two: SmallVec<Event, MAX_EVENTS> = [snv(Base::A), snv(Base::C)].into_iter().collect();
        let cluster = Cluster::new(two).expect("two at one anchor");
        assert_eq!(cluster.haplotypes(), [HapIdx(0), HapIdx(1), HapIdx(2)]);
        assert!(!Cluster::admits(cluster.events(), &snv(Base::G)));
    }

    /// A third event at one anchor, or a fourth in all, makes no cluster.
    #[test]
    fn a_cluster_holds_two_events_per_anchor_and_three_in_all() {
        let at = |anchor, n| Event { anchor: pos0(anchor), allele: deletion(n).into() };
        let cluster = |events: &[Event]| Cluster::new(events.iter().cloned().collect());
        assert!(cluster(&[at(10, 1), at(10, 2), at(20, 1)]).is_some());
        assert!(cluster(&[at(10, 1), at(10, 2), at(10, 3)]).is_none());
        assert!(cluster(&[at(10, 1), at(20, 1), at(30, 1), at(40, 1)]).is_none());
        assert!(!Cluster::admits(&[at(10, 1), at(10, 2)], &at(10, 3)));
        assert!(Cluster::admits(&[at(10, 1), at(10, 2)], &at(20, 1)));
    }

    /// Two events on one chromosome: reads carrying both make the pair of
    /// that haplotype with itself, and each event two copies.
    #[test]
    fn two_events_apart_can_be_one_haplotype() {
        let at =
            |anchor, allele: IndelAllele| Event { anchor: pos0(anchor), allele: allele.into() };
        let events = [at(10, deletion(1)), at(20, deletion(2))].into_iter().collect();
        let both = [-9.0, -5.0, -5.0, -1.0, NO, NO, NO, NO];
        let reads: Vec<ReadScores> =
            std::iter::repeat_n(ReadScores::new(None, MappingQuality::from_byte(60), both), 8)
                .collect();
        let cluster = Cluster::new(events).expect("a cluster");
        let s = HaplotypeScores::new(cluster, reads.into()).expect("a verdict");
        assert_eq!(s.verdict.called_pair(), pair(3, 3));
        assert_eq!(s.verdict.copies_of(ev(0)), 2);
        assert_eq!(s.verdict.copies_of(ev(1)), 2);
        assert!(s.verdict.margin(ev(0)).expect("scored").is_present());
        assert_eq!(s.reads_over(&[ev(1)]).ad.as_ref(), [0, 8]);
    }

    /// Reads of a second, uncalled deletion at the anchor are neither the
    /// reference nor the called deletion, so the called record's `AD` leaves
    /// them out; the record writing both alleles counts them.
    #[test]
    fn reads_of_another_allele_at_the_anchor_are_not_ref_support() {
        let s = scores(
            2,
            &[([-1.0, -9.0, -9.0], 10), ([-9.0, -1.0, -9.0], 10), ([-9.0, -9.0, -1.0], 4)],
        )
        .expect("a verdict");
        assert_eq!(s.verdict.copies_of(ev(1)), 0);
        assert_eq!(s.reads_over(&[ev(0)]), RecordReads { ad: [10, 10].into(), depth: 24 });
        assert_eq!(
            s.reads_over(&[ev(0), ev(1)]),
            RecordReads { ad: [10, 10, 4].into(), depth: 24 }
        );
    }

    /// The carrier of an event on one chromosome comes first, whatever else
    /// the two carry; an event on both or neither does not split the pair.
    #[test]
    fn a_pair_splits_on_an_event_one_chromosome_carries() {
        assert_eq!(pair(0, 3).split(ev(1)), Some((HapIdx(3), HapIdx(0))));
        assert_eq!(pair(1, 2).split(ev(1)), Some((HapIdx(2), HapIdx(1))));
        assert_eq!(pair(1, 2).split(ev(0)), Some((HapIdx(1), HapIdx(2))));
        assert_eq!(pair(3, 3).split(ev(0)), None);
        assert_eq!(pair(0, 0).split(ev(0)), None);
    }

    #[test]
    fn the_mean_of_two_likelihoods_is_stable_and_total() {
        assert!((log10_mean(-1.0, -1.0) - -1.0).abs() < 1e-12);
        assert!(
            (log10_mean(-400.0, f64::NEG_INFINITY) - (-400.0 - std::f64::consts::LOG10_2)).abs()
                < 1e-9
        );
        assert_eq!(log10_mean(f64::NEG_INFINITY, f64::NEG_INFINITY), f64::NEG_INFINITY);
    }

    #[test]
    fn the_ratio_is_signed_infinite_or_absent() {
        let r = read([-1.0, -3.0, NO]);
        assert_eq!(r.log10_ratio(HapIdx::REFERENCE, alt(0)), Some(2.0));
        assert_eq!(r.log10_ratio(alt(0), HapIdx::REFERENCE), Some(-2.0));
        assert_eq!(r.log10_ratio(HapIdx::REFERENCE, alt(1)), Some(f32::INFINITY));
        assert_eq!(read([NO, NO, NO]).log10_ratio(HapIdx::REFERENCE, alt(0)), None);
        assert!(EventIdx::new(MAX_EVENTS).is_none());
    }

    #[test]
    fn reads_split_between_two_haplotypes_are_a_het() {
        let v = verdict(scores(1, &[(REF_READ, 4), (ALT1_READ, 4)]));
        assert_eq!(v.best, pair(0, 1));
        assert_eq!(v.called_pair(), pair(0, 1));
        assert!(v.genotype_quality.as_int() > 10, "{:?}", v.genotype_quality);
    }

    #[test]
    fn reads_of_one_alt_are_a_homozygote_and_of_two_a_compound_het() {
        assert_eq!(verdict(scores(1, &[(ALT1_READ, 8)])).called_pair(), pair(1, 1));
        let compound = verdict(scores(2, &[(ALT1_READ, 4), (ALT2_READ, 4)]));
        assert_eq!(compound.called_pair(), pair(1, 2));
    }

    /// A record's likelihoods do not depend on where its event sits among the
    /// cluster's: three substitutions, two in cis and one in trans, scored in
    /// position order and reversed, give each the same `GL`.
    #[test]
    fn a_record_s_likelihoods_do_not_depend_on_its_neighbours_side() {
        let anchors = [10, 18, 26];
        // Each read kind by the anchors it carries, and how many of it.
        let kinds: [(&[u32], usize); 3] = [(&[10, 18], 12), (&[26], 9), (&[], 5)];
        let likelihoods = |order: [u32; 3]| -> Vec<GenotypeLikelihoods> {
            let events = order
                .iter()
                .map(|&anchor| Event { anchor: pos0(anchor), allele: Allele::Base(Base::C) })
                .collect();
            let carried = |hap: u8, anchor: u32| {
                order.iter().position(|&a| a == anchor).is_some_and(|i| hap & (1 << i) != 0)
            };
            let reads: Box<[ReadScores]> = kinds
                .iter()
                .flat_map(|&(own, n)| {
                    let mut log10 = [NO; MAX_HAPLOTYPES];
                    for (hap, slot) in (0u8..).zip(&mut log10) {
                        let wrong = anchors
                            .iter()
                            .filter(|&&a| carried(hap, a) != own.contains(&a))
                            .count();
                        *slot = -1.0 - 3.0 * wrong as f32;
                    }
                    std::iter::repeat_n(
                        ReadScores::new(None, MappingQuality::from_byte(60), log10),
                        n,
                    )
                })
                .collect();
            let cluster = Cluster::new(events).expect("a cluster");
            let v = verdict(HaplotypeScores::new(cluster, reads));
            anchors
                .iter()
                .map(|anchor| {
                    let event = ev(order.iter().position(|a| a == anchor).expect("scored"));
                    v.genotype_likelihoods(&[event], |other| v.is_called(other))
                        .expect("likelihoods")
                })
                .collect()
        };
        let forward = likelihoods(anchors);
        let reversed = likelihoods([26, 18, 10]);
        for (a, b) in forward.iter().zip(&reversed) {
            assert_eq!(a.phred(), b.phred());
            for (x, y) in a.0.iter().zip(&b.0) {
                assert!((x - y).abs() < 1e-9, "{a:?} {b:?}");
            }
        }
    }

    /// An absent third haplotype is never part of a pair, so its `-inf` costs
    /// nothing.
    #[test]
    fn a_missing_haplotype_does_not_take_part() {
        assert!(scores(1, &[([-1.0, -1.0, NO], 1)]).is_some());
    }

    /// Reads that say nothing leave the reference as the best pair, with
    /// five of six pairs as good: Phred 1.
    #[test]
    fn a_tie_goes_to_the_reference() {
        let v = verdict(scores(2, &[([-1.0, -1.0, -1.0], 3)]));
        assert_eq!(v.best, pair(0, 0));
        assert_eq!(v.genotype_quality.as_int(), 1);
    }

    #[test]
    fn no_reads_no_verdict() {
        assert!(scores(1, &[]).is_none());
        assert!(scores(1, &[([NO, NO, NO], 1)]).is_none());
    }

    /// The margin says how sure the reads are an allele is there at all, so a
    /// het and a hom-alt both clear it and the reference does not.
    #[test]
    fn the_margin_measures_presence_whatever_the_zygosity() {
        let het = verdict(scores(1, &[(REF_READ, 2), (ALT1_READ, 2)]));
        let hom = verdict(scores(1, &[(ALT1_READ, 4)]));
        let reference = verdict(scores(1, &[(REF_READ, 4)]));
        for present in [&het, &hom] {
            assert!(present.margin(ev(0)).expect("scored").is_present());
        }
        let absent = reference.margin(ev(0)).expect("scored");
        assert!(!absent.is_present() && absent < Margin(-0.5), "{absent:?}");
        assert_eq!(reference.called_pair(), pair(0, 0));
        assert!(reference.margin(ev(2)).is_none());
    }

    /// A read only the alternative explains rules out every pair without it,
    /// but not the artefact account, so the margin is strong and finite.
    #[test]
    fn an_impossible_reference_leaves_a_finite_margin() {
        let v = verdict(scores(1, &[([NO, -1.0, NO], 1), ([-1.0, -1.0, NO], 1)]));
        let margin = v.margin(ev(0)).expect("scored");
        assert!(margin.is_present() && margin.0.is_finite(), "{margin:?}");
        assert!(*margin.qual() < 255.0);
    }

    /// One read the alternative explains about three times better than the
    /// reference: the best pair is the homozygote, but the allele is not
    /// called, and the genotype falls back to the reference.
    #[test]
    fn a_best_pair_allele_below_the_threshold_is_not_called() {
        let v = verdict(scores(1, &[([-1.5, -1.0, NO], 1)]));
        let margin = v.margin(ev(0)).expect("scored");
        assert!(margin > Margin(0.0) && !margin.is_present(), "{margin:?}");
        assert_eq!(v.best, pair(1, 1));
        assert_eq!(v.called_pair(), pair(0, 0));
    }

    /// A compound het whose second allele the reads barely support keeps the
    /// first, as the reads have it without the second: a homozygote when no
    /// read carries the reference, a het when one does. Either way the
    /// pair's quality says the written pair is not the one the reads like
    /// best.
    #[test]
    fn an_uncalled_half_of_a_compound_leaves_the_other_as_the_reads_have_it() {
        let hom = verdict(scores(2, &[([-9.0, -1.0, -9.0], 2), ([-9.0, -1.5, -1.0], 4)]));
        assert_eq!(hom.best, pair(1, 2));
        assert!(hom.margin(ev(0)).expect("scored").is_present());
        let second = hom.margin(ev(1)).expect("scored");
        assert!(second > Margin(0.0) && !second.is_present(), "{second:?}");
        assert_eq!(hom.called_pair(), pair(1, 1));
        assert!(hom.genotype_quality.as_int() < 3, "{:?}", hom.genotype_quality);

        let het = verdict(scores(
            2,
            &[([-1.0, -9.0, -2.0], 1), ([-9.0, -1.0, -9.0], 2), ([-9.0, -1.5, -1.0], 2)],
        ));
        assert_eq!(het.best, pair(1, 2));
        assert!(!het.margin(ev(1)).expect("scored").is_present());
        assert_eq!(het.called_pair(), pair(0, 1));
        assert!(het.genotype_quality.as_int() < 3, "{:?}", het.genotype_quality);
    }

    /// A read each haplotype explains alone, and the other not at all.
    const ONLY_REF: [f32; 3] = [-1.0, -9.0, -9.0];
    const ONLY_ALT1: [f32; 3] = [-9.0, -1.0, -9.0];
    const ONLY_ALT2: [f32; 3] = [-9.0, -9.0, -1.0];

    /// Seven reads of an allele among seventy-three are what an artefact
    /// supplies, not a genotype: the reads are likelier as the reference with
    /// a tenth of them from the allele than as a het, so it is not called.
    /// Four in eight are a het.
    #[test]
    fn a_low_fraction_allele_at_depth_is_an_artefact() {
        let deep = verdict(scores(1, &[(ONLY_REF, 66), (ONLY_ALT1, 7)]));
        assert_eq!(deep.best, pair(0, 0), "the best account is the reference plus an artefact");
        let margin = deep.margin(ev(0)).expect("scored");
        assert!(margin < Margin(0.0), "{margin:?}");
        assert_eq!(deep.called_pair(), pair(0, 0));

        let shallow = verdict(scores(1, &[(ONLY_REF, 4), (ONLY_ALT1, 4)]));
        assert!(shallow.margin(ev(0)).expect("scored").is_present());
        assert_eq!(shallow.called_pair(), pair(0, 1));
    }

    /// An artefact is tried on top of every genotype without it, so a het of
    /// another allele does not shelter it.
    #[test]
    fn an_artefact_beside_a_real_allele_is_not_called() {
        let v = verdict(scores(2, &[(ONLY_REF, 33), (ONLY_ALT2, 33), (ONLY_ALT1, 7)]));
        assert!(!v.margin(ev(0)).expect("scored").is_present());
        assert!(v.margin(ev(1)).expect("scored").is_present());
        assert_eq!(v.called_pair(), pair(0, 2));
    }

    #[test]
    fn a_read_prefers_a_haplotype_ten_times_likelier() {
        let r = read([-1.0, -2.0, -1.5]);
        assert!(r.prefers(HapIdx::REFERENCE, alt(0)));
        assert!(!r.prefers(HapIdx::REFERENCE, alt(1)), "only about three times likelier");
        assert!(!r.prefers(alt(0), HapIdx::REFERENCE));
        assert!(read([NO, -1.0, NO]).prefers(alt(0), HapIdx::REFERENCE));
        assert!(!read([NO, NO, NO]).prefers(alt(0), HapIdx::REFERENCE));
    }

    /// A read counts for its likeliest allele when that is informative; a
    /// read torn between two counts for neither, and not for the depth.
    #[test]
    fn a_read_counts_for_its_informative_allele() {
        let torn = ([-4.0, -1.0, -1.1], 3);
        let leaning = ([-4.0, -1.0, -1.3], 1);
        let s = scores(2, &[(REF_READ, 4), (ALT1_READ, 2), (ALT2_READ, 5), torn, leaning])
            .expect("a verdict");
        assert_eq!(s.reads_over(&[ev(0), ev(1)]), RecordReads { ad: [4, 3, 5].into(), depth: 12 });

        // With one allele scored, the third haplotype does not compete.
        let one = scores(1, &[(REF_READ, 1), ([-4.0, -1.0, 0.0], 2)]).expect("a verdict");
        assert_eq!(one.reads_over(&[ev(0)]), RecordReads { ad: [1, 2].into(), depth: 3 });
    }

    /// A read over one of two linked events scores the haplotypes with and
    /// without the other alike: it counts for the event it covers, and for
    /// neither record of the other.
    #[test]
    fn a_read_over_one_event_counts_for_its_record() {
        let at =
            |anchor, allele: IndelAllele| Event { anchor: pos0(anchor), allele: allele.into() };
        let events = [at(10, deletion(1)), at(40, deletion(2))].into_iter().collect();
        // By haplotype: reference, first, second, both.
        let read = |log10: [f32; 4]| {
            let [a, b, c, d] = log10;
            ReadScores::new(None, MappingQuality::from_byte(60), [a, b, c, d, NO, NO, NO, NO])
        };
        let reads: Vec<ReadScores> = [
            std::iter::repeat_n(read([-9.0, -1.0, -9.0, -1.0]), 6),
            std::iter::repeat_n(read([-9.0, -9.0, -1.0, -1.0]), 5),
            std::iter::repeat_n(read([-9.0, -5.0, -5.0, -1.0]), 2),
        ]
        .into_iter()
        .flatten()
        .collect();
        let s = HaplotypeScores::new(Cluster::new(events).expect("a cluster"), reads.into())
            .expect("a verdict");
        assert_eq!(s.reads_over(&[ev(0)]), RecordReads { ad: [0, 8].into(), depth: 8 });
        assert_eq!(s.reads_over(&[ev(1)]), RecordReads { ad: [0, 7].into(), depth: 7 });
    }

    #[test]
    fn log10_sum_is_stable_and_total() {
        assert!((log10_sum(-1.0, -1.0) - (-1.0 + std::f64::consts::LOG10_2)).abs() < 1e-12);
        assert!((log10_sum(-400.0, f64::NEG_INFINITY) - -400.0).abs() < 1e-9);
        assert_eq!(log10_sum(f64::NEG_INFINITY, f64::NEG_INFINITY), f64::NEG_INFINITY);
    }

    /// The best pair is called, at its own quality.
    #[test]
    fn a_called_best_pair_keeps_its_quality() {
        let v = verdict(scores(1, &[(REF_READ, 4), (ALT1_READ, 4)]));
        assert_eq!(v.called_pair(), v.best);
        assert!(v.genotype_quality.as_int() > 10, "{:?}", v.genotype_quality);
    }

    #[test]
    fn qual_is_ten_times_the_margin_held_to_a_byte() {
        let qual = |m: f64| *Margin(m).qual();
        assert_eq!(qual(1.0), 10.0);
        assert_eq!(qual(2.34), 23.0);
        assert_eq!(qual(-3.0), 0.0);
        assert_eq!(qual(f64::NEG_INFINITY), 0.0);
        assert_eq!(qual(1e9), 255.0);
    }

    /// Presence is decided inside the best pair, and a margin above zero
    /// already implies it: the best pair carrying an allele beats every
    /// pair without it, so it is the best pair overall.
    #[hegel::test]
    fn a_present_allele_is_always_in_the_called_genotype(tc: TestCase) {
        let reads = tc.draw(gs::integers::<usize>().min_value(1).max_value(11));
        let rows: Vec<[f32; 3]> = tc.draw(
            gs::vecs(gs::arrays(gs::floats::<f32>().min_value(-30.0).max_value_exclusive(0.0)))
                .min_size(reads)
                .max_size(reads),
        );
        let alts = tc.draw(gs::integers::<usize>().min_value(1).max_value(2));
        let v = verdict(HaplotypeScores::at_column(
            (1..=alts).map(deletion).collect(),
            rows.into_iter().map(read).collect(),
        ));
        for event in (0..alts).filter_map(EventIdx::new) {
            let margin = v.margin(event).expect("scored");
            if margin > Margin(0.0) {
                assert!(v.best.carries(event), "{v:?}");
            }
            assert_eq!(margin.is_present(), v.called_pair().carries(event), "{v:?}");
        }
    }
}
