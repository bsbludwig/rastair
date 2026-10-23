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
//! (`Event::compatible`: not two at one anchor, no deletion over another
//! event's anchor), so at most eight haplotypes and alignments per read. The
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
//! `Artefact` mixed in at `ARTEFACT_FRACTIONS` less
//! `ARTEFACT_PRIOR_LOG10`: an event the pair lacks, or the reference for a
//! pair without it. That keeps 7 alt reads in 73 from being a genotype while 4
//! in 8 stay a het, and one stray reference read from making a homozygote a
//! het. Margins are per event, the best account whose pair carries it over the
//! best whose pair does not; an event is called at [`Margin::PRESENT`] and
//! QUAL is ten times the margin. [`Verdict::called_pair`], the likeliest pair
//! carrying exactly the called events, is the cluster's one genotype, so
//! "present but not in GT" cannot happen while no anchor has more than two
//! events. `quality` is the called pair's flat-prior quality over the pure
//! pairs, low when it is not the best pair; the phasing gate reads it, and it
//! is the `GC` of a compound and of a `0/0` record.
//!
//! **A read's likelihood under a pair** ([`ReadScores`]):
//!
//! - Each read is mixed with its mismapping chance `10^(-MAPQ/10)`, explained
//!   by its best haplotype (`ReadScores::mapq_mixed`), so a stray low-MAPQ
//!   read costs any pair at most `MAPQ / 10` log10.
//! - Every CpG is scored at `UNKNOWN_BETA` (½), which cancels wherever both
//!   haplotypes of a pair carry it. In an all-substitution cluster a CpG one
//!   haplotype has and another lacks (a de-novo `T>C`, a `C>T` at a reference
//!   CpG) is profiled instead: its carrier is scored at the converting strand's
//!   level as the pair's only carrier ([`ReadScores::log10_alone`]) and as both
//!   chromosomes (`log10_doubled`), and a pair reads the doubled score when the
//!   partner carries all the haplotype's divergent CpGs. At ½ every converted
//!   OT read of a methylated hom-alt de-novo CpG costs the homozygote.
//!   Shortcuts: clusters with an indel stay at ½; `log10_alone` is what
//!   assignment and phasing votes read, having no pair; a haplotype whose
//!   partner carries only some of its CpGs is scored alone at all of them;
//!   `Partner::at` judges the partner from the non-carriers' consensus, else
//!   the reference.
//!
//! **`AD` and `DP`** come from the same read set: a read is assigned to its
//! likeliest haplotype when that beats the runner-up by `INFORMATIVE_LOG10`
//! (GATK's informative-read threshold). It is far below `Margin::PRESENT` on
//! purpose: the pair-HMM sums a tract gap over every placement, so one read of
//! a 1 bp deletion in a 13-base homopolymer is only ~0.9 log10 likelier under
//! the deletion. [`HaplotypeScores::support_over`] counts a record's events;
//! `AD` sums to `DP`.
//!
//! **Dosage.** [`Verdict::dosage`] projects the called pair onto one event: its
//! copies, and a GQ of Phred(1 − P(those copies)) under a flat prior over the
//! pairs consistent with the presence calls, each at its best account; `None`
//! where the called pair lacks the event. Every indel record writes that
//! genotype. An SNV's presence is the forests' call, so its record writes
//! [`Verdict::dosage_given`]: the same projection over the pairs that carry
//! at least one of the column's called alts and none of its others. It carries no genotype prior (finer than the likelihood's
//! tie margins are calibrated), and it is not a per-site marginal: that one's
//! "without" is a free best-of-N over the cluster's other alleles, so at a hom
//! indel beside a stutter or twin-spelling candidate every read is neutral and
//! one copy wins. Where no read without a linked het deletion reaches an SNV,
//! `(ref, del+SNV)` and `(SNV, del+SNV)` tie and the earlier in GL order, the
//! het, is written at GQ 5 or less: the model's honest tie, as is such a call
//! on MAPQ 0-19 reads.
//!
//! **TAPS twin spellings.** A deletion next to a methylated CpG has two
//! spellings whose haplotypes differ by a `C`/`T` at the CpG; converted OT
//! reads match the wrong one for free, so at ½ it sometimes collects them. A
//! cluster calls one of the two; which one is still open.

use crate::call::{phasing::FragmentId, pileup::indels::IndelAllele};
use seqair_types::{Base, Phred, SmallVec};

mod scoring;
pub use scoring::{HaplotypeRead, HaplotypeScratch, Prepared, prepare, score};

/// Indel events a cluster scores at most: every subset of them is a
/// haplotype, so this bounds the alignments per read at `2^MAX_EVENTS`.
pub const MAX_EVENTS: usize = 3;
/// One haplotype per subset of the events.
pub const MAX_HAPLOTYPES: usize = 1 << MAX_EVENTS;
/// Events one anchor contributes at most: two at one anchor never share a
/// chromosome, so a third could not join a pair with the other two.
const MAX_PER_ANCHOR: usize = 2;

/// A three-wide score row, reference and two alternatives at one anchor,
/// widened to [`MAX_HAPLOTYPES`]: with both events at one anchor the
/// haplotypes are `{}`, `{0}` and `{1}`, at indices 0, 1 and 2.
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
    pub fn shift(&self) -> usize {
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
    pub anchor: u32,
    pub allele: Allele,
}

impl Event {
    /// Whether one chromosome can carry both: not at one anchor, and neither
    /// a deletion over the other's anchor, which for a substitution is the
    /// base it changes.
    fn compatible(&self, other: &Self) -> bool {
        let deletes_over = |a: &Self, b: &Self| match &a.allele {
            Allele::Indel(IndelAllele::Deletion(deleted)) => {
                let end = u64::from(a.anchor).saturating_add(deleted.len() as u64);
                u64::from(b.anchor) > u64::from(a.anchor) && u64::from(b.anchor) <= end
            }
            Allele::Indel(IndelAllele::Insertion(_)) | Allele::Base(_) => false,
        };
        self.anchor != other.anchor && !deletes_over(self, other) && !deletes_over(other, self)
    }
}

/// The indel events scored together, and the haplotypes they make.
///
/// At most [`MAX_EVENTS`] events and `MAX_PER_ANCHOR` per anchor: with no
/// more than two at an anchor, any set of events at distinct anchors fits on
/// one haplotype and a second event at an anchor goes on the other, which is
/// what lets [`Verdict`] always find a pair carrying exactly the events it
/// calls present.
#[derive(Debug, Clone)]
pub struct Cluster {
    /// In record order: most reads first, then anchor, then the allele, as
    /// `cluster_candidates` ranks them. A column's record covers the events
    /// anchored at it.
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
    /// The one event of a biallelic site.
    pub const FIRST: Self = Self(0);

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
    pub mapq: u8,
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
    pub fn new(fragment: Option<FragmentId>, mapq: u8, log10: [f32; MAX_HAPLOTYPES]) -> Self {
        Self { fragment, mapq, log10_alone: log10, log10_doubled: log10 }
    }

    /// The likelihood under `hap`, its divergent `CpG`s at the level profiled
    /// as the pair's only carrier. shortcut: whatever the partner, so
    /// assignment (AD/DP) and phasing votes, which have no pair, read this.
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
    #[must_use]
    pub fn log10_ratio(&self, a: HapIdx, b: HapIdx) -> Option<f32> {
        let (a, b) = (self.log10_alone(a), self.log10_alone(b));
        (a.is_finite() || b.is_finite()).then_some(a - b)
    }

    /// Whether the read is at least [`Margin::PRESENT`] likelier under `hap`
    /// than under `over`: the same "ten times likelier" an event is called on.
    #[must_use]
    pub fn prefers(&self, hap: HapIdx, over: HapIdx) -> bool {
        self.log10_ratio(hap, over).is_some_and(|ratio| f64::from(ratio) >= Margin::PRESENT.0)
    }

    /// `log10 P(read | pair)`, the read drawn from either haplotype of the pair
    /// with equal chance.
    fn log10_under(&self, cluster: &Cluster, pair: HapPair) -> f64 {
        self.mapq_mixed(cluster, self.log10_unmixed(cluster, pair))
    }

    fn log10_unmixed(&self, cluster: &Cluster, pair: HapPair) -> f64 {
        let (a, b) = pair.haplotypes();
        log10_mean(
            f64::from(self.log10_beside(cluster, a, b)),
            f64::from(self.log10_beside(cluster, b, a)),
        )
    }

    /// `log10` under a genotype, mixed with the read's mismapping chance
    /// `10^(-MAPQ/10)`: a mismapped read is explained by whichever haplotype
    /// fits it best, so it costs no genotype more than `MAPQ / 10` log10.
    fn mapq_mixed(&self, cluster: &Cluster, log10: f64) -> f64 {
        let mismapped = 10f64.powf(-f64::from(self.mapq) / 10.0);
        let best = cluster
            .haplotypes()
            .iter()
            .map(|&hap| f64::from(self.log10_alone(hap).max(self.log10_doubled(hap))))
            .fold(f64::NEG_INFINITY, f64::max);
        log10_sum((1.0 - mismapped).log10() + log10, mismapped.log10() + best)
    }

    /// `log10 P(read | pair contaminated by hap)`: the read drawn from `hap`
    /// with chance `fraction`, else from the pair.
    fn log10_under_contaminated(
        &self,
        cluster: &Cluster,
        pair: HapPair,
        hap: HapIdx,
        fraction: f64,
    ) -> f64 {
        let contaminated = log10_sum(
            self.log10_unmixed(cluster, pair) + (1.0 - fraction).log10(),
            f64::from(self.log10_alone(hap)) + fraction.log10(),
        );
        self.mapq_mixed(cluster, contaminated)
    }
}

/// The fractions at which an event is tried as an artefact: reads drawn
/// from its haplotype this often on top of a genotype without it, as
/// stutter, a mis-mapped paralog or a sample contaminant would supply them. A
/// true event has to beat the likeliest of these by [`Margin::PRESENT`] on
/// top of [`ARTEFACT_PRIOR_LOG10`], which is what keeps seven alt reads in
/// seventy-three from being a genotype while four in eight still are. One
/// fraction: on chr12 (26x and 73x) adding 0.02 or 0.05 beside it changed no
/// call, and 0.2 cost insertions at 26x.
const ARTEFACT_FRACTIONS: [f64; 1] = [0.1];

/// How much likelier a read must be under its best haplotype than under the
/// runner-up to be assigned to it for `AD` and `DP` (GATK's informative-read
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
    /// Only reads that at least one haplotype can explain.
    pub reads: Box<[ReadScores]>,
    pub verdict: Verdict,
    /// Reads assigned to each haplotype, by [`HapIdx::index`].
    support: [u32; MAX_HAPLOTYPES],
}

impl HaplotypeScores {
    /// `None` without reads, or when no pair can explain them all: there is
    /// then nothing to call.
    #[must_use]
    pub fn new(cluster: Cluster, reads: Box<[ReadScores]>) -> Option<Self> {
        let verdict = Verdict::of(&cluster, &reads)?;
        let mut support = [0u32; MAX_HAPLOTYPES];
        for hap in reads.iter().filter_map(|read| assignment(&cluster.haplotypes, read)) {
            if let Some(count) = support.get_mut(hap.index()) {
                *count = count.saturating_add(1);
            }
        }
        Some(Self { cluster, reads, verdict, support })
    }

    /// [`Self::new`] with every allele at column 0, for fixtures.
    #[cfg(test)]
    #[must_use]
    pub fn at_column(
        alleles: SmallVec<IndelAllele, MAX_EVENTS>,
        reads: Box<[ReadScores]>,
    ) -> Option<Self> {
        let events =
            alleles.into_iter().map(|allele| Event { anchor: 0, allele: allele.into() }).collect();
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
    pub fn events_at(&self, pos: u32) -> impl Iterator<Item = EventIdx> + '_ {
        self.cluster
            .events
            .iter()
            .enumerate()
            .filter(move |(_, event)| event.anchor == pos)
            .filter_map(|(position, _)| EventIdx::new(position))
    }

    /// The indel events anchored at `pos`, in record order: the ones an
    /// indel record writes.
    pub fn indel_events_at(&self, pos: u32) -> impl Iterator<Item = EventIdx> + '_ {
        self.events_at(pos)
            .filter(|&e| self.event(e).is_some_and(|e| e.allele.as_indel().is_some()))
    }

    /// The substitution events anchored at `pos`, in record order.
    pub fn base_events_at(&self, pos: u32) -> impl Iterator<Item = EventIdx> + '_ {
        self.events_at(pos)
            .filter(|&e| self.event(e).is_some_and(|e| e.allele.as_indel().is_none()))
    }

    /// The substitution of `base` anchored at `pos`, if the cluster scored it.
    pub fn base_event_at(&self, pos: u32, base: Base) -> Option<EventIdx> {
        self.base_events_at(pos)
            .find(|&e| self.event(e).is_some_and(|e| e.allele == Allele::Base(base)))
    }

    /// `P(0, 1 or 2 copies of event | reads)` under a flat prior over the
    /// pairs on their own, artefact accounts aside.
    #[must_use]
    pub fn copy_posteriors(&self, event: EventIdx) -> [f64; 3] {
        let pairs: SmallVec<(HapPair, f64), 36> = HapPair::all(&self.cluster.haplotypes)
            .map(|pair| {
                (pair, self.reads.iter().map(|read| read.log10_under(&self.cluster, pair)).sum())
            })
            .collect();
        let best = pairs.iter().map(|&(_, l)| l).fold(f64::NEG_INFINITY, f64::max);
        let mut by_copies = [0.0; 3];
        if !best.is_finite() {
            return by_copies;
        }
        for (pair, log10) in pairs {
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

    /// Reads assigned to `hap`: its `AD`. Reads assigned to no haplotype are
    /// left out of every count, so the supports add up to [`Self::depth`].
    #[must_use]
    pub fn support(&self, hap: HapIdx) -> u32 {
        self.support.get(hap.index()).copied().unwrap_or(0)
    }

    /// A record's `AD` over `events`: reads assigned to a haplotype carrying
    /// none of them, then to one carrying each, by the first that matches.
    #[must_use]
    pub fn support_over(&self, events: &[EventIdx]) -> SmallVec<u32, 3> {
        let mut counts: SmallVec<u32, 3> = std::iter::repeat_n(0, events.len() + 1).collect();
        for &hap in &self.cluster.haplotypes {
            let slot = events.iter().position(|&event| hap.carries(event)).map_or(0, |i| i + 1);
            if let Some(count) = counts.get_mut(slot) {
                *count = count.saturating_add(self.support(hap));
            }
        }
        counts
    }

    /// Reads assigned to some haplotype: the column's `DP`, the reads that
    /// say anything about the indels.
    #[must_use]
    pub fn depth(&self) -> u32 {
        self.support.iter().fold(0, |depth, &count| depth.saturating_add(count))
    }
}

/// The haplotype a read is assigned to: the likeliest, when it beats the
/// runner-up by [`INFORMATIVE_LOG10`]; `None` for a read that cannot tell
/// them apart, as one ending before the indel cannot.
fn assignment(haplotypes: &[HapIdx], read: &ReadScores) -> Option<HapIdx> {
    let mut best: Option<(HapIdx, f32)> = None;
    let mut runner_up = f32::NEG_INFINITY;
    for &hap in haplotypes {
        let log10 = read.log10_alone(hap);
        match best {
            Some((_, top)) if log10 <= top => runner_up = runner_up.max(log10),
            Some((_, top)) => {
                runner_up = top;
                best = Some((hap, log10));
            }
            None => best = Some((hap, log10)),
        }
    }
    let (hap, top) = best?;
    (top.is_finite() && top - runner_up >= INFORMATIVE_LOG10).then_some(hap)
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
        // Pairs come in GL order, each pure before its contaminated forms, so
        // ties go to the earliest and reads that say nothing make the
        // reference pair.
        let accounts: SmallVec<Account, MAX_ACCOUNTS> = HapPair::all(haplotypes)
            .flat_map(|pair| {
                let pure = Account {
                    pair,
                    artefact: None,
                    log10: reads.iter().map(|read| read.log10_under(cluster, pair)).sum(),
                };
                // The reference is a contaminant too, of a pair without it: a
                // stray reference read at a homozygous alt is the same kind of
                // artefact as a stray alt read at a homozygous reference.
                let artefacts = event_indices()
                    .filter(move |event| !pair.carries(*event))
                    .map(Artefact::Event)
                    .chain((pair.low != HapIdx::REFERENCE).then_some(Artefact::Reference));
                let contaminated = artefacts.flat_map(move |artefact| {
                    let hap = artefact.haplotype();
                    ARTEFACT_FRACTIONS.iter().map(move |&fraction| {
                        let log10: f64 = reads
                            .iter()
                            .map(|read| read.log10_under_contaminated(cluster, pair, hap, fraction))
                            .sum();
                        Account {
                            pair,
                            artefact: Some(artefact),
                            log10: log10 - ARTEFACT_PRIOR_LOG10,
                        }
                    })
                });
                std::iter::once(pure).chain(contaminated)
            })
            .collect();
        let best = earliest_max(accounts.iter()).filter(|account| account.log10.is_finite())?.pair;

        // Flat prior over the pairs on their own, artefact accounts aside:
        // P(not this pair) = others / (1 + others). The phasing gate on
        // `quality` was set on this, so artefacts do not enter it yet.
        let pure = || accounts.iter().filter(|account| account.artefact.is_none());
        let quality_of = |called: &Account| {
            let others: f64 = pure()
                .filter(|other| other.pair != called.pair)
                .map(|other| 10f64.powf(other.log10 - called.log10))
                .sum();
            phred_of_log10(others.log10() - others.ln_1p() / std::f64::consts::LN_10)
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
        let present =
            |event: EventIdx| margins.get(event.position()).copied().map(Margin::is_present);
        let consistent = |pair: HapPair| {
            event_indices()
                .all(|event| present(event).is_some_and(|present| present == pair.carries(event)))
        };
        let called = earliest_max(accounts.iter().filter(|account| consistent(account.pair)))?;
        // The called account may be contaminated; its quality is its pair's
        // on its own.
        let called_pure = pure().find(|account| account.pair == called.pair)?;
        let pair_scores = HapPair::all(haplotypes)
            .map(|pair| {
                let best = accounts
                    .iter()
                    .filter(|account| account.pair == pair)
                    .map(|account| account.log10)
                    .fold(f64::NEG_INFINITY, f64::max);
                (pair, best)
            })
            .collect();
        Some(Self {
            best,
            called: called.pair,
            quality: quality_of(called_pure),
            margins,
            pair_scores,
        })
    }
}

/// Pairs over [`MAX_HAPLOTYPES`] haplotypes, each on its own and with each of
/// up to [`MAX_EVENTS`] events it lacks, and the reference, as an artefact at
/// every fraction.
const MAX_ACCOUNTS: usize = MAX_PAIRS * (1 + (MAX_EVENTS + 1) * ARTEFACT_FRACTIONS.len());

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
/// of an `artefact` haplotype mixed in at an [`ARTEFACT_FRACTIONS`] entry, its
/// prior already paid in `log10`.
struct Account {
    pair: HapPair,
    artefact: Option<Artefact>,
    log10: f64,
}

/// The likeliest account, the earliest of equals.
fn earliest_max<'a>(accounts: impl Iterator<Item = &'a Account>) -> Option<&'a Account> {
    accounts.fold(None, |best, account| match best {
        Some(best) if best.log10 >= account.log10 => Some(best),
        _ => Some(account),
    })
}

/// How much likelier a cluster's reads are with an event present than without
/// it: `log10` of the best account whose pair carries it over the best account
/// without it, a pair without the event or such a pair with the event's reads
/// as an artefact at one of [`ARTEFACT_FRACTIONS`], whatever its zygosity.
/// `-inf` when nothing carrying the event explains the reads, never `+inf` (an
/// artefact account always explains them) and never NaN, as the best account
/// always does.
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
    /// event's haplotype at one of [`ARTEFACT_FRACTIONS`] instead.
    pub best: HapPair,
    /// The genotype the records write: the likeliest pair carrying exactly
    /// the events called present, which is the best pair unless an event in
    /// it misses [`Margin::PRESENT`].
    called: HapPair,
    /// Phred of the chance the called pair is wrong, under a flat prior over
    /// the pairs on their own (artefact accounts aside); low whenever it is
    /// not the best pair.
    pub quality: Phred,
    /// By [`EventIdx`].
    margins: SmallVec<Margin, MAX_EVENTS>,
    /// Each pair's best account, pure or with an artefact, in `GL` order.
    pair_scores: SmallVec<(HapPair, f64), MAX_PAIRS>,
}

/// How many chromosomes of the called pair carry an event.
#[derive(Debug, Clone, Copy)]
pub struct Dosage {
    /// 1 or 2.
    pub copies: u8,
    /// Phred of the chance the other count is right, under a flat prior over
    /// the pairs consistent with the presence calls.
    pub gq: Phred,
}

impl Verdict {
    /// `event`'s margin; `None` for an event not scored.
    #[must_use]
    pub fn margin(&self, event: EventIdx) -> Option<Margin> {
        self.margins.get(event.position()).copied()
    }

    /// The genotype the records write: the likeliest pair over the events
    /// called present.
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

    /// `event`'s copies in the called pair, and their GQ: the posterior of
    /// the other copy count under a flat prior over the pairs consistent with
    /// the presence calls, each at its best account. So every record of the
    /// cluster writes one genotype, the called pair. `None` when the called
    /// pair lacks the event.
    ///
    /// No genotype prior: a het prior of 0.184 log10 (GIAB's het:hom ratio)
    /// is finer than the likelihood's tie margins are calibrated, so it does
    /// not compose with them; it cost both SNV and indel records on chr12.
    #[must_use]
    pub fn dosage(&self, event: EventIdx) -> Option<Dosage> {
        let present = |e: EventIdx| self.margin(e).is_some_and(Margin::is_present);
        self.dosage_of(event, self.called, |pair| {
            (0..self.margins.len())
                .all(|i| EventIdx::new(i).is_some_and(|e| present(e) == pair.carries(e)))
        })
    }

    /// [`Self::dosage`] over the pairs `admits` takes instead of those the
    /// margins call: the copies in the likeliest of them, `None` where it
    /// lacks the event. SNV records admit the pairs the forests' calls allow.
    #[must_use]
    pub fn dosage_given(
        &self,
        event: EventIdx,
        admits: impl Fn(HapPair) -> bool + Copy,
    ) -> Option<Dosage> {
        let called = self.pair_scores.iter().filter(|(pair, _)| admits(*pair)).fold(
            None,
            |best: Option<(HapPair, f64)>, &(pair, log10)| match best {
                Some((_, top)) if top >= log10 => best,
                _ => Some((pair, log10)),
            },
        )?;
        self.dosage_of(event, called.0, admits)
    }

    /// `event`'s copies in `called`, and their GQ over the pairs `admits`
    /// takes.
    fn dosage_of(
        &self,
        event: EventIdx,
        called: HapPair,
        admits: impl Fn(HapPair) -> bool + Copy,
    ) -> Option<Dosage> {
        let copies = called.copies_of(event);
        if copies == 0 {
            return None;
        }
        let scores = || self.pair_scores.iter().filter(move |(pair, _)| admits(*pair));
        let best = scores().map(|&(_, log10)| log10).fold(f64::NEG_INFINITY, f64::max);
        if !best.is_finite() {
            return None;
        }
        let (mut same, mut other) = (0.0, 0.0);
        for &(pair, log10) in scores() {
            let weight = 10f64.powf(log10 - best);
            if pair.copies_of(event) == copies {
                same += weight;
            } else {
                other += weight;
            }
        }
        let gq = phred_of_log10((other / (same + other)).log10());
        Some(Dosage { copies, gq })
    }
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

    pub(crate) const NO: f32 = f32::NEG_INFINITY;
    pub(crate) const REF_READ: [f32; 3] = [-1.0, -4.0, -4.0];
    pub(crate) const ALT1_READ: [f32; 3] = [-4.0, -1.0, -4.0];
    pub(crate) const ALT2_READ: [f32; 3] = [-4.0, -4.0, -1.0];

    /// `n` anonymous reads of each `log10` row, at MAPQ 60.
    pub(crate) fn rows(rows: &[([f32; 3], usize)]) -> Vec<ReadScores> {
        rows.iter()
            .flat_map(|&(log10, n)| std::iter::repeat_n(ReadScores::new(None, 60, row(log10)), n))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::{fixtures::*, *};
    use proptest::prelude::*;
    use seqair_types::Base;

    fn deletion(n: usize) -> IndelAllele {
        IndelAllele::Deletion(std::iter::repeat_n(Base::T, n).collect())
    }

    fn read(log10: [f32; 3]) -> ReadScores {
        ReadScores::new(None, 60, row(log10))
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

    /// A strong reference read at MAPQ 10 costs the homozygote of the alt at
    /// most its MAPQ once mixed; unmixed it costs the full three log10.
    #[test]
    fn a_low_mapq_read_costs_a_pair_at_most_its_mapq() {
        let cluster =
            Cluster::new([Event { anchor: 0, allele: deletion(1).into() }].into_iter().collect())
                .expect("a cluster");
        let read = ReadScores::new(None, 10, row([-1.0, -4.0, NO]));
        let hom = pair(1, 1);
        let unmixed = read.log10_unmixed(&cluster, hom);
        assert!((unmixed + 4.0).abs() < 1e-9);
        let mixed = read.log10_under(&cluster, hom);
        assert!((-2.0..-1.0).contains(&mixed), "{mixed}");
        let reference = read.log10_under(&cluster, pair(0, 0));
        assert!((reference + 1.0).abs() < 0.05, "{reference}");
        let unmapped = ReadScores::new(None, 0, row([-1.0, -4.0, NO]));
        assert!((unmapped.log10_under(&cluster, hom) + 1.0).abs() < 1e-9);
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
    fn dosage_is_the_best_carrying_pair_at_its_best_account() {
        let dosage =
            |reads: &[([f32; 3], usize)]| verdict(scores(1, reads)).dosage(ev(0)).expect("present");
        let stray = dosage(&[([0.0, -4.5, NO], 1), ([-4.5, 0.0, NO], 12)]);
        assert_eq!(stray.copies, 2);
        let het = dosage(&[([0.0, -4.5, NO], 6), ([-4.5, 0.0, NO], 6)]);
        assert_eq!(het.copies, 1);
        assert!(*het.gq > 30.0, "{:?}", het.gq);
        let hom = dosage(&[([-4.5, 0.0, NO], 3)]);
        let deeper = dosage(&[([-4.5, 0.0, NO], 12)]);
        assert_eq!((hom.copies, deeper.copies), (2, 2));
        assert!(*deeper.gq > *hom.gq, "{:?} {:?}", hom.gq, deeper.gq);
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
        let snv = |anchor, base| Event { anchor, allele: Allele::Base(base) };
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

        let mut read = ReadScores::new(None, 60, [-1.0, -2.0, -3.0, -4.0, NO, NO, NO, NO]);
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
        let at = |anchor, allele: IndelAllele| Event { anchor, allele: allele.into() };
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
        let snv = |anchor| Event { anchor, allele: Allele::Base(Base::A) };
        let del = Event { anchor: 10, allele: deletion(3).into() };
        let valid = |events: &[Event]| -> Vec<u8> {
            HapIdx::all_valid(events).iter().map(|h| h.0).collect()
        };
        for inside in 11..=13 {
            assert_eq!(valid(&[del.clone(), snv(inside)]), [0, 1, 2], "SNV at {inside}");
        }
        for outside in [9, 14, 30] {
            assert_eq!(valid(&[del.clone(), snv(outside)]), [0, 1, 2, 3], "SNV at {outside}");
        }
        assert_eq!(valid(&[del.clone(), snv(10)]), [0, 1, 2], "one anchor");
        let ins = Event {
            anchor: 12,
            allele: IndelAllele::Insertion([Base::G].into_iter().collect()).into(),
        };
        assert_eq!(valid(&[snv(12), ins.clone()]), [0, 1, 2], "one anchor");
        assert_eq!(valid(&[snv(13), ins]), [0, 1, 2, 3]);
        assert_eq!(valid(&[snv(5), snv(6)]), [0, 1, 2, 3]);
    }

    /// Two substitutions at one base are a multi-allelic site: never one
    /// chromosome, and at most two of them.
    #[test]
    fn substitutions_at_one_base_follow_the_anchor_limits() {
        let snv = |base| Event { anchor: 7, allele: Allele::Base(base) };
        let two: SmallVec<Event, MAX_EVENTS> = [snv(Base::A), snv(Base::C)].into_iter().collect();
        let cluster = Cluster::new(two).expect("two at one anchor");
        assert_eq!(cluster.haplotypes(), [HapIdx(0), HapIdx(1), HapIdx(2)]);
        assert!(!Cluster::admits(cluster.events(), &snv(Base::G)));
    }

    /// A third event at one anchor, or a fourth in all, makes no cluster.
    #[test]
    fn a_cluster_holds_two_events_per_anchor_and_three_in_all() {
        let at = |anchor, n| Event { anchor, allele: deletion(n).into() };
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
        let at = |anchor, allele: IndelAllele| Event { anchor, allele: allele.into() };
        let events = [at(10, deletion(1)), at(20, deletion(2))].into_iter().collect();
        let both = [-9.0, -5.0, -5.0, -1.0, NO, NO, NO, NO];
        let reads: Vec<ReadScores> =
            std::iter::repeat_n(ReadScores::new(None, 60, both), 8).collect();
        let cluster = Cluster::new(events).expect("a cluster");
        let s = HaplotypeScores::new(cluster, reads.into()).expect("a verdict");
        assert_eq!(s.verdict.called_pair(), pair(3, 3));
        assert_eq!(s.verdict.copies_of(ev(0)), 2);
        assert_eq!(s.verdict.copies_of(ev(1)), 2);
        assert!(s.verdict.margin(ev(0)).expect("scored").is_present());
        assert_eq!(s.support_over(&[ev(1)]).as_ref(), [0, 8]);
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
        assert!(v.quality.as_int() > 10, "{:?}", v.quality);
    }

    #[test]
    fn reads_of_one_alt_are_a_homozygote_and_of_two_a_compound_het() {
        assert_eq!(verdict(scores(1, &[(ALT1_READ, 8)])).called_pair(), pair(1, 1));
        let compound = verdict(scores(2, &[(ALT1_READ, 4), (ALT2_READ, 4)]));
        assert_eq!(compound.called_pair(), pair(1, 2));
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
        assert_eq!(v.quality.as_int(), 1);
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
    /// read carries the reference, a het when one does. Either way `GC` says
    /// the written pair is not the one the reads like best.
    #[test]
    fn an_uncalled_half_of_a_compound_leaves_the_other_as_the_reads_have_it() {
        let hom = verdict(scores(2, &[([-9.0, -1.0, -9.0], 2), ([-9.0, -1.5, -1.0], 4)]));
        assert_eq!(hom.best, pair(1, 2));
        assert!(hom.margin(ev(0)).expect("scored").is_present());
        let second = hom.margin(ev(1)).expect("scored");
        assert!(second > Margin(0.0) && !second.is_present(), "{second:?}");
        assert_eq!(hom.called_pair(), pair(1, 1));
        assert!(hom.quality.as_int() < 3, "{:?}", hom.quality);

        let het = verdict(scores(
            2,
            &[([-1.0, -9.0, -2.0], 1), ([-9.0, -1.0, -9.0], 2), ([-9.0, -1.5, -1.0], 2)],
        ));
        assert_eq!(het.best, pair(1, 2));
        assert!(!het.margin(ev(1)).expect("scored").is_present());
        assert_eq!(het.called_pair(), pair(0, 1));
        assert!(het.quality.as_int() < 3, "{:?}", het.quality);
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

    /// A read counts for its likeliest haplotype when that is informative;
    /// a read torn between two counts for neither, and not for the depth.
    #[test]
    fn support_counts_reads_assigned_to_a_haplotype() {
        let torn = ([-4.0, -1.0, -1.1], 3);
        let leaning = ([-4.0, -1.0, -1.3], 1);
        let s = scores(2, &[(REF_READ, 4), (ALT1_READ, 2), (ALT2_READ, 5), torn, leaning])
            .expect("a verdict");
        assert_eq!(s.support(HapIdx::REFERENCE), 4);
        assert_eq!(s.support(alt(0)), 3);
        assert_eq!(s.support(alt(1)), 5);
        assert_eq!(s.depth(), 12);

        // With one allele scored, the third haplotype does not compete.
        let one = scores(1, &[(REF_READ, 1), ([-4.0, -1.0, 0.0], 2)]).expect("a verdict");
        assert_eq!(one.support(alt(0)), 2);
        assert_eq!(one.depth(), 3);
    }

    #[test]
    fn log10_sum_is_stable_and_total() {
        assert!((log10_sum(-1.0, -1.0) - (-1.0 + std::f64::consts::LOG10_2)).abs() < 1e-12);
        assert!((log10_sum(-400.0, f64::NEG_INFINITY) - -400.0).abs() < 1e-9);
        assert_eq!(log10_sum(f64::NEG_INFINITY, f64::NEG_INFINITY), f64::NEG_INFINITY);
    }

    /// The best pair is called, and `GC` is its own confidence.
    #[test]
    fn a_called_best_pair_keeps_its_quality() {
        let v = verdict(scores(1, &[(REF_READ, 4), (ALT1_READ, 4)]));
        assert_eq!(v.called_pair(), v.best);
        assert!(v.quality.as_int() > 10, "{:?}", v.quality);
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

    proptest! {
        /// Presence is decided inside the best pair, and a margin above zero
        /// already implies it: the best pair carrying an allele beats every
        /// pair without it, so it is the best pair overall.
        #[test]
        fn a_present_allele_is_always_in_the_called_genotype(
            rows in prop::collection::vec(prop::array::uniform3(-30.0f32..0.0), 1..12),
            alts in 1usize..=2,
        ) {
            let v = verdict(HaplotypeScores::at_column(
                (1..=alts).map(deletion).collect(),
                rows.into_iter().map(read).collect(),
            ));
            for event in (0..alts).filter_map(EventIdx::new) {
                let margin = v.margin(event).expect("scored");
                if margin > Margin(0.0) {
                    prop_assert!(v.best.carries(event), "{v:?}");
                }
                prop_assert_eq!(margin.is_present(), v.called_pair().carries(event), "{:?}", v);
            }
        }
    }
}
