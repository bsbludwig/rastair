//! Which haplotype a read came from, and what that says at a column TAPS
//! conversion makes ambiguous.
//!
//! A reference C read as T on the top strand is either a methylated C or a C>T
//! variant; split by chromosome the two differ, because a variant piles its T
//! reads onto one haplotype and methylation need not. At a CpG the split is
//! that CpG's methylation per haplotype.
//!
//! A fragment's haplotype is read off the heterozygous sites, SNVs and the
//! pair-HMM's indels alike, from the same votes the solver phases them on: a
//! phased site places it in its block's
//! frame, an unphased one in a one-site set of its own. A column is never
//! scored on the votes it cast itself, so what it reports is what the
//! surrounding heterozygotes say.

use super::{
    FragmentId, GtAllele, PhaseObservation, PhaseSet, RecordPhases, observations_of,
    problem::{SiteVotes, het_alleles, indel_site, snv_site},
};
use crate::{
    call::variant_calling::GenotypeTag,
    metrics::{AltCall, PileupMetrics, methylation::CpgSide},
    utils::{ByHaplotype, ByStrand, Haplotype},
};
use rustc_hash::FxHashMap;
use seqair_types::{Base, Pos0, SmallVec, Strand};
use std::{
    cmp::{Ordering, Reverse},
    num::{NonZeroU8, NonZeroUsize},
};
use tracing::debug;

/// Whether one haplotype carries a CpG at a column, and so whether its
/// conversion rate is a methylation level at all.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum CpgPresence {
    /// It carries the pair's unconverted base.
    Present,
    /// It does not, so its reads convert by construction.
    Absent,
    /// Nothing at this column could say.
    #[default]
    Unresolved,
}

impl CpgPresence {
    /// The `HPX` flag, `None` where the column could not say.
    #[must_use]
    pub const fn flag(self) -> Option<i32> {
        match self {
            Self::Present => Some(1),
            Self::Absent => Some(0),
            Self::Unresolved => None,
        }
    }

    /// A CpG is on a haplotype only where both of its bases are.
    const fn and(self, other: Self) -> Self {
        match (self, other) {
            (Self::Absent, _) | (_, Self::Absent) => Self::Absent,
            (Self::Present, Self::Present) => Self::Present,
            _ => Self::Unresolved,
        }
    }
}

/// Fragments of one phase set, each counted once by the majority of its CpG
/// reads, so that a molecule crossing several CpGs is one observation.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct FragmentCounts {
    pub unmethylated: u32,
    pub methylated: u32,
}

impl FragmentCounts {
    fn count(&mut self, verdict: Methylation) {
        let slot = match verdict {
            Methylation::Methylated => &mut self.methylated,
            Methylation::Unmethylated => &mut self.unmethylated,
        };
        *slot = slot.saturating_add(1);
    }
}

/// A column's reads on its TAPS-confounded strand, counted per haplotype.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HaplotypeSplit {
    /// The set the haplotypes belong to, written out as `HPS`: a `PS` this file
    /// emits, or an unphased heterozygote's own position, in which case
    /// haplotype 1 carries the first allele of that site's `GT`.
    pub set: PhaseSet,
    /// Reads showing the unconverted base (`C` at a C/T column, `G` at G/A).
    pub unconverted: ByHaplotype<u32>,
    /// Reads showing the converted base (`T` at a C/T column, `A` at G/A).
    pub converted: ByHaplotype<u32>,
    /// Written out as `HPX`.
    pub cpg: ByHaplotype<CpgPresence>,
    /// The whole set's fragment counts, the same on every record of the set.
    /// Written out as `HPF`.
    pub fragments: ByHaplotype<FragmentCounts>,
}

impl HaplotypeSplit {
    /// `HPC`: unconverted then converted on haplotype 1, then on haplotype 2.
    #[must_use]
    pub const fn counts(&self) -> [u32; 4] {
        [
            self.unconverted.first,
            self.converted.first,
            self.unconverted.second,
            self.converted.second,
        ]
    }

    /// `HPF`, in `HPC`'s order: unmethylated then methylated fragments on
    /// haplotype 1, then on haplotype 2.
    #[must_use]
    pub const fn fragment_counts(&self) -> [u32; 4] {
        [
            self.fragments.first.unmethylated,
            self.fragments.first.methylated,
            self.fragments.second.unmethylated,
            self.fragments.second.methylated,
        ]
    }

    fn depth(&self) -> u32 {
        self.counts().into_iter().fold(0, u32::saturating_add)
    }
}

/// A heterozygous site's right to place fragments on a haplotype.
#[derive(Debug, Clone, Copy)]
struct Anchor {
    set: PhaseSet,
    /// The allele haplotype 1 carries.
    first: GtAllele,
    kind: AnchorKind,
}

impl Anchor {
    /// A phased site anchors its block; an unphased one a set of its own, filed
    /// under its column's position so it cannot collide with a block's.
    ///
    /// Such a set is named only by position, so a column mints at most one, and
    /// none once any of its records is phased: `minted` carries that across
    /// the column's sites.
    fn of(sites: &SiteVotes, phases: RecordPhases, minted: &mut bool) -> Option<Self> {
        let [lower, higher] = sites.site.alleles;
        match phases.get(sites.site.kind) {
            Some(call) if call.first == lower || call.first == higher => {
                Some(Self { set: call.set, first: call.first, kind: AnchorKind::Block })
            }
            // A call naming neither allele would silently invert half a block.
            Some(_) => None,
            None if phases.snv.is_none() && phases.indel.is_none() && !*minted => {
                *minted = true;
                Some(Self {
                    set: PhaseSet::at(sites.site.pos)?,
                    first: lower,
                    kind: AnchorKind::Singleton { confounds: sites.confounded_on },
                })
            }
            None => None,
        }
    }

    fn haplotype_of(&self, allele: GtAllele) -> Haplotype {
        if allele == self.first { Haplotype::First } else { Haplotype::Second }
    }
}

/// Which kind of anchor a placement came from; a singleton also carries the
/// strand its own alleles are confounded on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AnchorKind {
    Block,
    Singleton { confounds: Option<Strand> },
}

/// Every fragment's votes per phase set, and each column's own share of them.
#[derive(Debug, Default)]
struct HaplotypeAssignment {
    votes: FxHashMap<FragmentId, Votes>,
    /// Indexed like the segment's columns.
    contributed: Vec<Option<ColumnContribution>>,
}

impl HaplotypeAssignment {
    fn contribution(&self, column: usize) -> Option<&ColumnContribution> {
        self.contributed.get(column).and_then(Option::as_ref)
    }

    /// At most one placement per set the fragment reached, with `exclude`'s
    /// own votes left out.
    fn placements(
        &self,
        fragment: FragmentId,
        exclude: Option<&ColumnContribution>,
    ) -> SmallVec<Placement, 2> {
        let Some(votes) = self.votes.get(&fragment) else {
            return SmallVec::new();
        };
        votes.placements(|set| exclude.map_or(0, |column| column.own_vote(fragment, set)))
    }
}

/// The votes one column cast, so they can be taken back out when that column
/// is the one being scored. A column holding an SNV and an indel site can
/// vote in two sets.
#[derive(Debug, Default)]
struct ColumnContribution {
    votes: FxHashMap<(FragmentId, PhaseSet), i32>,
}

impl ColumnContribution {
    fn add(&mut self, fragment: FragmentId, set: PhaseSet, signed: i32) {
        let slot = self.votes.entry((fragment, set)).or_default();
        *slot = slot.saturating_add(signed);
    }

    fn own_vote(&self, fragment: FragmentId, set: PhaseSet) -> i32 {
        self.votes.get(&(fragment, set)).copied().unwrap_or(0)
    }
}

/// A fragment's haplotype, and the set it means anything in.
#[derive(Debug, Clone, Copy)]
struct Placement {
    set: PhaseSet,
    haplotype: Haplotype,
    kind: AnchorKind,
}

impl Placement {
    /// Did this placement's own anchor decide where a read on this strand
    /// landed?
    ///
    /// A singleton confounded on a strand places only the templates carrying
    /// one of its alleles there, so its split on that strand is one-sided by
    /// construction, and selected on the anchor's own methylation. A block
    /// pools several hets and keeps every placement.
    fn one_sided_at(self, strand: Strand) -> bool {
        match self.kind {
            AnchorKind::Block => false,
            AnchorKind::Singleton { confounds } => confounds == Some(strand),
        }
    }
}

/// Votes for one fragment, one entry per set its reads reached.
#[derive(Debug, Default)]
struct Votes(SmallVec<SetVotes, 2>);

#[derive(Debug, Clone, Copy)]
struct SetVotes {
    set: PhaseSet,
    kind: AnchorKind,
    /// Positive leans towards haplotype 1, negative towards haplotype 2.
    score: i32,
}

impl Votes {
    fn add(&mut self, set: PhaseSet, kind: AnchorKind, signed: i32) {
        match self.0.iter_mut().find(|votes| votes.set == set) {
            Some(votes) => votes.score = votes.score.saturating_add(signed),
            None => self.0.push(SetVotes { set, kind, score: signed }),
        }
    }

    /// A tie is no evidence, so a fragment whose reads contradict each other in
    /// a set has no haplotype there; `own` is taken out of each set first.
    fn placements(&self, own: impl Fn(PhaseSet) -> i32) -> SmallVec<Placement, 2> {
        self.0
            .iter()
            .filter_map(|votes| {
                let score = votes.score.saturating_sub(own(votes.set));
                (score != 0).then_some(Placement {
                    set: votes.set,
                    haplotype: if score > 0 { Haplotype::First } else { Haplotype::Second },
                    kind: votes.kind,
                })
            })
            .collect()
    }
}

/// The sites of a column that may anchor a split: its SNV record, then its
/// indel record, whose votes are the pair-HMM's read assignments. The SNV
/// record comes first, so it is the one an unphased column's set is named for.
fn anchor_sites(pileup: &PileupMetrics, index: usize) -> impl Iterator<Item = SiteVotes> {
    snv_site(pileup, index).into_iter().chain(indel_site(pileup, index))
}

/// Read each fragment's haplotype off the heterozygous sites the segment
/// emits.
fn assign_haplotypes(
    pileups: &[PileupMetrics],
    phases: &[RecordPhases],
    is_core: &impl Fn(u64) -> bool,
) -> HaplotypeAssignment {
    let mut votes: FxHashMap<FragmentId, Votes> = FxHashMap::default();
    let mut contributed = Vec::with_capacity(pileups.len());

    for (index, pileup) in pileups.iter().enumerate() {
        let contribution = is_core(u64::from(pileup.pos))
            .then(|| {
                let phases = phases.get(index).copied().unwrap_or_default();
                column_votes(pileup, index, phases, &mut votes)
            })
            .filter(|own| !own.votes.is_empty());
        contributed.push(contribution);
    }

    HaplotypeAssignment { votes, contributed }
}

/// Cast one column's votes, returning them as its contribution.
fn column_votes(
    pileup: &PileupMetrics,
    index: usize,
    phases: RecordPhases,
    votes: &mut FxHashMap<FragmentId, Votes>,
) -> ColumnContribution {
    let mut own = ColumnContribution::default();
    let mut minted = false;
    for sites in anchor_sites(pileup, index) {
        let Some(anchor) = Anchor::of(&sites, phases, &mut minted) else {
            continue;
        };
        for vote in &sites.votes {
            let weight = i32::from(vote.cost);
            let signed = match anchor.haplotype_of(sites.site.allele(vote.allele)) {
                Haplotype::First => weight,
                Haplotype::Second => -weight,
            };
            votes.entry(vote.fragment).or_default().add(anchor.set, anchor.kind, signed);
            own.add(vote.fragment, anchor.set, signed);
        }
    }
    own
}

/// Split every core column the assignment reaches by haplotype, then settle
/// the CpG partners and count fragments per set.
///
/// `phases` holds the segment's phase calls, indexed like `pileups`, and
/// receives the splits. Only sites the segment emits anchor, so a set is a
/// `PS` this file writes or a het it writes unphased.
pub(super) fn annotate_haplotype_split(
    pileups: &[PileupMetrics],
    phases: &mut [RecordPhases],
    is_core: impl Fn(u64) -> bool,
) {
    let assignment = assign_haplotypes(pileups, phases, &is_core);

    let mut split_columns = 0usize;
    for (index, (pileup, phase)) in pileups.iter().zip(phases.iter_mut()).enumerate() {
        if !is_core(u64::from(pileup.pos)) {
            continue;
        }
        phase.split = best_tally(pileup, &assignment, index)
            .and_then(|tally| confounded_split(pileup, &tally));
        split_columns += usize::from(phase.split.is_some());
    }

    propagate_partner_presence(pileups, phases);
    // Last: a fragment only counts where both haplotypes carry the CpG, which
    // the partner pass settles at the `G` side of a heterozygous one.
    count_fragments(pileups, phases, &assignment);

    debug!(
        fragments = assignment.votes.len(),
        columns = split_columns,
        "Split columns by haplotype"
    );
}

/// Carry each CpG side's presence across to the other base of the CpG.
///
/// A `C>T` destroys the CpG on one chromosome, but the `G` a base later is a
/// `G` on both, so left to itself that column would report a cytosine the
/// variant took away. A haplotype carries the CpG only where both sides say
/// so; a side that cannot say, or says it in another set's frame, leaves the
/// partner unresolved.
///
/// `pileups` must be sorted by position, as a segment's columns are.
fn propagate_partner_presence(pileups: &[PileupMetrics], phases: &mut [RecordPhases]) {
    debug_assert!(pileups.is_sorted_by_key(|pileup| pileup.pos));
    // Gathered first, so a partner that is itself a source never reads what
    // its neighbour just wrote.
    let partners: Vec<(usize, CpgSide, PartnerPresence)> = pileups
        .iter()
        .zip(phases.iter())
        .filter_map(|(pileup, phase)| {
            if !pileup.tags.in_cpg() {
                return None;
            }
            let side = CpgSide::of_base(pileup.reference_base)?;
            let presence = PartnerPresence::of(pileup, phase.split, side)?;
            let partner = side.partner_position(Pos0::new(pileup.pos)?)?;
            // The trims leave gaps, so the neighbour is searched for.
            let index =
                pileups.binary_search_by_key(&partner.as_u32(), |pileup| pileup.pos).ok()?;
            Some((index, side, presence))
        })
        .collect();

    for (index, side, presence) in partners {
        let (Some(target), Some(phase)) = (pileups.get(index), phases.get_mut(index)) else {
            continue;
        };
        if !target.tags.in_cpg()
            || CpgSide::of_base(target.reference_base).is_none_or(|other| other == side)
        {
            continue;
        }
        if let Some(split) = phase.split.as_mut() {
            let incoming = presence.in_set(split.set);
            split.cpg = ByHaplotype::new(
                split.cpg.first.and(incoming.first),
                split.cpg.second.and(incoming.second),
            );
        }
    }
}

/// What one side of a CpG that is not homozygous reference says about the
/// haplotypes carrying the whole CpG.
#[derive(Debug, Clone, Copy)]
enum PartnerPresence {
    /// The same on both haplotypes, so it holds in every set.
    Everywhere(ByHaplotype<CpgPresence>),
    /// A heterozygote's resolved label, meaningful only in its own set.
    InSet(PhaseSet, ByHaplotype<CpgPresence>),
}

impl PartnerPresence {
    fn of(pileup: &PileupMetrics, split: Option<HaplotypeSplit>, side: CpgSide) -> Option<Self> {
        match pileup.pos_metrics.extended.genotype?.genotype {
            GenotypeTag::HomRef => None,
            GenotypeTag::HomAlt(n) => Some(Self::Everywhere(homozygous_alt(pileup, n, side))),
            GenotypeTag::RefHet(_) | GenotypeTag::AltHet(..) => Some(match split {
                Some(split) => Self::InSet(split.set, split.cpg),
                None => Self::Everywhere(both(CpgPresence::Unresolved)),
            }),
        }
    }

    /// The presence in `set`'s frame.
    fn in_set(self, set: PhaseSet) -> ByHaplotype<CpgPresence> {
        match self {
            Self::Everywhere(presence) => presence,
            Self::InSet(own, presence) if own == set => presence,
            Self::InSet(..) => both(CpgPresence::Unresolved),
        }
    }
}

/// What a read, or a fragment's majority, says about methylation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Methylation {
    Methylated,
    Unmethylated,
}

/// One fragment's CpG reads inside one set.
#[derive(Debug, Clone, Copy)]
struct FragmentVote {
    /// `None` once two columns placed the fragment on different haplotypes.
    haplotype: Option<Haplotype>,
    converted: u32,
    unconverted: u32,
}

impl FragmentVote {
    const fn new(haplotype: Haplotype) -> Self {
        Self { haplotype: Some(haplotype), converted: 0, unconverted: 0 }
    }

    fn observe(&mut self, haplotype: Haplotype, read: Methylation) {
        if self.haplotype != Some(haplotype) {
            self.haplotype = None;
        }
        let slot = match read {
            Methylation::Methylated => &mut self.converted,
            Methylation::Unmethylated => &mut self.unconverted,
        };
        *slot = slot.saturating_add(1);
    }

    /// The majority, where there is one; a tie carries no direction.
    fn verdict(self) -> Option<(Haplotype, Methylation)> {
        let verdict = match self.converted.cmp(&self.unconverted) {
            Ordering::Greater => Methylation::Methylated,
            Ordering::Less => Methylation::Unmethylated,
            Ordering::Equal => return None,
        };
        Some((self.haplotype?, verdict))
    }
}

/// Count each set's CpG evidence by fragment and write it onto every column of
/// the set.
///
/// The reads are the ones `HPC` counts, narrowed to CpG columns both
/// haplotypes carry.
fn count_fragments(
    pileups: &[PileupMetrics],
    phases: &mut [RecordPhases],
    assignment: &HaplotypeAssignment,
) {
    let mut votes: FxHashMap<(PhaseSet, FragmentId), FragmentVote> = FxHashMap::default();

    for (index, (pileup, phase)) in pileups.iter().zip(phases.iter()).enumerate() {
        let (Some(split), Some(observations)) = (phase.split, observations_of(pileup)) else {
            continue;
        };
        // A haplotype without the CpG converts by construction, and counting
        // only the other one would test composition rather than methylation.
        if !pileup.tags.in_cpg() || split.cpg != both(CpgPresence::Present) {
            continue;
        }
        let Some(side) = confounded_side(pileup) else {
            continue;
        };
        let own = assignment.contribution(index);

        for observation in observations {
            let Some(strand) = usable_strand(observation) else {
                continue;
            };
            if strand != side.strand() {
                continue;
            }
            let read = if observation.base == side.mod_base() {
                Methylation::Methylated
            } else if observation.base == side.unmod_base() {
                Methylation::Unmethylated
            } else {
                continue;
            };
            let Some(placement) = assignment
                .placements(observation.fragment, own)
                .into_iter()
                .find(|placement| placement.set == split.set)
            else {
                continue;
            };
            if placement.one_sided_at(strand) {
                continue;
            }
            votes
                .entry((split.set, observation.fragment))
                .or_insert_with(|| FragmentVote::new(placement.haplotype))
                .observe(placement.haplotype, read);
        }
    }

    let mut counted: FxHashMap<PhaseSet, ByHaplotype<FragmentCounts>> = FxHashMap::default();
    for ((set, _), vote) in votes {
        if let Some((haplotype, verdict)) = vote.verdict() {
            counted.entry(set).or_default().get_mut(haplotype).count(verdict);
        }
    }

    for split in phases.iter_mut().filter_map(|phase| phase.split.as_mut()) {
        split.fragments = counted.get(&split.set).copied().unwrap_or_default();
    }
}

/// Every read one set placed at a column, by haplotype, base and strand.
#[derive(Debug, Clone, Copy)]
struct HaplotypeTally {
    set: PhaseSet,
    kind: AnchorKind,
    /// Per haplotype, indexed by [`Base::known_index`].
    placed: ByHaplotype<[ByStrand<u32>; Base::KNOWN.len()]>,
}

impl HaplotypeTally {
    fn new(set: PhaseSet, kind: AnchorKind) -> Self {
        Self { set, kind, placed: ByHaplotype::default() }
    }

    fn add(&mut self, haplotype: Haplotype, base: Base, strand: Strand) {
        let Some(slot) = base
            .known_index()
            .and_then(|index| self.placed.get_mut(haplotype).get_mut(index))
            .and_then(|counts| counts.get_mut(strand))
        else {
            return;
        };
        *slot = slot.saturating_add(1);
    }

    fn count(&self, haplotype: Haplotype, base: Base, strand: Strand) -> u32 {
        base.known_index()
            .and_then(|index| self.placed.get(haplotype).get(index))
            .and_then(|counts| counts.get(strand))
            .copied()
            .unwrap_or(0)
    }

    /// Placed reads over every base and both strands.
    fn depth(&self) -> u32 {
        [&self.placed.first, &self.placed.second]
            .into_iter()
            .flatten()
            .map(ByStrand::total)
            .fold(0, u32::saturating_add)
    }

    /// The confounded-strand counts `HPC` reports, without the presence.
    fn confounded(&self, side: CpgSide) -> HaplotypeSplit {
        let on = |haplotype, base| self.count(haplotype, base, side.strand());
        HaplotypeSplit {
            set: self.set,
            unconverted: ByHaplotype::new(
                on(Haplotype::First, side.unmod_base()),
                on(Haplotype::Second, side.unmod_base()),
            ),
            converted: ByHaplotype::new(
                on(Haplotype::First, side.mod_base()),
                on(Haplotype::Second, side.mod_base()),
            ),
            cpg: ByHaplotype::default(),
            fragments: ByHaplotype::default(),
        }
    }

    /// One haplotype's reads on the strand that tells the pair's bases apart.
    fn label(&self, haplotype: Haplotype, side: CpgSide) -> LabelCounts {
        let on = |base| self.count(haplotype, base, side.unconfounded_strand());
        LabelCounts { converted: on(side.mod_base()), unconverted: on(side.unmod_base()) }
    }
}

/// One haplotype's share of the label strand at a column heterozygous for the
/// confounded pair; there the two counts are alleles, not conversion.
#[derive(Debug, Clone, Copy)]
struct LabelCounts {
    converted: u32,
    unconverted: u32,
}

impl LabelCounts {
    /// Fewer converted reads than this may be sequencing errors.
    const MIN_READS: u32 = 2;

    const fn total(self) -> u32 {
        self.converted.saturating_add(self.unconverted)
    }

    /// At least four reads in five carry the converted base.
    const fn carries_converted(self) -> bool {
        self.converted >= Self::MIN_READS
            && self.converted.saturating_mul(5) >= self.total().saturating_mul(4)
    }

    /// At most one read in five does; no reads at all counts, since a
    /// heterozygote's other allele is the one left over.
    const fn carries_unconverted(self) -> bool {
        self.converted.saturating_mul(5) <= self.total()
    }
}

/// The pair a column is TAPS-confounded on, if it has anything to say in it.
///
/// A CpG always does, converted reads or not: an unmethylated CpG is the null.
/// A de-novo CpG's partner never does: its reads are kept only where they
/// look heterozygous, so its split would be a selected sample.
/// Elsewhere the pair's other base has to be a called alt.
fn confounded_side(pileup: &PileupMetrics) -> Option<CpgSide> {
    let side = CpgSide::of_base(pileup.reference_base)?;
    let tags = &pileup.tags;
    if tags.cpg || tags.denovo_cpg {
        return Some(side);
    }
    if tags.denovo_cpg_partner {
        return None;
    }
    let partner = side.other_base(pileup.reference_base);
    pileup
        .alts
        .iter()
        .any(|alt| alt.base == partner && alt.call == AltCall::RealVariant)
        .then_some(side)
}

/// The strand an observation can be counted on; a read with no strand or no
/// known base cannot be.
fn usable_strand(observation: &PhaseObservation) -> Option<Strand> {
    observation.base.known_index()?;
    observation.strand.ok()
}

/// Tally every set that placed a read at this column, on the evidence of every
/// other column.
fn column_tallies(
    observations: &[PhaseObservation],
    assignment: &HaplotypeAssignment,
    own: Option<&ColumnContribution>,
) -> SmallVec<HaplotypeTally, 2> {
    let mut tallies: SmallVec<HaplotypeTally, 2> = SmallVec::new();
    for observation in observations {
        let Some(strand) = usable_strand(observation) else {
            continue;
        };
        for placement in assignment.placements(observation.fragment, own) {
            if placement.one_sided_at(strand) {
                continue;
            }
            let tally = match tallies.iter().position(|tally| tally.set == placement.set) {
                Some(index) => tallies.get_mut(index),
                None => {
                    tallies.push(HaplotypeTally::new(placement.set, placement.kind));
                    tallies.last_mut()
                }
            };
            if let Some(tally) = tally {
                tally.add(placement.haplotype, observation.base, strand);
            }
        }
    }
    tallies
}

/// The one set a column reports in.
///
/// A set that can fill `HPC` first, then the most placed reads, then a block
/// over a singleton, then the lower set id.
fn best_tally(
    pileup: &PileupMetrics,
    assignment: &HaplotypeAssignment,
    column: usize,
) -> Option<HaplotypeTally> {
    let observations = observations_of(pileup)?;
    let confounding = confounded_side(pileup);
    let own = assignment.contribution(column);
    column_tallies(observations, assignment, own).into_iter().max_by_key(|tally| {
        let answers_confounded = confounding.is_some_and(|side| tally.confounded(side).depth() > 0);
        let is_block = tally.kind == AnchorKind::Block;
        (answers_confounded, tally.depth(), is_block, Reverse(tally.set))
    })
}

/// What `HPC` and `HPX` report, or `None` where the column has no ambiguity to
/// resolve or no placed read on the confounded strand.
fn confounded_split(pileup: &PileupMetrics, tally: &HaplotypeTally) -> Option<HaplotypeSplit> {
    let side = confounded_side(pileup)?;
    let counts = tally.confounded(side);
    (counts.depth() > 0)
        .then(|| HaplotypeSplit { cpg: cpg_presence(pileup, tally, side), ..counts })
}

/// Which haplotype has a cytosine here to methylate: the one carrying the
/// pair's unconverted base, from the genotype for a homozygote and from the
/// label strand for a heterozygote.
fn cpg_presence(
    pileup: &PileupMetrics,
    tally: &HaplotypeTally,
    side: CpgSide,
) -> ByHaplotype<CpgPresence> {
    if !pileup.tags.in_cpg() {
        return both(CpgPresence::Absent);
    }
    let Some(genotype) = pileup.pos_metrics.extended.genotype else {
        return both(CpgPresence::Unresolved);
    };
    match genotype.genotype {
        GenotypeTag::HomRef => carrying(side, pileup.reference_base),
        GenotypeTag::HomAlt(n) => homozygous_alt(pileup, n, side),
        GenotypeTag::RefHet(_) | GenotypeTag::AltHet(..) => label_haplotypes(pileup, tally, side),
    }
}

/// Presence on both haplotypes of a homozygote for `base`.
fn carrying(side: CpgSide, base: Base) -> ByHaplotype<CpgPresence> {
    both(if base == side.unmod_base() { CpgPresence::Present } else { CpgPresence::Absent })
}

/// Presence on both haplotypes of a homozygote for its `n`-th alt.
fn homozygous_alt(pileup: &PileupMetrics, n: NonZeroU8, side: CpgSide) -> ByHaplotype<CpgPresence> {
    pileup
        .alts
        .get(NonZeroUsize::from(n).get().saturating_sub(1))
        .map_or_else(|| both(CpgPresence::Unresolved), |alt| carrying(side, alt.base))
}

/// Label the haplotypes of a heterozygote whose alleles are the confounded
/// pair itself, from the strand that reads them apart.
///
/// Any other heterozygote stays unresolved, as does a label strand too thin
/// or contradictory to decide: a mislabel would turn a real methylation level
/// into a missing one and a missing one into a fabricated 1.0.
fn label_haplotypes(
    pileup: &PileupMetrics,
    tally: &HaplotypeTally,
    side: CpgSide,
) -> ByHaplotype<CpgPresence> {
    let unresolved = both(CpgPresence::Unresolved);
    let Some(alleles) = het_alleles(pileup) else {
        return unresolved;
    };
    if CpgSide::of_alleles(alleles.map(|allele| allele.base)) != Some(side) {
        return unresolved;
    }

    let (first, second) =
        (tally.label(Haplotype::First, side), tally.label(Haplotype::Second, side));
    let (carrier, other, labelled) = match first.converted.cmp(&second.converted) {
        Ordering::Greater => {
            (first, second, ByHaplotype::new(CpgPresence::Absent, CpgPresence::Present))
        }
        Ordering::Less => {
            (second, first, ByHaplotype::new(CpgPresence::Present, CpgPresence::Absent))
        }
        Ordering::Equal => return unresolved,
    };
    if carrier.carries_converted() && other.carries_unconverted() { labelled } else { unresolved }
}

const fn both(presence: CpgPresence) -> ByHaplotype<CpgPresence> {
    ByHaplotype::new(presence, presence)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        call::phasing::{
            PhaseCall,
            fixtures::{
                ALT_1, ALT_2, call, column, confirming, graded, het, insertion, obs, with_indels,
            },
        },
        vcf::InCpG,
    };
    use seqair_types::Base::*;

    /// A column and the phase call its SNV record carries.
    type Column = (PileupMetrics, Option<PhaseCall>);

    fn phased(column: PileupMetrics, set: u32, first: u8) -> Column {
        (column, Some(call(set, first)))
    }

    const fn free(column: PileupMetrics) -> Column {
        (column, None)
    }

    /// Mark a column as a reference CpG, which `RecordTags::default` does not,
    /// on both the tag the split reads and the context the phasing gate reads.
    fn cpg(mut column: PileupMetrics) -> PileupMetrics {
        column.tags.cpg = true;
        column.pos_metrics.cpg = match column.reference_base {
            C => InCpG::C,
            G => InCpG::G,
            _ => InCpG::No,
        };
        column
    }

    /// Mark a column as a CpG a variant creates.
    fn denovo(mut column: PileupMetrics) -> PileupMetrics {
        column.tags.denovo_cpg = true;
        column
    }

    fn split_input(columns: &[Column]) -> (Vec<PileupMetrics>, Vec<RecordPhases>) {
        columns
            .iter()
            .cloned()
            .map(|(column, snv)| (column, RecordPhases { snv, ..RecordPhases::default() }))
            .unzip()
    }

    fn assign(columns: &[Column]) -> HaplotypeAssignment {
        let (pileups, phases) = split_input(columns);
        assign_haplotypes(&pileups, &phases, &|_| true)
    }

    /// A segment after the split.
    struct Annotated {
        phases: Vec<RecordPhases>,
    }

    impl Annotated {
        fn split(&self, index: usize) -> Option<HaplotypeSplit> {
            self.phases.get(index)?.split
        }

        /// The presence pair, in the order `HPX` writes it.
        fn presence(&self, index: usize) -> Option<[CpgPresence; 2]> {
            let cpg = self.split(index)?.cpg;
            Some([cpg.first, cpg.second])
        }

        fn set(&self, index: usize) -> Option<i32> {
            self.split(index).map(|split| split.set.as_i32())
        }

        fn split_columns(&self) -> usize {
            self.phases.iter().filter(|phase| phase.split.is_some()).count()
        }
    }

    fn annotate(columns: &[Column]) -> Annotated {
        let (pileups, mut phases) = split_input(columns);
        annotate_haplotype_split(&pileups, &mut phases, |_| true);
        Annotated { phases }
    }

    /// The same columns read off the other strand: every base complemented
    /// and OT swapped with OB, which turns a C/T fixture into its G/A mirror.
    fn mirror(columns: &[Column]) -> Vec<Column> {
        let strand = |strand| match strand {
            Strand::OT => Strand::OB,
            Strand::OB => Strand::OT,
            Strand::Unknown => Strand::Unknown,
        };
        columns
            .iter()
            .cloned()
            .map(|(mut column, call)| {
                column.reference_base = column.reference_base.inverse();
                column.pos_metrics.cpg = match column.pos_metrics.cpg {
                    InCpG::C => InCpG::G,
                    InCpG::G => InCpG::C,
                    InCpG::No => InCpG::No,
                };
                for alt in &mut column.alts {
                    alt.base = alt.base.inverse();
                }
                if let Some(ColumnPhase::Observed(observations)) = column.phase.as_deref_mut() {
                    for observation in observations.iter_mut() {
                        observation.base = observation.base.inverse();
                        observation.strand = strand(observation.strand);
                    }
                }
                (column, call)
            })
            .collect()
    }

    use crate::call::phasing::ColumnPhase;

    /// A reference CpG the caller only saw methylation evidence at.
    fn methylation_candidate(pos: u32, observations: &[PhaseObservation]) -> Column {
        free(cpg(column(
            pos,
            C,
            &[(T, AltCall::MethylationEvidenceOnly { for_base: C })],
            Some(GenotypeTag::HomRef),
            observations,
        )))
    }

    /// The bottom-strand mirror of [`methylation_candidate`].
    fn bottom_strand_candidate(pos: u32, observations: &[PhaseObservation]) -> Column {
        free(cpg(column(
            pos,
            G,
            &[(A, AltCall::MethylationEvidenceOnly { for_base: G })],
            Some(GenotypeTag::HomRef),
            observations,
        )))
    }

    /// Every set the fragment was placed in, ascending.
    fn all_placed(assignment: &HaplotypeAssignment, fragment: u64) -> Vec<(i32, Haplotype)> {
        let Some(id) = FragmentId::new(fragment) else {
            return Vec::new();
        };
        let mut placements: Vec<_> = assignment
            .placements(id, None)
            .iter()
            .map(|placement| (placement.set.as_i32(), placement.haplotype))
            .collect();
        placements.sort_by_key(|&(set, _)| set);
        placements
    }

    /// The fragment's placement where it has exactly one.
    fn placed(assignment: &HaplotypeAssignment, fragment: u64) -> Option<(i32, Haplotype)> {
        match all_placed(assignment, fragment).as_slice() {
            [placement] => Some(*placement),
            _ => None,
        }
    }

    /// An A/C heterozygote at 1000: `hap1` fragments carry the reference `A`,
    /// `hap2` the `C`. Neither allele is one conversion produces, so it places
    /// every template of both strands.
    fn labelled_anchor(hap1: &[(u64, Strand)], hap2: &[(u64, Strand)]) -> Column {
        let observations: Vec<PhaseObservation> = hap1
            .iter()
            .map(|&(fragment, strand)| obs(fragment, A, strand))
            .chain(hap2.iter().map(|&(fragment, strand)| obs(fragment, C, strand)))
            .collect();
        free(het(1000, A, C, &observations))
    }

    /// Four top-strand fragments, 1 and 2 on haplotype 1 (`A`), 3 and 4 on
    /// haplotype 2 (`G`).
    fn anchor_column() -> PileupMetrics {
        het(
            1000,
            A,
            G,
            &[
                obs(1, A, Strand::OT),
                obs(2, A, Strand::OT),
                obs(3, G, Strand::OT),
                obs(4, G, Strand::OT),
            ],
        )
    }

    /// [`anchor_column`] unphased, in a set of its own at 1001.
    fn singleton_anchor() -> Column {
        free(anchor_column())
    }

    /// [`anchor_column`], phased in the block at 1001.
    fn block_anchor() -> Column {
        phased(anchor_column(), 1001, 0)
    }

    /// A fragment takes the haplotype of the allele it carries, and the phase
    /// call says which allele is haplotype 1.
    #[test]
    fn a_fragment_takes_the_haplotype_of_the_allele_it_carries() {
        for (first, carriers) in
            [(0, [Haplotype::First, Haplotype::Second]), (1, [Haplotype::Second, Haplotype::First])]
        {
            let columns = [phased(
                het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OB)]),
                1001,
                first,
            )];

            let assignment = assign(&columns);

            let [reference, alt] = carriers;
            assert_eq!(placed(&assignment, 1), Some((1001, reference)), "first allele {first}");
            assert_eq!(placed(&assignment, 2), Some((1001, alt)), "first allele {first}");
        }
    }

    /// At a C/T site in a CpG the top-strand `T` may be a methylated `C`.
    #[test]
    fn a_taps_confounded_observation_places_nothing() {
        let columns = [phased(
            cpg(het(
                1000,
                C,
                T,
                &[
                    obs(1, T, Strand::OT),
                    obs(2, T, Strand::OB),
                    obs(3, C, Strand::OT),
                    obs(9, T, Strand::OB),
                ],
            )),
            1001,
            0,
        )];

        let assignment = assign(&columns);

        assert_eq!(placed(&assignment, 1), None, "an OT T at a C/T CpG site is ambiguous");
        assert_eq!(placed(&assignment, 2), Some((1001, Haplotype::Second)));
        assert_eq!(placed(&assignment, 3), Some((1001, Haplotype::First)));
    }

    /// A het the solver could not link still defines two local haplotypes, in
    /// a set of its own whose haplotype 1 carries the first allele of its `GT`.
    #[test]
    fn an_unphased_heterozygote_anchors_a_set_of_its_own() {
        let columns = [free(het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]))];

        let assignment = assign(&columns);

        assert_eq!(placed(&assignment, 1), Some((1001, Haplotype::First)), "the reference allele");
        assert_eq!(placed(&assignment, 2), Some((1001, Haplotype::Second)));
    }

    /// Only a site the segment emits anchors, so every set is one the file
    /// names.
    #[test]
    fn a_site_outside_the_core_anchors_nothing() {
        let (pileups, phases) =
            split_input(&[free(het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]))]);

        let assignment = assign_haplotypes(&pileups, &phases, &|pos| pos != 1000);

        assert!(assignment.votes.is_empty());
    }

    /// The pair-HMM's read assignments at a het indel place fragments exactly
    /// as an SNV's bases do, and an unphased one names its own set.
    #[test]
    fn a_heterozygous_indel_anchors_the_split() {
        // Fragments 1000-1005 favour the reference haplotype, 1010-1015 the
        // insertion; see `confirming`.
        let indel = with_indels(
            column(1000, A, &[], Some(GenotypeTag::HomRef), &[]),
            &[insertion(&[C])],
            &confirming([0, 1]),
        );
        let methylated: Vec<_> = (1000..1003)
            .map(|fragment| obs(fragment, T, Strand::OT))
            .chain((1010..1013).map(|fragment| obs(fragment, C, Strand::OT)))
            .collect();
        let columns = [free(indel), methylation_candidate(1040, &methylated)];

        let annotated = annotate(&columns);

        assert_eq!(annotated.set(1), Some(1001));
        assert_eq!(
            annotated.split(1).map(|split| split.counts()),
            Some([0, 3, 3, 0]),
            "the reference haplotype is methylated, the insertion's is not"
        );
    }

    /// Frames are independent, so a fragment reaching two sets has a haplotype
    /// in each, whether the second is a block or a singleton.
    #[test]
    fn a_fragment_reaching_two_sets_is_placed_in_each() {
        let second_block = phased(het(1040, A, G, &[obs(1, A, Strand::OT)]), 1041, 1);
        let singleton = free(het(1040, A, G, &[obs(1, A, Strand::OT)]));
        for (second, expected) in [(second_block, Haplotype::Second), (singleton, Haplotype::First)]
        {
            let columns = [phased(het(1000, A, G, &[obs(1, A, Strand::OT)]), 1001, 0), second];

            assert_eq!(
                all_placed(&assign(&columns), 1),
                [(1001, Haplotype::First), (1041, expected)]
            );
        }
    }

    /// Two reads of one template disagreeing is a sequencing error or a qname
    /// collision; a tie leaves the template with no haplotype.
    #[test]
    fn a_self_contradicting_fragment_is_dropped() {
        let columns = [
            phased(het(1000, A, G, &[graded(1, A, Strand::OT, 30)]), 1001, 0),
            phased(het(1040, A, G, &[graded(1, G, Strand::OT, 30)]), 1001, 0),
        ];
        assert_eq!(placed(&assign(&columns), 1), None);

        let decided = [
            phased(het(1000, A, G, &[graded(1, A, Strand::OT, 30)]), 1001, 0),
            phased(het(1040, A, G, &[graded(1, G, Strand::OT, 10)]), 1001, 0),
        ];
        assert_eq!(placed(&assign(&decided), 1), Some((1001, Haplotype::First)));
    }

    /// A C/T heterozygote in a CpG places every bottom-strand template but
    /// only the `C`-carrying top-strand ones, so on the top strand its split
    /// would be one-sided by construction; the bottom-strand column is
    /// unaffected, and a block keeps both.
    #[test]
    fn a_confounded_singleton_only_anchors_the_other_strands_column() {
        let anchor = |pos| {
            cpg(het(
                pos,
                C,
                T,
                &[obs(1, C, Strand::OB), obs(2, T, Strand::OB), obs(9, T, Strand::OB)],
            ))
        };
        let columns = [
            free(anchor(1000)),
            bottom_strand_candidate(1040, &[obs(1, G, Strand::OB), obs(2, A, Strand::OB)]),
            methylation_candidate(1080, &[obs(1, C, Strand::OT), obs(2, T, Strand::OT)]),
        ];

        let annotated = annotate(&columns);
        assert_eq!(annotated.split_columns(), 1);
        let split = annotated.split(1).expect("the bottom-strand column is resolved");
        assert_eq!(split.counts(), [1, 0, 0, 1], "the C haplotype is unconverted here");
        assert_eq!(annotated.split(2), None, "and the top-strand column is not");

        let blocked = [
            phased(anchor(1000), 1001, 0),
            methylation_candidate(1080, &[obs(1, C, Strand::OT), obs(2, T, Strand::OT)]),
        ];
        assert!(annotate(&blocked).split(1).is_some(), "the same reads, anchored by a block");
    }

    /// Covering a column equally well, a block beats a singleton: its
    /// haplotypes are the ones `PS` names.
    #[test]
    fn a_block_wins_a_depth_tie_against_a_singleton() {
        let columns = [
            free(het(960, A, G, &[obs(3, A, Strand::OT), obs(4, G, Strand::OT)])),
            phased(het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]), 1001, 0),
            methylation_candidate(
                1040,
                &[
                    obs(1, T, Strand::OT),
                    obs(2, C, Strand::OT),
                    obs(3, T, Strand::OT),
                    obs(4, C, Strand::OT),
                ],
            ),
        ];

        let split = annotate(&columns).split(2).expect("split");
        assert_eq!(split.set.as_i32(), 1001, "the block, not the singleton at 961");
        assert_eq!(split.depth(), 2, "and only the two fragments it placed");
    }

    /// A deeper set that only reaches the label strand cannot fill `HPC`, so a
    /// shallower one reaching the confounded strand takes the column.
    #[test]
    fn a_set_that_can_fill_hpc_outranks_a_deeper_one() {
        let columns = [
            phased(
                het(
                    960,
                    A,
                    C,
                    &[
                        obs(3, A, Strand::OB),
                        obs(4, A, Strand::OB),
                        obs(5, C, Strand::OB),
                        obs(6, C, Strand::OB),
                    ],
                ),
                961,
                0,
            ),
            phased(het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]), 1001, 0),
            methylation_candidate(
                1040,
                &[
                    obs(1, T, Strand::OT),
                    obs(2, C, Strand::OT),
                    obs(3, C, Strand::OB),
                    obs(4, C, Strand::OB),
                    obs(5, C, Strand::OB),
                    obs(6, C, Strand::OB),
                ],
            ),
        ];

        let split = annotate(&columns).split(2).expect("split");
        assert_eq!(split.set.as_i32(), 1001, "not the set at 961, which placed four reads");
        assert_eq!(split.counts(), [0, 1, 1, 0]);
    }

    /// The set that placed more of the column's reads takes it.
    #[test]
    fn the_better_covering_set_takes_the_column() {
        let columns = [
            free(het(960, A, G, &[obs(3, A, Strand::OT), obs(4, G, Strand::OT)])),
            phased(
                het(
                    1000,
                    A,
                    G,
                    &[
                        obs(1, A, Strand::OT),
                        obs(2, G, Strand::OT),
                        obs(5, A, Strand::OT),
                        obs(6, G, Strand::OT),
                    ],
                ),
                1001,
                0,
            ),
            methylation_candidate(
                1040,
                &[
                    obs(1, T, Strand::OT),
                    obs(2, C, Strand::OT),
                    obs(3, T, Strand::OT),
                    obs(4, C, Strand::OT),
                    obs(5, T, Strand::OT),
                    obs(6, C, Strand::OT),
                ],
            ),
        ];

        let split = annotate(&columns).split(2).expect("split");
        assert_eq!(split.set.as_i32(), 1001, "the block reached four of the six fragments");
        assert_eq!(split.depth(), 4);
    }

    #[test]
    fn a_compound_het_places_both_of_its_alts() {
        let site = phased(
            column(
                1000,
                A,
                &[(G, AltCall::RealVariant), (C, AltCall::RealVariant)],
                Some(GenotypeTag::AltHet(ALT_1, ALT_2)),
                &[obs(1, G, Strand::OT), obs(2, C, Strand::OT), obs(3, A, Strand::OT)],
            ),
            1001,
            2,
        );

        let assignment = assign(&[site]);

        assert_eq!(placed(&assignment, 1), Some((1001, Haplotype::Second)), "alt 1 is haplotype 2");
        assert_eq!(placed(&assignment, 2), Some((1001, Haplotype::First)), "alt 2 was named");
        assert_eq!(placed(&assignment, 3), None, "the reference is not an allele here");
    }

    /// A reference C whose T reads all come from one chromosome, next to a
    /// heterozygote that says which, phased or not; unphased, the set is the
    /// heterozygote's own 1-based position and its reference is haplotype 1.
    #[test]
    fn a_confounded_column_splits_its_top_strand_reads() {
        for anchor in [block_anchor(), singleton_anchor()] {
            let columns = [
                anchor,
                methylation_candidate(
                    1040,
                    &[
                        obs(1, T, Strand::OT),
                        obs(2, T, Strand::OT),
                        obs(3, C, Strand::OT),
                        obs(4, C, Strand::OT),
                    ],
                ),
            ];

            let annotated = annotate(&columns);
            assert_eq!(annotated.split_columns(), 1);

            let split = annotated.split(1).expect("the confounded column is split");
            assert_eq!(split.set.as_i32(), 1001);
            assert_eq!(split.counts(), [0, 2, 2, 0], "every T on haplotype 1, every C on 2");
            assert_eq!(annotated.presence(1), Some([CpgPresence::Present, CpgPresence::Present]));
        }
    }

    /// The other strand reads the same alleles without ambiguity.
    #[test]
    fn only_the_confounded_strand_is_counted() {
        let columns = [
            phased(het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]), 1001, 0),
            methylation_candidate(1040, &[obs(1, C, Strand::OB), obs(2, A, Strand::OB)]),
        ];

        assert_eq!(annotate(&columns).split(1), None);
    }

    /// An unmethylated CpG is the null; leaving it out would make the surface a
    /// sample of the methylated positions.
    #[test]
    fn an_unmethylated_cpg_is_split_with_no_converted_reads() {
        let columns = [
            free(het(
                1000,
                A,
                C,
                &[
                    obs(1, A, Strand::OT),
                    obs(2, A, Strand::OT),
                    obs(3, C, Strand::OT),
                    obs(4, C, Strand::OT),
                ],
            )),
            free(cpg(column(
                1040,
                C,
                &[],
                Some(GenotypeTag::HomRef),
                &[
                    obs(1, C, Strand::OT),
                    obs(2, C, Strand::OT),
                    obs(3, C, Strand::OT),
                    obs(4, C, Strand::OT),
                ],
            ))),
        ];

        let annotated = annotate(&columns);
        assert_eq!(annotated.split_columns(), 1);
        let split = annotated.split(1).expect("an unmethylated CpG is still split");
        assert_eq!(split.counts(), [2, 0, 2, 0], "two unconverted reads on each haplotype");
    }

    /// Outside a CpG only a called `C>T` has anything to resolve; `T` reads
    /// the caller put down to errors are not one.
    #[test]
    fn a_reference_c_outside_a_cpg_needs_a_called_partner_allele() {
        let columns = [
            free(het(1000, A, C, &[obs(1, A, Strand::OT), obs(2, C, Strand::OT)])),
            free(column(
                1040,
                C,
                &[(T, AltCall::ReadError)],
                Some(GenotypeTag::HomRef),
                &[obs(1, T, Strand::OT), obs(2, C, Strand::OT)],
            )),
        ];

        assert_eq!(annotate(&columns).split(1), None);
    }

    /// A split with no placed reads would read as evidence of an even spread.
    #[test]
    fn a_column_no_placed_fragment_reaches_is_not_split() {
        let columns = [
            phased(het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]), 1001, 0),
            methylation_candidate(1040, &[obs(7, T, Strand::OT), obs(8, T, Strand::OT)]),
        ];

        assert_eq!(annotate(&columns).split(1), None);
    }

    /// Methylation, told apart from a variant: converted reads on both
    /// chromosomes.
    #[test]
    fn methylation_spreads_across_both_haplotypes() {
        let columns = [
            block_anchor(),
            methylation_candidate(
                1040,
                &[
                    obs(1, T, Strand::OT),
                    obs(2, C, Strand::OT),
                    obs(3, T, Strand::OT),
                    obs(4, C, Strand::OT),
                ],
            ),
        ];

        assert_eq!(annotate(&columns).split(1).map(|split| split.counts()), Some([1, 1, 1, 1]));
    }

    /// A heterozygote's own alleles separate the reads crossing it, so its own
    /// votes are left out when it is scored: alone, it places none of its reads.
    #[test]
    fn a_heterozygote_alone_places_none_of_its_own_reads() {
        let columns = [free(het(
            1000,
            A,
            G,
            &[
                obs(1, A, Strand::OT),
                obs(2, A, Strand::OT),
                obs(3, G, Strand::OB),
                obs(4, G, Strand::OB),
            ],
        ))];

        assert_eq!(annotate(&columns).split_columns(), 0, "its own reads are not its own evidence");
        assert_eq!(assign(&columns).votes.len(), 4, "the votes are still there");
    }

    /// Here the block's other heterozygote places fragment 2, which this
    /// column's own vote would have cancelled.
    #[test]
    fn a_confounded_column_is_left_out_of_its_own_split() {
        let columns = [
            phased(het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]), 1001, 0),
            phased(het(1040, C, T, &[obs(1, T, Strand::OT), obs(2, C, Strand::OT)]), 1001, 0),
        ];

        let annotated = annotate(&columns);
        assert_eq!(annotated.split_columns(), 1);
        assert_eq!(
            annotated.split(1).map(|split| split.counts()),
            Some([0, 1, 1, 0]),
            "the T on haplotype 1 and the C on haplotype 2, both said by the column at 1000"
        );
    }

    /// Fragments 1-4 and 9-10 are read on the top strand at 1040, 5-8 on the
    /// bottom. The anchor puts 1-4 and 7-8 on haplotype 1 (`A`), the rest on
    /// haplotype 2 (`C`).
    fn heterozygous_c_to_t_cpg() -> [Column; 2] {
        [
            labelled_anchor(
                &[
                    (1, Strand::OT),
                    (2, Strand::OT),
                    (3, Strand::OT),
                    (4, Strand::OT),
                    (7, Strand::OB),
                    (8, Strand::OB),
                ],
                &[(5, Strand::OB), (6, Strand::OB), (9, Strand::OT), (10, Strand::OT)],
            ),
            free(cpg(het(
                1040,
                C,
                T,
                &[
                    obs(1, T, Strand::OT),
                    obs(2, T, Strand::OT),
                    obs(3, C, Strand::OT),
                    obs(4, C, Strand::OT),
                    obs(7, C, Strand::OB),
                    obs(8, C, Strand::OB),
                    obs(5, T, Strand::OB),
                    obs(6, T, Strand::OB),
                    obs(9, T, Strand::OT),
                    obs(10, T, Strand::OT),
                ],
            ))),
        ]
    }

    /// A heterozygous `C>T` CpG: the bottom strand, where the alleles read
    /// apart, says which haplotype has the `T` and so no cytosine. At a `G>A`
    /// the label comes off the top strand.
    #[test]
    fn a_heterozygous_cpg_variant_marks_the_haplotype_that_lost_its_cytosine() {
        for columns in [heterozygous_c_to_t_cpg().to_vec(), mirror(&heterozygous_c_to_t_cpg())] {
            let annotated = annotate(&columns);

            let split = annotated.split(1).expect("split");
            assert_eq!(split.counts(), [2, 2, 0, 2], "the counts stay raw");
            assert_eq!(annotated.presence(1), Some([CpgPresence::Present, CpgPresence::Absent]));
        }
    }

    /// A de-novo CpG sits on the *alt* chromosome; the same rule finds it
    /// because it asks the pair for its unconverted base, at a `T>C` and at
    /// its `A>G` mirror.
    #[test]
    fn a_de_novo_cpg_is_present_on_the_haplotype_that_creates_it() {
        let columns = [
            labelled_anchor(
                &[(1, Strand::OT), (2, Strand::OT), (7, Strand::OB), (8, Strand::OB)],
                &[
                    (3, Strand::OT),
                    (4, Strand::OT),
                    (5, Strand::OB),
                    (6, Strand::OB),
                    (9, Strand::OT),
                    (10, Strand::OT),
                ],
            ),
            free(denovo(het(
                1040,
                T,
                C,
                &[
                    obs(1, T, Strand::OT),
                    obs(2, T, Strand::OT),
                    obs(7, T, Strand::OB),
                    obs(8, T, Strand::OB),
                    obs(5, C, Strand::OB),
                    obs(6, C, Strand::OB),
                    obs(3, T, Strand::OT),
                    obs(4, T, Strand::OT),
                    obs(9, C, Strand::OT),
                    obs(10, C, Strand::OT),
                ],
            ))),
        ];

        for columns in [columns.to_vec(), mirror(&columns)] {
            let annotated = annotate(&columns);

            assert_eq!(annotated.split(1).map(|split| split.counts()), Some([0, 2, 2, 2]));
            assert_eq!(annotated.presence(1), Some([CpgPresence::Absent, CpgPresence::Present]));
        }
    }

    #[test]
    fn a_homozygous_reference_cpg_has_a_cytosine_on_both_haplotypes() {
        let columns = [
            labelled_anchor(
                &[(1, Strand::OT), (2, Strand::OT)],
                &[(3, Strand::OT), (4, Strand::OT)],
            ),
            methylation_candidate(
                1040,
                &[
                    obs(1, T, Strand::OT),
                    obs(2, C, Strand::OT),
                    obs(3, T, Strand::OT),
                    obs(4, T, Strand::OT),
                ],
            ),
        ];

        let annotated = annotate(&columns);

        assert_eq!(annotated.presence(1), Some([CpgPresence::Present, CpgPresence::Present]));
        assert_eq!(annotated.split(1).map(|split| split.counts()), Some([1, 1, 0, 2]));
    }

    /// The `G` a base later is a `G` on both chromosomes, and still has no CpG
    /// on either.
    #[test]
    fn a_homozygous_c_to_t_destroys_the_cpg_on_both_haplotypes() {
        let columns = [
            labelled_anchor(
                &[(1, Strand::OT), (2, Strand::OT), (7, Strand::OB)],
                &[(3, Strand::OT), (4, Strand::OT), (8, Strand::OB)],
            ),
            free(cpg(column(
                1040,
                C,
                &[(T, AltCall::RealVariant)],
                Some(GenotypeTag::HomAlt(ALT_1)),
                &[
                    obs(1, T, Strand::OT),
                    obs(2, T, Strand::OT),
                    obs(3, T, Strand::OT),
                    obs(4, T, Strand::OT),
                ],
            ))),
            bottom_strand_candidate(1041, &[obs(7, G, Strand::OB), obs(8, A, Strand::OB)]),
        ];

        let annotated = annotate(&columns);

        let absent = Some([CpgPresence::Absent, CpgPresence::Absent]);
        assert_eq!(annotated.presence(1), absent);
        assert_eq!(annotated.presence(2), absent, "the partner G has no CpG either");
    }

    /// A `C>A` at a CpG gets no split of its own when nothing else places its
    /// reads, but it anchors its partner as a singleton; that partner is on a
    /// chromosome the `A` took the CpG from.
    #[test]
    fn a_heterozygous_cpg_without_a_split_leaves_its_partner_unresolved() {
        let columns = [
            free(cpg(het(
                1040,
                C,
                A,
                &[
                    obs(1, C, Strand::OT),
                    obs(2, A, Strand::OT),
                    obs(3, C, Strand::OB),
                    obs(4, A, Strand::OB),
                ],
            ))),
            bottom_strand_candidate(1041, &[obs(3, G, Strand::OB), obs(4, A, Strand::OB)]),
        ];

        let annotated = annotate(&columns);

        assert_eq!(annotated.split(0), None);
        assert_eq!(annotated.set(1), Some(1041), "the heterozygote anchors its partner");
        assert_eq!(annotated.presence(1), Some([CpgPresence::Unresolved, CpgPresence::Unresolved]));
    }

    /// One read each way on the label strand is a coin toss, and one read is
    /// below the two the label needs; the raw counts are written anyway. The
    /// `C>T` has still destroyed the CpG on one chromosome, just not a known
    /// one, so its `G` must not claim both either.
    #[test]
    fn a_label_strand_that_cannot_decide_leaves_both_haplotypes_unresolved() {
        let contested = |first: Base, second: Base| {
            [
                labelled_anchor(
                    &[(1, Strand::OT), (7, Strand::OB)],
                    &[(5, Strand::OB), (9, Strand::OT)],
                ),
                free(cpg(het(
                    1040,
                    C,
                    T,
                    &[
                        obs(1, T, Strand::OT),
                        obs(9, T, Strand::OT),
                        obs(7, first, Strand::OB),
                        obs(5, second, Strand::OB),
                    ],
                ))),
                bottom_strand_candidate(1041, &[obs(7, G, Strand::OB), obs(5, A, Strand::OB)]),
            ]
        };
        let unresolved = Some([CpgPresence::Unresolved, CpgPresence::Unresolved]);

        let split_evenly = annotate(&contested(T, T));
        assert_eq!(split_evenly.presence(1), unresolved, "one T on each haplotype names neither");
        assert_eq!(split_evenly.set(2), split_evenly.set(1));
        assert_eq!(split_evenly.presence(2), unresolved, "and the partner inherits that");

        let single_read = annotate(&contested(C, T));
        assert_eq!(single_read.presence(1), unresolved, "one read is below the two needed");
    }

    /// A CpG's `G` a base after a heterozygous `C>T` is a `G` on both
    /// chromosomes, and learns from its partner which one lost the CpG; a
    /// partner in another set learns nothing.
    #[test]
    fn a_heterozygous_cpg_tells_its_partner_position_which_haplotype_lost_it() {
        let partner = |extra: &[PhaseObservation]| {
            let mut observations = vec![
                obs(7, G, Strand::OB),
                obs(8, A, Strand::OB),
                obs(5, G, Strand::OB),
                obs(6, A, Strand::OB),
            ];
            observations.extend_from_slice(extra);
            bottom_strand_candidate(1041, &observations)
        };

        let [anchor, heterozygous] = heterozygous_c_to_t_cpg();
        let one_set = annotate(&[anchor, heterozygous, partner(&[])]);

        assert_eq!(one_set.presence(1), Some([CpgPresence::Present, CpgPresence::Absent]));
        assert_eq!(one_set.set(2), one_set.set(1));
        assert_eq!(
            one_set.presence(2),
            Some([CpgPresence::Present, CpgPresence::Absent]),
            "the G has no CpG on the chromosome whose C became a T"
        );

        let [anchor, heterozygous] = heterozygous_c_to_t_cpg();
        let two_sets = annotate(&[
            anchor,
            heterozygous,
            partner(&[obs(11, G, Strand::OB), obs(12, A, Strand::OB)]),
            free(het(
                1080,
                A,
                C,
                &[
                    obs(7, A, Strand::OB),
                    obs(8, A, Strand::OB),
                    obs(11, A, Strand::OB),
                    obs(5, C, Strand::OB),
                    obs(6, C, Strand::OB),
                    obs(12, C, Strand::OB),
                ],
            )),
        ]);

        assert_ne!(two_sets.set(2), two_sets.set(1));
        assert_eq!(
            two_sets.presence(2),
            Some([CpgPresence::Unresolved, CpgPresence::Unresolved]),
            "copying across sets would name a chromosome at random"
        );
    }

    /// Outside a CpG a called `T>C` is split on the same C/T pair as a `C>T`,
    /// its reference being the converted base, and neither haplotype has a CpG.
    #[test]
    fn a_t_to_c_column_is_split_on_the_c_t_pair() {
        let columns = [
            phased(het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]), 1001, 0),
            free(column(
                1040,
                T,
                &[(C, AltCall::RealVariant)],
                None,
                &[obs(1, C, Strand::OT), obs(2, T, Strand::OT)],
            )),
        ];

        let annotated = annotate(&columns);

        assert_eq!(annotated.split(1).map(|split| split.counts()), Some([1, 0, 0, 1]));
        assert_eq!(annotated.presence(1), Some([CpgPresence::Absent, CpgPresence::Absent]));
    }

    /// A de-novo CpG's partner `G` is kept for phasing only where its reads
    /// look heterozygous, so a split there would be a selected sample.
    #[test]
    fn a_de_novo_cpg_partner_is_not_split_methylated_or_not() {
        let partner = |bases: [Base; 4]| {
            let observations: Vec<_> =
                (1..).zip(bases).map(|(fragment, base)| obs(fragment, base, Strand::OB)).collect();
            let mut column = column(1041, G, &[], Some(GenotypeTag::HomRef), &observations);
            column.tags.denovo_cpg_partner = true;
            free(column)
        };

        for bases in [[G, G, G, G], [A, A, G, G]] {
            let columns = [block_anchor(), partner(bases)];
            assert_eq!(annotate(&columns).split(1), None, "{bases:?}");
        }
    }

    /// Each molecule is seen at two CpG columns of one set, so `HPC` counts it
    /// twice and `HPF` once; a molecule whose reads tie does not vote.
    #[test]
    fn a_fragment_votes_once_per_set_by_the_majority_of_its_cpgs() {
        let columns = [
            block_anchor(),
            methylation_candidate(
                1040,
                &[
                    obs(1, T, Strand::OT),
                    obs(2, C, Strand::OT),
                    obs(3, T, Strand::OT),
                    obs(4, C, Strand::OT),
                ],
            ),
            methylation_candidate(
                1080,
                &[
                    obs(1, T, Strand::OT),
                    obs(2, C, Strand::OT),
                    obs(3, C, Strand::OT),
                    obs(4, C, Strand::OT),
                ],
            ),
        ];

        let annotated = annotate(&columns);
        assert_eq!(annotated.split_columns(), 2);

        let split = annotated.split(1).expect("split");
        assert_eq!(split.counts(), [1, 1, 1, 1]);
        assert_eq!(split.fragment_counts(), [1, 1, 1, 0]);
        assert_eq!(
            annotated.split(2).map(|split| split.fragment_counts()),
            Some([1, 1, 1, 0]),
            "the count belongs to the set"
        );
    }

    /// A haplotype with no cytosine converts every read by construction, and
    /// the other haplotype's molecules alone would make the set's table test
    /// which chromosome has the CpG, so neither side of the column votes.
    #[test]
    fn a_cpg_on_one_haplotype_adds_no_fragments() {
        let annotated = annotate(&heterozygous_c_to_t_cpg());

        let split = annotated.split(1).expect("the column still reports its reads");
        assert_eq!(split.fragment_counts(), [0, 0, 0, 0]);
    }

    #[test]
    fn a_fragment_reaching_two_sets_is_counted_in_both() {
        let columns = [
            phased(het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]), 1001, 0),
            methylation_candidate(1020, &[obs(1, T, Strand::OT), obs(2, C, Strand::OT)]),
            phased(het(1040, A, G, &[obs(1, A, Strand::OT), obs(3, G, Strand::OT)]), 1041, 0),
            methylation_candidate(1060, &[obs(1, T, Strand::OT), obs(3, C, Strand::OT)]),
        ];

        let annotated = annotate(&columns);

        let counts = |index| annotated.split(index).map(|s| (s.set.as_i32(), s.fragment_counts()));
        assert_eq!(counts(1), Some((1001, [0, 1, 1, 0])));
        assert_eq!(counts(3), Some((1041, [0, 1, 1, 0])));
    }
}
