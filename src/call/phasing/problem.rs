//! Turning a segment's columns into a phasing problem.
//!
//! Two lists come out of a segment: the *sites*, heterozygous columns whose two
//! alleles have to be split across haplotypes, and the *fragments*, templates
//! whose alleles at several sites share a chromosome. A fragment carries one
//! bit per site saying which allele it saw; at a C/T or G/A heterozygote whose
//! cytosine is in a CpG, observations TAPS conversion alone explains are
//! dropped, or methylation would read as the other allele.

use super::{
    Cost, FragmentId, GtAllele, MIN_ALT_READS, MIN_MAPQ, PhaseObservation, SiteKind,
    observations_of,
};
use crate::{
    call::{
        pair_hmm::{HapIdx, Margin, ReadScores},
        variant_calling::{GenotypeTag, indel_calling::IndelRecords},
    },
    metrics::{
        AltCall, PileupMetrics,
        methylation::{CpgSide, cpg_origin},
    },
};
use rustc_hash::FxHashMap;
use seqair_types::{Base, Pos0, SmallVec, Strand};
use std::{
    num::NonZeroU8,
    ops::{Index, IndexMut},
};
use tracing::{debug, warn};

/// A site of one segment's [`Sites`].
///
/// `u32` rather than `usize` keeps an [`AlleleSupport`] at eight bytes, so a
/// fragment's supports stay inline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SiteIndex(u32);

impl SiteIndex {
    /// A token for a problem a test assembled by hand.
    #[cfg(test)]
    pub(crate) const fn new(index: u32) -> Self {
        Self(index)
    }

    #[must_use]
    pub const fn get(self) -> usize {
        self.0 as usize
    }
}

/// A segment's heterozygous sites, and the only thing that mints a
/// [`SiteIndex`].
pub type Sites = SiteMap<PhaseSite>;

/// A value per site, indexed by [`SiteIndex`].
///
/// Indexing is total by construction rather than by type: tokens come only from
/// [`SiteMap::push`] on the segment's [`Sites`], which is only appended to, and
/// every other map is shaped after it. Nothing stops a token from one segment
/// reaching another segment's map; nothing in the pipeline does that either.
#[derive(Debug, Clone)]
pub struct SiteMap<T>(Vec<T>);

// Derived, it would demand `T: Default`.
impl<T> Default for SiteMap<T> {
    fn default() -> Self {
        Self(Vec::new())
    }
}

impl<T> SiteMap<T> {
    /// A map with one `value` per site of `shape`.
    pub fn filled<U>(shape: &SiteMap<U>, value: T) -> Self
    where
        T: Clone,
    {
        Self(vec![value; shape.len()])
    }

    pub fn from_fn<U>(shape: &SiteMap<U>, f: impl FnMut(SiteIndex) -> T) -> Self {
        Self(shape.indices().map(f).collect())
    }

    /// `None` once a segment holds more sites than a [`SiteIndex`] can name.
    fn push(&mut self, value: T) -> Option<SiteIndex> {
        let index = SiteIndex(u32::try_from(self.0.len()).ok()?);
        self.0.push(value);
        Some(index)
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[cfg(test)]
    pub(crate) fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    #[cfg(test)]
    pub(crate) fn iter(&self) -> std::slice::Iter<'_, T> {
        self.0.iter()
    }

    /// Every site, ascending.
    pub fn indices(&self) -> impl Iterator<Item = SiteIndex> + use<T> {
        // `push` refuses any index past `u32::MAX`, so this drops nothing.
        (0..self.0.len()).filter_map(|index| u32::try_from(index).ok().map(SiteIndex))
    }
}

/// A problem a test assembled by hand.
#[cfg(test)]
impl<T> From<Vec<T>> for SiteMap<T> {
    fn from(values: Vec<T>) -> Self {
        Self(values)
    }
}

impl<T> Index<SiteIndex> for SiteMap<T> {
    type Output = T;

    #[expect(
        clippy::indexing_slicing,
        reason = "every map is shaped after the Sites that minted the token"
    )]
    fn index(&self, site: SiteIndex) -> &T {
        &self.0[site.get()]
    }
}

impl<T> IndexMut<SiteIndex> for SiteMap<T> {
    #[expect(
        clippy::indexing_slicing,
        reason = "every map is shaped after the Sites that minted the token"
    )]
    fn index_mut(&mut self, site: SiteIndex) -> &mut T {
        &mut self.0[site.get()]
    }
}

/// One of the two alleles a heterozygous site carries, as a position in
/// [`PhaseSite::alleles`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiteAllele {
    First,
    Second,
}

impl SiteAllele {
    #[must_use]
    pub const fn flipped(self) -> Self {
        match self {
            Self::First => Self::Second,
            Self::Second => Self::First,
        }
    }
}

/// One allele of a heterozygous SNV genotype, and the base a read shows for it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HetAllele {
    pub allele: GtAllele,
    pub base: Base,
}

/// A heterozygous column's two alleles, ascending.
///
/// Every allele must be a called real variant with a known base, and the two
/// must read differently, or an observation could not be told apart.
fn het_alleles(pileup: &PileupMetrics) -> Option<[HetAllele; 2]> {
    let alt = |n: NonZeroU8| {
        let alt = pileup.alts.get(usize::from(n.get()) - 1)?;
        (alt.call == AltCall::RealVariant && alt.base.known_index().is_some())
            .then_some(HetAllele { allele: GtAllele::Alt(n), base: alt.base })
    };

    let alleles = match pileup.pos_metrics.extended.genotype?.genotype {
        GenotypeTag::RefHet(n) => {
            let reference = pileup.reference_base;
            reference.known_index()?;
            [HetAllele { allele: GtAllele::Ref, base: reference }, alt(n)?]
        }
        GenotypeTag::AltHet(m, n) => {
            debug_assert!(m < n, "a compound het names its alts ascending");
            [alt(m)?, alt(n)?]
        }
        GenotypeTag::HomRef | GenotypeTag::HomAlt(_) => return None,
    };

    let [first, second] = alleles;
    (first.base != second.base).then_some(alleles)
}

/// A heterozygous record a solver can phase.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PhaseSite {
    /// Where the column sits in the slice [`build_problem`] was given.
    pub pileup: usize,
    /// Which of the column's records this is.
    pub kind: SiteKind,
    pub pos: Pos0,
    /// The genotype's two alleles, ascending.
    pub alleles: [GtAllele; 2],
    /// Why the site may join a block but not link two, if it is one of those.
    pub guard: Option<Guard>,
}

/// Why a site is kept out of the solve and only attached to a block after it:
/// each marks sites whose reads are wrong together often enough that one
/// such site bridging two blocks is a switch error between them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Guard {
    /// One of the SNV's alleles is carried by reads of RMS MAPQ below
    /// [`MIN_ALLELE_MAPQ`].
    AlleleMapq,
    /// Over [`MAX_DEPTH_RATIO`] times the depth expected of the segment.
    Depth,
    /// [`MAX_NEIGHBOURS`] or more other sites within [`NEIGHBOURHOOD`] bases.
    HetDensity,
}

/// RMS MAPQ both of an SNV's alleles need: guarding below it halves the
/// switch errors outside chr12.
const MIN_ALLELE_MAPQ: f64 = 50.0;
/// Depth, relative to the segment's expected depth, above which a site is
/// likely in a collapsed repeat.
const MAX_DEPTH_RATIO: f64 = 2.0;

/// Bases either side of a site [`Guard::HetDensity`] looks.
const NEIGHBOURHOOD: u64 = 100;
/// Other sites within [`NEIGHBOURHOOD`] that make a site's region too dense
/// to trust with a link.
const MAX_NEIGHBOURS: usize = 5;

/// The run's mapped reads per base over the segment's.
///
/// Scales the segment's median depth to the depth expected of it, so a
/// segment lying wholly inside a collapsed repeat is judged by the run's
/// depth and not its own. Without it (no index counts) the median is used
/// as it is, which misses repeats longer than a segment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DepthScale(pub f64);

/// Each site's [`Guard`], the first that applies; it needs the whole segment's
/// sites, for the median depth and the neighbours.
fn guards(
    sites: &Sites,
    pileups: &[PileupMetrics],
    scale: Option<DepthScale>,
) -> SiteMap<Option<Guard>> {
    let depth = |site: &PhaseSite| pileups.get(site.pileup).map_or(0, |p| p.pos_metrics.depth);
    let mut depths: Vec<u32> = sites.0.iter().map(depth).collect();
    depths.sort_unstable();
    let median = depths.get(depths.len() / 2).copied().unwrap_or(0);
    let expected = f64::from(median) * scale.map_or(1.0, |DepthScale(scale)| scale);
    let positions: Vec<u64> = sites.0.iter().map(|site| site.pos.as_u64()).collect();

    SiteMap::from_fn(sites, |index| {
        let site = &sites[index];
        let pileup = pileups.get(site.pileup)?;
        let allele_mapq = |allele: GtAllele| match allele {
            GtAllele::Ref => Some(*pileup.ref_metrics.mapq),
            GtAllele::Alt(n) => {
                pileup.alts.get(usize::from(n.get()) - 1).map(|alt| *alt.metrics.mapq)
            }
        };
        let pos = site.pos.as_u64();
        // The window holds the site itself.
        let neighbours = positions.partition_point(|&p| p <= pos + NEIGHBOURHOOD).saturating_sub(
            positions.partition_point(|&p| p < pos.saturating_sub(NEIGHBOURHOOD)) + 1,
        );

        if site.kind == SiteKind::Snv
            && site
                .alleles
                .into_iter()
                .any(|allele| allele_mapq(allele).is_none_or(|mapq| mapq < MIN_ALLELE_MAPQ))
        {
            Some(Guard::AlleleMapq)
        } else if f64::from(depth(site)) > MAX_DEPTH_RATIO * expected {
            Some(Guard::Depth)
        } else if neighbours >= MAX_NEIGHBOURS {
            Some(Guard::HetDensity)
        } else {
            None
        }
    })
}

impl PhaseSite {
    #[must_use]
    pub const fn allele(&self, allele: SiteAllele) -> GtAllele {
        let [first, second] = self.alleles;
        match allele {
            SiteAllele::First => first,
            SiteAllele::Second => second,
        }
    }
}

/// One read's say at a site: which allele its fragment shows there, and what
/// a haplotype is charged for disagreeing.
#[derive(Debug, Clone, Copy)]
struct Vote {
    fragment: FragmentId,
    allele: SiteAllele,
    cost: Cost,
}

/// An SNV site's alleles as the bases reads show for them.
struct BaseAlleles {
    bases: [Base; 2],
    /// The CpG cytosine one allele is, when conversion makes it read as the
    /// other: a C/T site's top strand or a G/A site's bottom strand.
    methylable: Option<CpgSide>,
}

impl BaseAlleles {
    fn of(pileup: &PileupMetrics, alleles: [HetAllele; 2]) -> Self {
        let bases = alleles.map(|allele| allele.base);
        // Outside a CpG, cytosine methylation is too rare in mammals to give
        // up a strand's reads for.
        let methylable = [CpgSide::C, CpgSide::G].into_iter().find(|&side| {
            bases.contains(&side.unmod_base())
                && bases.contains(&side.mod_base())
                && cpg_origin(pileup, side).is_some()
        });
        Self { bases, methylable }
    }

    fn vote(&self, observation: &PhaseObservation) -> Option<Vote> {
        let allele = self.allele_of(observation.base)?;
        if self.taps_confounds(observation.base, observation.strand) {
            return None;
        }
        Some(Vote { fragment: observation.fragment, allele, cost: Cost::of(observation.qual) })
    }

    /// Which allele an observed base is, or `None` when it is neither.
    fn allele_of(&self, base: Base) -> Option<SiteAllele> {
        match self.bases.map(|allele| allele == base) {
            [true, _] => Some(SiteAllele::First),
            [_, true] => Some(SiteAllele::Second),
            _ => None,
        }
    }

    /// Could conversion alone have produced this read? One with no strand
    /// counts, since nothing rules conversion out for it.
    fn taps_confounds(&self, base: Base, strand: Strand) -> bool {
        self.methylable.is_some_and(|side| {
            base == side.mod_base() && strand.ok().is_none_or(|strand| strand == side.strand())
        })
    }
}

/// A column's SNV record as a site, if it is a heterozygote with
/// [`MIN_ALT_READS`] votes on each allele, and the votes its reads' bases cast.
fn snv_site(pileup: &PileupMetrics, pileup_idx: usize) -> Option<(PhaseSite, Vec<Vote>)> {
    let alleles = het_alleles(pileup)?;
    let site = PhaseSite {
        pileup: pileup_idx,
        kind: SiteKind::Snv,
        pos: Pos0::new(pileup.pos)?,
        alleles: alleles.map(|allele| allele.allele),
        guard: None,
    };
    let observations = observations_of(pileup).filter(|observations| !observations.is_empty())?;
    let bases = BaseAlleles::of(pileup, alleles);
    let votes: Vec<Vote> =
        observations.iter().filter_map(|observation| bases.vote(observation)).collect();
    // A site with one read behind an allele is phased by that read's coin
    // toss, and is mostly a false het: the reference allele's reads place it.
    let backed = |allele| votes.iter().filter(|vote| vote.allele == allele).count();
    [SiteAllele::First, SiteAllele::Second]
        .into_iter()
        .all(|allele| backed(allele) >= MIN_ALT_READS as usize)
        .then_some((site, votes))
}

/// How sure the pair-HMM must be of an indel column's genotype (Phred, before
/// the cap at 99 it is written with) for its record to be a phase site. A
/// genotype it doubts is mostly a homozygote, an artefact or a stutter-ridden
/// tract: its reads' votes are near coin tosses, and one such site bridging two
/// blocks is a switch error between their SNVs, however weak its edges.
/// Measured on chr12 NA12878 with `MIN_READS_PER_INDEL_ALLELE`; see AGENTS.md.
const MIN_INDEL_GENOTYPE_QUALITY: f64 = 100.0;

/// Votes each allele of an indel site needs at [`Margin::PRESENT`] or more:
/// reads the pair-HMM finds ten times likelier under it than under the other,
/// that can vote (a fragment, the phasing MAPQ). The flat-prior genotype
/// quality treats every read's likelihoods as exact, so a few reads with
/// large gap penalties make a confident het out of thin evidence.
const MIN_READS_PER_INDEL_ALLELE: usize = 5;

impl Cost {
    /// What a vote by a read `margin` likelier under its allele costs.
    #[expect(clippy::cast_possible_truncation, reason = "a vote's ratio is an f32 too")]
    fn of_margin(margin: Margin) -> Self {
        Self::of_log10_ratio(margin.log10() as f32)
    }
}

/// A column's indel record as a site, and the votes its scored reads cast.
///
/// The pair-HMM's called genotype must be heterozygous, a `0/1` or the
/// compound het, and *be* its best pair at [`MIN_INDEL_GENOTYPE_QUALITY`] or
/// better, each allele carried by [`MIN_READS_PER_INDEL_ALLELE`] votes. A
/// `0/1` whose best pair holds a second, uncalled allele is no site: its
/// reads are split three ways. The alleles are numbered as the record writes
/// them, a compound het's `1/2` in call order.
fn indel_site(pileup: &PileupMetrics, pileup_idx: usize) -> Option<(PhaseSite, Vec<Vote>)> {
    let data = pileup.indel_data.as_deref()?;
    let scored = data.scored.as_ref()?;
    let scores = &scored.haplotypes;
    let verdict = &scores.verdict;
    let called = verdict.called_pair();
    if called != verdict.best || !called.is_heterozygous() {
        return None;
    }
    // The site's two alleles are the called pair's two haplotypes, named by
    // the record: the one carrying the event and the one not.
    let [(low, low_hap), (high, high_hap)] = match &scored.records {
        IndelRecords::Compound([first, second]) => {
            let GenotypeTag::AltHet(m, n) = first.genotype else { return None };
            let (first_hap, second_hap) = called.split(first.event)?;
            if !second_hap.carries(second.event) {
                return None;
            }
            [(GtAllele::Alt(m), first_hap), (GtAllele::Alt(n), second_hap)]
        }
        IndelRecords::Separate(calls) => {
            let call = calls.iter().find(|call| call.genotype.is_heterozygous())?;
            let GenotypeTag::RefHet(n) = call.genotype else { return None };
            let (carrier, other) = called.split(call.event)?;
            [(GtAllele::Ref, other), (GtAllele::Alt(n), carrier)]
        }
    };
    if *verdict.quality < MIN_INDEL_GENOTYPE_QUALITY {
        return None;
    }
    let votes: Vec<Vote> =
        scores.reads.iter().filter_map(|read| scored_vote(read, [low_hap, high_hap])).collect();
    let present = Cost::of_margin(Margin::PRESENT);
    let carried_by = |allele: SiteAllele| {
        votes.iter().filter(|vote| vote.allele == allele && vote.cost >= present).count()
    };
    if carried_by(SiteAllele::First).min(carried_by(SiteAllele::Second))
        < MIN_READS_PER_INDEL_ALLELE
    {
        return None;
    }
    let site = PhaseSite {
        pileup: pileup_idx,
        kind: SiteKind::Indel,
        pos: Pos0::new(pileup.pos)?,
        alleles: [low, high],
        guard: None,
    };
    Some((site, votes))
}

/// A read's vote for the haplotype the pair-HMM finds it likelier to come
/// from, charged by how much likelier. A read without a template, below the
/// phasing MAPQ, or with no preference worth a Phred point casts none. Base
/// quality is not asked: the pair-HMM has weighed it.
fn scored_vote(read: &ReadScores, [first, second]: [HapIdx; 2]) -> Option<Vote> {
    let fragment = read.fragment?;
    if read.mapq < MIN_MAPQ {
        return None;
    }
    let ratio = read.log10_ratio(first, second)?;
    let allele = if ratio > 0.0 { SiteAllele::First } else { SiteAllele::Second };
    let cost = Cost::of_log10_ratio(ratio);
    (cost > Cost(0)).then_some(Vote { fragment, allele, cost })
}

/// A fragment's allele at one site.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlleleSupport {
    pub site: SiteIndex,
    pub allele: SiteAllele,
    /// What a haplotype is charged for disagreeing with this observation.
    pub cost: Cost,
}

/// One template's alleles across the sites it covers; both mates share it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fragment {
    pub id: FragmentId,
    /// Ascending by site, at most one entry per site.
    pub supports: SmallVec<AlleleSupport, 4>,
}

/// A segment's phasing problem: what to phase, and the evidence for it.
#[derive(Debug, Clone, Default)]
pub struct PhaseProblem {
    /// Every heterozygous column that carried observations, ascending,
    /// including ones no fragment reaches.
    pub sites: Sites,
    /// Templates spanning at least two sites, in order of their first site.
    pub fragments: Vec<Fragment>,
}

/// Collect a segment's heterozygous sites and the fragments linking them.
#[must_use]
pub fn build_problem(pileups: &[PileupMetrics], scale: Option<DepthScale>) -> PhaseProblem {
    let mut builder = Builder::default();

    'columns: for (pileup_idx, pileup) in pileups.iter().enumerate() {
        let sites = [snv_site(pileup, pileup_idx), indel_site(pileup, pileup_idx)];
        for (site, votes) in sites.into_iter().flatten() {
            let Some(index) = builder.sites.push(site) else {
                warn!(
                    sites = builder.sites.len(),
                    "Too many heterozygous columns in one segment to phase; the rest stay unphased"
                );
                break 'columns;
            };
            for vote in votes {
                builder.record(index, vote);
            }
        }
    }

    let mut problem = builder.finish();
    let guards = guards(&problem.sites, pileups, scale);
    for index in problem.sites.indices() {
        problem.sites[index].guard = guards[index];
    }
    debug!(
        sites = problem.sites.len(),
        fragments = problem.fragments.len(),
        "Built phasing problem"
    );
    problem
}

#[derive(Debug, Default)]
struct Builder {
    sites: Sites,
    fragments: Vec<Fragment>,
    /// Where each template sits in `fragments`; insertion order is column
    /// order, so the result does not depend on the hasher.
    index: FxHashMap<FragmentId, usize>,
}

impl Builder {
    fn record(&mut self, site_index: SiteIndex, Vote { fragment: id, allele, cost }: Vote) {
        let next = self.fragments.len();
        let fragment_idx = *self.index.entry(id).or_insert(next);
        if fragment_idx == next {
            self.fragments.push(Fragment { id, supports: SmallVec::new() });
        }
        let Some(fragment) = self.fragments.get_mut(fragment_idx) else {
            return;
        };
        // The overlap dedup leaves one read per template in a column, so a
        // second observation of a site can only come from a `-F` that admits
        // supplementary alignments, or a qname collision; the first one stands.
        // Observations arrive site by site, so it is the entry just pushed.
        if fragment.supports.last().is_some_and(|previous| previous.site == site_index) {
            return;
        }
        fragment.supports.push(AlleleSupport { site: site_index, allele, cost });
    }

    fn finish(mut self) -> PhaseProblem {
        self.fragments.retain(|fragment| fragment.supports.len() > 1);
        PhaseProblem { sites: self.sites, fragments: self.fragments }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        call::{
            pair_hmm::ReadScores,
            phasing::fixtures::{
                ALT_1, ALT_2, backing, column, confirming, graded, het, insertion, mapq, obs,
                scored, with_indels,
            },
            variant_calling::GenotypeTag,
        },
        metrics::AltCall,
        vcf::InCpG,
    };
    use seqair_types::Base::*;

    /// A het whose reference base is one half of a reference CpG.
    fn cpg_het(
        pos: u32,
        reference: Base,
        alt: Base,
        cpg: InCpG,
        observations: &[PhaseObservation],
    ) -> PileupMetrics {
        let mut column = het(pos, reference, alt, observations);
        column.pos_metrics.cpg = cpg;
        column
    }

    fn supports(problem: &PhaseProblem, fragment: u64) -> Vec<(usize, SiteAllele, u8)> {
        let id = FragmentId::new(fragment).expect("non-zero test fragment id");
        problem
            .fragments
            .iter()
            .find(|f| f.id == id)
            .map(|f| f.supports.iter().map(|s| (s.site.get(), s.allele, s.cost.0)).collect())
            .unwrap_or_default()
    }

    #[test]
    fn a_template_spanning_two_hets_becomes_one_fragment() {
        let columns = [
            het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OB)]),
            het(1040, C, G, &[obs(1, G, Strand::OT), obs(2, C, Strand::OB)]),
        ];

        let problem = build_problem(&columns, None);

        assert_eq!(problem.sites.len(), 2);
        assert_eq!(problem.fragments.len(), 2);
        assert_eq!(
            supports(&problem, 1),
            [(0, SiteAllele::First, 26), (1, SiteAllele::Second, 26)],
            "fragment 1 carries reference then alt"
        );
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::Second, 26), (1, SiteAllele::First, 26)],
            "fragment 2 is the opposite haplotype"
        );
    }

    #[test]
    fn a_fragment_reaching_one_site_is_dropped() {
        let columns = [
            het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]),
            het(1040, C, G, &[obs(1, G, Strand::OT)]),
        ];

        let problem = build_problem(&columns, None);

        assert_eq!(problem.fragments.len(), 1, "fragment 2 links nothing");
        assert_eq!(supports(&problem, 2), []);
    }

    #[test]
    fn a_site_no_fragment_reaches_is_still_a_site() {
        let columns = [
            het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]),
            het(1040, C, G, &[obs(3, G, Strand::OT), obs(4, C, Strand::OT)]),
        ];

        let problem = build_problem(&columns, None);

        assert_eq!(problem.sites.len(), 2);
        assert!(problem.fragments.is_empty(), "no template spans both, so nothing links");
    }

    #[test]
    fn only_heterozygous_columns_are_sites() {
        let observations = [obs(1, A, Strand::OT), obs(2, G, Strand::OT)];
        let columns = [
            column(1000, A, &[], Some(GenotypeTag::HomRef), &observations),
            column(
                1040,
                A,
                &[(G, AltCall::RealVariant)],
                Some(GenotypeTag::HomAlt(ALT_1)),
                &observations,
            ),
            column(1080, A, &[(G, AltCall::RealVariant)], None, &observations),
        ];

        assert!(build_problem(&columns, None).sites.is_empty());
    }

    #[test]
    fn a_het_whose_allele_was_not_called_a_real_variant_is_not_a_site() {
        let observations = [obs(1, C, Strand::OB), obs(2, T, Strand::OB)];
        let columns = [column(
            1000,
            C,
            &[(T, AltCall::MethylationEvidenceOnly { for_base: C })],
            Some(GenotypeTag::RefHet(ALT_1)),
            &observations,
        )];

        assert!(build_problem(&columns, None).sites.is_empty());
    }

    /// One read is a coin toss: the site would be phased on it and the
    /// reference reads, and such a site is mostly a false het.
    #[test]
    fn each_allele_needs_two_votes() {
        let site = |observations: &[PhaseObservation]| {
            let columns = [
                column(
                    1000,
                    A,
                    &[(G, AltCall::RealVariant)],
                    Some(GenotypeTag::RefHet(ALT_1)),
                    observations,
                ),
                het(1040, A, G, &[obs(1, A, Strand::OT)]),
            ];
            build_problem(&columns, None).sites.len()
        };
        let reference = [obs(1, A, Strand::OT), obs(2, A, Strand::OT), obs(3, A, Strand::OB)];

        assert_eq!(site(&[reference.as_slice(), &[obs(4, G, Strand::OT)]].concat()), 1);
        assert_eq!(
            site(&[reference.as_slice(), &[obs(4, G, Strand::OT), obs(5, G, Strand::OB)]].concat()),
            2
        );
        assert_eq!(site(&[obs(1, A, Strand::OT), obs(4, G, Strand::OT), obs(5, G, Strand::OB)]), 1);
    }

    /// The floor counts votes after the TAPS gate: at a C/T site in a CpG two
    /// top-strand `T`s are no vote at all.
    #[test]
    fn converted_reads_do_not_back_an_allele() {
        let site = |t_strand| {
            let observations = [
                obs(1, C, Strand::OT),
                obs(2, C, Strand::OB),
                obs(3, T, t_strand),
                obs(4, T, t_strand),
            ];
            let mut cpg = column(
                1000,
                C,
                &[(T, AltCall::RealVariant)],
                Some(GenotypeTag::RefHet(ALT_1)),
                &observations,
            );
            cpg.pos_metrics.cpg = InCpG::C;
            let columns = [cpg, het(1040, A, G, &[obs(1, A, Strand::OT)])];
            build_problem(&columns, None).sites.len()
        };

        assert_eq!(site(Strand::OT), 1);
        assert_eq!(site(Strand::OB), 2);
    }

    #[test]
    fn a_column_without_observations_is_not_a_site() {
        let columns = [het(1000, A, G, &[])];

        assert!(build_problem(&columns, None).sites.is_empty());
    }

    /// At a C/T site in a CpG the top-strand `T` reads may be methylated `C`s.
    #[test]
    fn a_converted_strand_read_is_not_evidence_at_a_c_t_site() {
        let columns = [
            cpg_het(
                1000,
                C,
                T,
                InCpG::C,
                &[
                    obs(1, T, Strand::OT), // may be a methylated C
                    obs(2, T, Strand::OB), // unambiguously the T allele
                    obs(3, C, Strand::OT),
                    obs(4, C, Strand::OB),
                ],
            ),
            het(1040, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT), obs(3, A, Strand::OT)]),
        ];

        let problem = build_problem(&columns, None);

        assert_eq!(supports(&problem, 1), [], "one site left is one site too few");
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::Second, 26), (1, SiteAllele::Second, 26)]
        );
        assert_eq!(supports(&problem, 3), [(0, SiteAllele::First, 26), (1, SiteAllele::First, 26)]);
    }

    /// Mirror image: at a G/A site an `A` is only evidence on an OT read.
    #[test]
    fn a_converted_strand_read_is_not_evidence_at_a_g_a_site() {
        let columns = [
            cpg_het(
                1000,
                G,
                A,
                InCpG::G,
                &[obs(1, A, Strand::OB), obs(2, A, Strand::OT), obs(3, G, Strand::OB)],
            ),
            het(1040, C, G, &[obs(1, C, Strand::OB), obs(2, G, Strand::OB), obs(3, C, Strand::OB)]),
        ];

        let problem = build_problem(&columns, None);

        assert_eq!(supports(&problem, 1), [], "the only OB A was dropped, leaving one site");
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::Second, 26), (1, SiteAllele::Second, 26)]
        );
        assert_eq!(supports(&problem, 3), [(0, SiteAllele::First, 26), (1, SiteAllele::First, 26)]);
    }

    /// A de-novo CpG het is a T/C site with the reference on the ambiguous
    /// side: an OT `T` may be a methylated alternative `C`.
    #[test]
    fn a_denovo_cpg_het_confounds_the_reference_allele() {
        let mut denovo =
            het(1000, T, C, &[obs(1, T, Strand::OT), obs(2, T, Strand::OB), obs(3, C, Strand::OT)]);
        denovo.context.after_1 = Some(G);
        let columns = [
            denovo,
            het(1040, A, G, &[obs(1, A, Strand::OT), obs(2, A, Strand::OT), obs(3, G, Strand::OT)]),
        ];

        let problem = build_problem(&columns, None);

        assert_eq!(supports(&problem, 1), [], "the OT T could be a methylated C");
        assert_eq!(supports(&problem, 2), [(0, SiteAllele::First, 26), (1, SiteAllele::First, 26)]);
        assert_eq!(
            supports(&problem, 3),
            [(0, SiteAllele::Second, 26), (1, SiteAllele::Second, 26)]
        );
    }

    #[test]
    fn an_unknown_strand_read_is_confounded_at_a_c_t_site() {
        let columns = [
            cpg_het(
                1000,
                C,
                T,
                InCpG::C,
                &[obs(1, T, Strand::Unknown), obs(2, C, Strand::Unknown)],
            ),
            het(1040, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]),
        ];

        let problem = build_problem(&columns, None);

        assert_eq!(supports(&problem, 1), [], "the T might be a methylated C");
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::First, 26), (1, SiteAllele::Second, 26)],
            "a C is unambiguous whatever strand it came from"
        );
    }

    /// Outside a CpG the cytosine is almost never methylated, so the top-strand
    /// `T` is the `T` allele.
    #[test]
    fn a_converted_strand_read_is_evidence_outside_a_cpg() {
        let columns = [
            het(1000, C, T, &[obs(1, T, Strand::OT), obs(2, C, Strand::OT)]),
            het(1040, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]),
        ];

        let problem = build_problem(&columns, None);

        assert_eq!(
            supports(&problem, 1),
            [(0, SiteAllele::Second, 26), (1, SiteAllele::First, 26)]
        );
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::First, 26), (1, SiteAllele::Second, 26)]
        );
    }

    /// A T reference becomes a CpG only where the base after it is a `G`.
    #[test]
    fn a_t_c_het_without_a_following_g_is_not_confounded() {
        let columns =
            [het(1000, T, C, &[obs(1, T, Strand::OT)]), het(1040, A, G, &[obs(1, A, Strand::OT)])];

        assert_eq!(
            supports(&build_problem(&columns, None), 1),
            [(0, SiteAllele::First, 26), (1, SiteAllele::First, 26)],
            "the OT T is the reference allele"
        );
    }

    #[test]
    fn a_base_that_is_neither_allele_is_skipped() {
        let columns = [
            het(1000, A, G, &[obs(1, C, Strand::OT), obs(2, A, Strand::OT)]),
            het(1040, A, G, &[obs(1, G, Strand::OT), obs(2, G, Strand::OT)]),
        ];

        let problem = build_problem(&columns, None);

        assert_eq!(supports(&problem, 1), [], "a third base supports neither allele");
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::First, 26), (1, SiteAllele::Second, 26)]
        );
    }

    /// Reads of one template reach a column once each, so a repeat is a
    /// supplementary alignment `-F` let through, and it must not let a
    /// template vote twice.
    #[test]
    fn a_template_observes_a_site_once() {
        let columns = [
            het(1000, A, G, &[graded(1, G, Strand::OT, 22), graded(1, A, Strand::OB, 37)]),
            het(1040, A, G, &[obs(1, A, Strand::OT)]),
        ];

        let problem = build_problem(&columns, None);

        assert_eq!(
            supports(&problem, 1),
            [(0, SiteAllele::Second, 22), (1, SiteAllele::First, 26)],
            "the first observation stands"
        );
    }

    #[test]
    fn cost_is_base_quality_capped() {
        let columns = [
            het(1000, A, G, &[graded(1, G, Strand::OT, 60), graded(2, A, Strand::OT, 25)]),
            het(1040, A, G, &[graded(1, A, Strand::OT, 41), graded(2, G, Strand::OT, 12)]),
        ];

        let problem = build_problem(&columns, None);

        assert_eq!(
            supports(&problem, 1),
            [(0, SiteAllele::Second, 26), (1, SiteAllele::First, 26)]
        );
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::First, 25), (1, SiteAllele::Second, 12)]
        );
    }

    /// A compound het has no reference allele, so observations of the
    /// reference base are noise.
    #[test]
    fn a_compound_het_phases_its_two_alts() {
        let alts = [(T, AltCall::RealVariant), (G, AltCall::RealVariant)];
        let columns = [
            column(
                1000,
                A,
                &alts,
                Some(GenotypeTag::AltHet(ALT_1, ALT_2)),
                &[
                    [obs(1, T, Strand::OB), obs(2, G, Strand::OT), obs(3, A, Strand::OT)]
                        .as_slice(),
                    &backing(1000, [T, G]),
                ]
                .concat(),
            ),
            het(1040, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT), obs(3, G, Strand::OT)]),
        ];

        let problem = build_problem(&columns, None);

        let site = &problem.sites[SiteIndex(0)];
        assert_eq!(site.alleles, [GtAllele::Alt(ALT_1), GtAllele::Alt(ALT_2)]);
        assert_eq!(supports(&problem, 1), [(0, SiteAllele::First, 26), (1, SiteAllele::First, 26)]);
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::Second, 26), (1, SiteAllele::Second, 26)]
        );
        assert_eq!(supports(&problem, 3), [], "the reference base is neither allele");
    }

    fn guards_of(columns: &[PileupMetrics]) -> Vec<Option<Guard>> {
        build_problem(columns, None).sites.iter().map(|site| site.guard).collect()
    }

    fn lone_het(pos: u32) -> PileupMetrics {
        het(pos, A, G, &[obs(u64::from(pos), A, Strand::OT)])
    }

    #[test]
    fn an_snv_allele_mapped_below_mapq_50_is_guarded() {
        let with_mapq = |reference: u8, alt: u8| {
            let mut column = lone_het(1000);
            column.ref_metrics.mapq = mapq(reference);
            if let Some(alt_allele) = column.alts.first_mut() {
                alt_allele.metrics.mapq = mapq(alt);
            }
            guards_of(&[column])
        };

        assert_eq!(with_mapq(50, 50), [None]);
        assert_eq!(with_mapq(49, 60), [Some(Guard::AlleleMapq)]);
        assert_eq!(with_mapq(60, 49), [Some(Guard::AlleleMapq)]);
    }

    #[test]
    fn a_site_over_twice_the_median_depth_is_guarded() {
        let with_depths = |depths: [u32; 3]| {
            let columns = depths.into_iter().zip([1000, 2000, 3000]).map(|(depth, pos)| {
                let mut column = lone_het(pos);
                column.pos_metrics.depth = depth;
                column
            });
            guards_of(&columns.collect::<Vec<_>>())
        };

        assert_eq!(with_depths([10, 10, 20]), [None, None, None]);
        assert_eq!(with_depths([10, 10, 21]), [None, None, Some(Guard::Depth)]);
    }

    #[test]
    fn the_depth_guard_scales_the_median_to_the_run() {
        let with_scale = |scale| {
            let columns: Vec<_> = [(10, 1000), (10, 2000), (25, 3000)]
                .into_iter()
                .map(|(depth, pos)| {
                    let mut column = lone_het(pos);
                    column.pos_metrics.depth = depth;
                    column
                })
                .collect();
            build_problem(&columns, scale).sites.iter().map(|site| site.guard).collect::<Vec<_>>()
        };

        assert_eq!(with_scale(None), [None, None, Some(Guard::Depth)]);
        assert_eq!(with_scale(Some(DepthScale(1.0))), [None, None, Some(Guard::Depth)]);
        // The segment holds three times the run's reads per base: its median
        // is a repeat's, and the run expects a third of it.
        assert_eq!(with_scale(Some(DepthScale(0.33))), [Some(Guard::Depth); 3]);
        assert_eq!(with_scale(Some(DepthScale(1.5))), [None; 3]);
    }

    #[test]
    fn a_site_with_five_others_within_100_bases_is_guarded() {
        let at = |positions: &[u32]| {
            guards_of(&positions.iter().map(|&pos| lone_het(pos)).collect::<Vec<_>>())
        };

        assert_eq!(at(&[1000, 1020, 1040, 1060, 1080]), [None; 5], "four others");
        assert_eq!(at(&[1000, 1020, 1040, 1060, 1080, 1100]), [Some(Guard::HetDensity); 6]);
        assert_eq!(
            at(&[1000, 1020, 1040, 1060, 1080, 1101]),
            [
                None,
                Some(Guard::HetDensity),
                Some(Guard::HetDensity),
                Some(Guard::HetDensity),
                Some(Guard::HetDensity),
                None
            ],
            "1000 and 1101 are 101 bases apart, so each has four others"
        );
    }

    #[test]
    fn an_indel_site_is_not_guarded_on_mapq() {
        let mut indel = ref_het_indel(1040, confirming([0, 1]));
        indel.ref_metrics.mapq = mapq(0);
        assert_eq!(guards_of(&[indel]), [None]);
    }

    /// The solver writes its answer back through this index.
    #[test]
    fn a_site_points_back_at_its_column() {
        let observations = [obs(1, A, Strand::OT), obs(2, G, Strand::OT)];
        let columns = [
            column(1000, A, &[], Some(GenotypeTag::HomRef), &observations),
            het(1040, A, G, &observations),
            column(1080, A, &[], Some(GenotypeTag::HomRef), &observations),
            het(1120, A, G, &observations),
        ];

        let problem = build_problem(&columns, None);

        let located: Vec<_> = problem.sites.iter().map(|s| (s.pileup, s.pos.as_u32())).collect();
        assert_eq!(located, [(1, 1040), (3, 1120)]);
    }

    const NO: f32 = f32::NEG_INFINITY;

    fn sites_of(problem: &PhaseProblem) -> Vec<(u32, SiteKind, [GtAllele; 2])> {
        problem.sites.iter().map(|s| (s.pos.as_u32(), s.kind, s.alleles)).collect()
    }

    /// An SNV het every fragment reaches, so a fragment's support list shows
    /// whether it also voted at the next site.
    fn anchor(fragments: std::ops::RangeInclusive<u64>) -> PileupMetrics {
        let observations: Vec<_> = fragments.map(|f| obs(f, A, Strand::OT)).collect();
        het(1000, A, G, &observations)
    }

    /// A `0/1` insertion at `pos` the pair-HMM calls from `reads`.
    fn ref_het_indel(pos: u32, reads: impl IntoIterator<Item = ReadScores>) -> PileupMetrics {
        let reads: Vec<_> = reads.into_iter().collect();
        with_indels(
            column(pos, A, &[], Some(GenotypeTag::HomRef), &[]),
            &[insertion(&[T, T])],
            &reads,
        )
    }

    /// The vote is the sign of the log10 ratio, its cost the ratio's Phred,
    /// rounded and capped; a read that cannot tell the two apart, one below
    /// the phasing MAPQ and one without a fragment cast none.
    #[test]
    fn an_indel_het_is_a_site_its_scores_vote_on() {
        let mut low_mapq = scored(5, [-1.0, -4.0, NO]);
        low_mapq.mapq = 19;
        let mut unnamed = scored(9, [-1.0, -4.0, NO]);
        unnamed.fragment = None;
        let reads = [
            scored(1, [-1.0, -3.0, NO]),
            scored(2, [-6.0, -1.0, NO]),
            scored(3, [NO, -1.0, NO]),
            low_mapq,
            scored(6, [-1.0, -1.04, NO]),
            scored(7, [-1.27, -1.0, NO]),
            scored(8, [-1.0, -1.0, NO]),
            unnamed,
        ];
        let columns =
            [anchor(1..=8), ref_het_indel(1040, reads.into_iter().chain(confirming([0, 1])))];

        let problem = build_problem(&columns, None);

        assert_eq!(
            sites_of(&problem),
            [
                (1000, SiteKind::Snv, [GtAllele::Ref, GtAllele::Alt(ALT_1)]),
                (1040, SiteKind::Indel, [GtAllele::Ref, GtAllele::Alt(ALT_1)]),
            ]
        );
        let at_indel = |fragment| -> Vec<_> {
            supports(&problem, fragment).into_iter().filter(|s| s.0 == 1).collect()
        };
        assert_eq!(at_indel(1), [(1, SiteAllele::First, 20)]);
        assert_eq!(at_indel(2), [(1, SiteAllele::Second, 26)], "capped");
        assert_eq!(
            at_indel(3),
            [(1, SiteAllele::Second, 26)],
            "an impossible haplotype is the cap"
        );
        assert_eq!(at_indel(5), [], "below the phasing MAPQ");
        assert_eq!(at_indel(6), [], "rounds to nothing");
        assert_eq!(at_indel(7), [(1, SiteAllele::Second, 3)]);
        assert_eq!(at_indel(8), [], "no preference");
        assert_eq!(problem.fragments.len(), 4, "fragments 1, 2, 3 and 7 link the two sites");
    }

    /// A het event of a cluster votes on the called pair's two haplotypes as
    /// they are, the cluster's other events included: the carrier here also
    /// carries the deletion at 1050, and the haplotype of the insertion alone
    /// explains no read, so scoring against it would leave the carriers
    /// without a vote.
    #[test]
    fn a_cluster_event_votes_on_the_haplotypes_of_the_called_pair() {
        use crate::call::{
            pair_hmm::{Cluster, Event, HaplotypeScores, MAX_HAPLOTYPES},
            pileup::indels::{IndelAllele, IndelCounts, IndelData, Scored},
        };
        use std::sync::Arc;
        let read = |fragment, hap: usize| {
            let mut log10 = [-4.0; MAX_HAPLOTYPES];
            if let Some(likeliest) = log10.get_mut(hap) {
                *likeliest = -1.0;
            }
            ReadScores { fragment: FragmentId::new(fragment), mapq: 60, log10 }
        };
        // Haplotype 0 is the reference, 3 carries both events.
        let reads: Vec<ReadScores> =
            (1..=6).map(|f| read(f, 0)).chain((7..=12).map(|f| read(f, 3))).collect();
        let events = [
            Event { anchor: 1040, allele: insertion(&[T, T]) },
            Event { anchor: 1050, allele: IndelAllele::Deletion([C, A].into_iter().collect()) },
        ];
        let cluster = Cluster::new(events.into_iter().collect()).expect("a cluster");
        let scores = HaplotypeScores::new(cluster, reads.into()).expect("a verdict");
        let mut indel = column(1040, A, &[], Some(GenotypeTag::HomRef), &[]);
        let scored = Some(Scored::new(Arc::new(scores), 1040));
        indel.indel_data = Some(Box::new(IndelData { counts: IndelCounts::default(), scored }));

        let problem = build_problem(&[anchor(1..=12), indel], None);

        assert_eq!(
            sites_of(&problem),
            [
                (1000, SiteKind::Snv, [GtAllele::Ref, GtAllele::Alt(ALT_1)]),
                (1040, SiteKind::Indel, [GtAllele::Ref, GtAllele::Alt(ALT_1)]),
            ]
        );
        let at_indel = |fragment| -> Vec<_> {
            supports(&problem, fragment).into_iter().filter(|s| s.0 == 1).collect()
        };
        assert_eq!(at_indel(1), [(1, SiteAllele::First, 26)]);
        assert_eq!(at_indel(7), [(1, SiteAllele::Second, 26)]);
    }

    /// The pair-HMM has to be sure of the column's het, or the site's votes
    /// are noise that can bridge two blocks the wrong way round; a
    /// homozygote is never a site.
    #[test]
    fn an_indel_het_the_pair_hmm_doubts_is_not_a_site() {
        let sites = |reads: Vec<ReadScores>| {
            build_problem(&[anchor(1..=1), ref_het_indel(1040, reads)], None).sites.len()
        };

        assert_eq!(sites(confirming([0, 1])), 2);
        assert_eq!(sites(confirming([1, 1]).repeat(4)), 1, "a confident homozygote");
        assert_eq!(
            sites(vec![scored(1, [-1.0, -4.0, NO]), scored(2, [-4.0, -1.0, NO])]),
            1,
            "a het, but at about Q21"
        );
    }

    /// However sure the pair-HMM is, each allele needs votes of its own: here
    /// four reads the reference cannot explain, against five of the
    /// reference, make a het at the top quality, and still only five make a
    /// site. A read that cannot vote, without a fragment or below the phasing
    /// MAPQ, does not count.
    #[test]
    fn an_indel_allele_needs_five_votes_of_its_own() {
        let sites = |alt_reads: Vec<ReadScores>| {
            let reads = (1..=5).map(|fragment| scored(fragment, [-1.0, -9.0, NO])).chain(alt_reads);
            build_problem(&[anchor(1..=1), ref_het_indel(1040, reads)], None).sites.len()
        };
        let alt = |fragment| scored(fragment, [NO, -1.0, NO]);
        let mut low_mapq = alt(15);
        low_mapq.mapq = 19;
        let unnamed = ReadScores { fragment: None, ..alt(16) };
        assert_eq!(sites((11..=14).map(alt).collect()), 1);
        assert_eq!(sites((11..=15).map(alt).collect()), 2);
        assert_eq!(sites((11..=14).map(alt).chain([low_mapq, unnamed]).collect()), 1);
    }

    /// The SNV record and the indel record of one column are separate sites,
    /// SNV first, and one fragment can vote at both.
    #[test]
    fn an_snv_and_an_indel_at_one_column_are_two_sites() {
        let both = with_indels(
            het(1000, A, G, &[obs(1, G, Strand::OT), obs(2, A, Strand::OT)]),
            &[insertion(&[T, T])],
            &[
                [scored(1, [-4.0, -1.0, NO]), scored(2, [-1.0, -4.0, NO])].as_slice(),
                &confirming([0, 1]),
            ]
            .concat(),
        );
        let columns = [both, het(1040, C, T, &[obs(1, T, Strand::OB)])];

        let problem = build_problem(&columns, None);

        assert_eq!(
            sites_of(&problem).into_iter().map(|(pos, kind, _)| (pos, kind)).collect::<Vec<_>>(),
            [(1000, SiteKind::Snv), (1000, SiteKind::Indel), (1040, SiteKind::Snv)]
        );
        assert!(
            problem.sites.iter().all(|site| site.pileup == usize::from(site.pos.as_u32() != 1000))
        );
        assert_eq!(
            supports(&problem, 1),
            [(0, SiteAllele::Second, 26), (1, SiteAllele::Second, 26), (2, SiteAllele::Second, 26)]
        );
        assert_eq!(supports(&problem, 2), [(0, SiteAllele::First, 26), (1, SiteAllele::First, 26)]);
    }

    /// A compound het's alleles are the record's `1/2`, each fragment voting
    /// for the haplotype that carries its allele.
    #[test]
    fn a_compound_het_indel_votes_on_its_record_alleles() {
        let (short, long) = (insertion(&[T]), insertion(&[T, T, T]));
        let indel = with_indels(
            column(1040, A, &[], Some(GenotypeTag::HomRef), &[]),
            &[short, long],
            &[
                [scored(1, [-9.0, -1.0, -9.0]), scored(2, [-9.0, -9.0, -1.0])].as_slice(),
                &confirming([1, 2]),
            ]
            .concat(),
        );

        let problem = build_problem(&[anchor(1..=2), indel], None);

        assert_eq!(
            sites_of(&problem).get(1),
            Some(&(1040, SiteKind::Indel, [GtAllele::Alt(ALT_1), GtAllele::Alt(ALT_2)]))
        );
        assert_eq!(supports(&problem, 1), [(0, SiteAllele::First, 26), (1, SiteAllele::First, 26)]);
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::First, 26), (1, SiteAllele::Second, 26)]
        );
    }
}
