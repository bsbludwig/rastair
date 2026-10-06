//! Turning a segment's columns into a phasing problem.
//!
//! Two lists come out of a segment: the *sites*, heterozygous columns whose two
//! alleles have to be split across haplotypes, and the *fragments*, templates
//! whose alleles at several sites share a chromosome. A fragment carries one
//! bit per site saying which allele it saw; at a C/T or G/A heterozygote whose
//! cytosine is in a CpG, observations TAPS conversion alone explains are
//! dropped, or methylation would read as the other allele.

use super::{Cost, FragmentId, GtAllele, PhaseObservation, observations_of};
use crate::{
    call::variant_calling::GenotypeTag,
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

/// One allele of a heterozygous genotype.
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

/// A heterozygous column a solver can phase.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PhaseSite {
    /// Where the column sits in the slice [`build_problem`] was given.
    pub pileup: usize,
    pub pos: Pos0,
    /// The genotype's two alleles, ascending.
    pub alleles: [HetAllele; 2],
    /// The CpG cytosine one allele is, when conversion makes it read as the
    /// other: a C/T site's top strand or a G/A site's bottom strand.
    pub methylable: Option<CpgSide>,
}

impl PhaseSite {
    fn of(pileup: &PileupMetrics, pileup_idx: usize) -> Option<Self> {
        let alleles = het_alleles(pileup)?;
        let bases = alleles.map(|allele| allele.base);
        // Outside a CpG, cytosine methylation is too rare in mammals to give
        // up a strand's reads for.
        let methylable = [CpgSide::C, CpgSide::G].into_iter().find(|&side| {
            bases.contains(&side.unmod_base())
                && bases.contains(&side.mod_base())
                && cpg_origin(pileup, side).is_some()
        });
        Some(Self { pileup: pileup_idx, pos: Pos0::new(pileup.pos)?, alleles, methylable })
    }

    #[must_use]
    pub const fn allele(&self, allele: SiteAllele) -> HetAllele {
        let [first, second] = self.alleles;
        match allele {
            SiteAllele::First => first,
            SiteAllele::Second => second,
        }
    }

    /// Which allele an observed base is, or `None` when it is neither.
    fn allele_of(&self, base: Base) -> Option<SiteAllele> {
        match self.alleles.map(|allele| allele.base == base) {
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
pub fn build_problem(pileups: &[PileupMetrics]) -> PhaseProblem {
    let mut builder = Builder::default();

    for (pileup_idx, pileup) in pileups.iter().enumerate() {
        let Some(observations) = observations_of(pileup) else {
            continue;
        };
        let Some(site) = PhaseSite::of(pileup, pileup_idx) else {
            continue;
        };
        let Some(index) = builder.sites.push(site) else {
            warn!(
                sites = builder.sites.len(),
                "Too many heterozygous columns in one segment to phase; the rest stay unphased"
            );
            break;
        };
        for observation in observations {
            builder.observe(index, observation);
        }
    }

    let problem = builder.finish();
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
    fn observe(&mut self, site_index: SiteIndex, observation: &PhaseObservation) {
        let site = &self.sites[site_index];
        let Some(allele) = site.allele_of(observation.base) else {
            return;
        };
        if site.taps_confounds(observation.base, observation.strand) {
            return;
        }
        let cost = Cost::of(observation.qual);

        let next = self.fragments.len();
        let fragment_idx = *self.index.entry(observation.fragment).or_insert(next);
        if fragment_idx == next {
            self.fragments.push(Fragment { id: observation.fragment, supports: SmallVec::new() });
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
            phasing::fixtures::{ALT_1, ALT_2, column, graded, het, obs},
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

        let problem = build_problem(&columns);

        assert_eq!(problem.sites.len(), 2);
        assert_eq!(problem.fragments.len(), 2);
        assert_eq!(
            supports(&problem, 1),
            [(0, SiteAllele::First, 30), (1, SiteAllele::Second, 30)],
            "fragment 1 carries reference then alt"
        );
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::Second, 30), (1, SiteAllele::First, 30)],
            "fragment 2 is the opposite haplotype"
        );
    }

    #[test]
    fn a_fragment_reaching_one_site_is_dropped() {
        let columns = [
            het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]),
            het(1040, C, G, &[obs(1, G, Strand::OT)]),
        ];

        let problem = build_problem(&columns);

        assert_eq!(problem.fragments.len(), 1, "fragment 2 links nothing");
        assert_eq!(supports(&problem, 2), []);
    }

    #[test]
    fn a_site_no_fragment_reaches_is_still_a_site() {
        let columns = [
            het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]),
            het(1040, C, G, &[obs(3, G, Strand::OT), obs(4, C, Strand::OT)]),
        ];

        let problem = build_problem(&columns);

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

        assert!(build_problem(&columns).sites.is_empty());
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

        assert!(build_problem(&columns).sites.is_empty());
    }

    #[test]
    fn a_column_without_observations_is_not_a_site() {
        let columns = [het(1000, A, G, &[])];

        assert!(build_problem(&columns).sites.is_empty());
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

        let problem = build_problem(&columns);

        assert_eq!(supports(&problem, 1), [], "one site left is one site too few");
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::Second, 30), (1, SiteAllele::Second, 30)]
        );
        assert_eq!(supports(&problem, 3), [(0, SiteAllele::First, 30), (1, SiteAllele::First, 30)]);
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

        let problem = build_problem(&columns);

        assert_eq!(supports(&problem, 1), [], "the only OB A was dropped, leaving one site");
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::Second, 30), (1, SiteAllele::Second, 30)]
        );
        assert_eq!(supports(&problem, 3), [(0, SiteAllele::First, 30), (1, SiteAllele::First, 30)]);
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

        let problem = build_problem(&columns);

        assert_eq!(supports(&problem, 1), [], "the OT T could be a methylated C");
        assert_eq!(supports(&problem, 2), [(0, SiteAllele::First, 30), (1, SiteAllele::First, 30)]);
        assert_eq!(
            supports(&problem, 3),
            [(0, SiteAllele::Second, 30), (1, SiteAllele::Second, 30)]
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

        let problem = build_problem(&columns);

        assert_eq!(supports(&problem, 1), [], "the T might be a methylated C");
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::First, 30), (1, SiteAllele::Second, 30)],
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

        let problem = build_problem(&columns);

        assert_eq!(
            supports(&problem, 1),
            [(0, SiteAllele::Second, 30), (1, SiteAllele::First, 30)]
        );
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::First, 30), (1, SiteAllele::Second, 30)]
        );
    }

    /// A T reference becomes a CpG only where the base after it is a `G`.
    #[test]
    fn a_t_c_het_without_a_following_g_is_not_confounded() {
        let columns =
            [het(1000, T, C, &[obs(1, T, Strand::OT)]), het(1040, A, G, &[obs(1, A, Strand::OT)])];

        assert_eq!(
            supports(&build_problem(&columns), 1),
            [(0, SiteAllele::First, 30), (1, SiteAllele::First, 30)],
            "the OT T is the reference allele"
        );
    }

    #[test]
    fn a_base_that_is_neither_allele_is_skipped() {
        let columns = [
            het(1000, A, G, &[obs(1, C, Strand::OT), obs(2, A, Strand::OT)]),
            het(1040, A, G, &[obs(1, G, Strand::OT), obs(2, G, Strand::OT)]),
        ];

        let problem = build_problem(&columns);

        assert_eq!(supports(&problem, 1), [], "a third base supports neither allele");
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::First, 30), (1, SiteAllele::Second, 30)]
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

        let problem = build_problem(&columns);

        assert_eq!(
            supports(&problem, 1),
            [(0, SiteAllele::Second, 22), (1, SiteAllele::First, 30)],
            "the first observation stands"
        );
    }

    #[test]
    fn cost_is_base_quality_capped() {
        let columns = [
            het(1000, A, G, &[graded(1, G, Strand::OT, 60), graded(2, A, Strand::OT, 25)]),
            het(1040, A, G, &[graded(1, A, Strand::OT, 41), graded(2, G, Strand::OT, 12)]),
        ];

        let problem = build_problem(&columns);

        assert_eq!(
            supports(&problem, 1),
            [(0, SiteAllele::Second, 40), (1, SiteAllele::First, 40)]
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
                &[obs(1, T, Strand::OB), obs(2, G, Strand::OT), obs(3, A, Strand::OT)],
            ),
            het(1040, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT), obs(3, G, Strand::OT)]),
        ];

        let problem = build_problem(&columns);

        let site = &problem.sites[SiteIndex(0)];
        assert_eq!(
            site.alleles,
            [
                HetAllele { allele: GtAllele::Alt(ALT_1), base: T },
                HetAllele { allele: GtAllele::Alt(ALT_2), base: G }
            ]
        );
        assert_eq!(supports(&problem, 1), [(0, SiteAllele::First, 30), (1, SiteAllele::First, 30)]);
        assert_eq!(
            supports(&problem, 2),
            [(0, SiteAllele::Second, 30), (1, SiteAllele::Second, 30)]
        );
        assert_eq!(supports(&problem, 3), [], "the reference base is neither allele");
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

        let problem = build_problem(&columns);

        let located: Vec<_> = problem.sites.iter().map(|s| (s.pileup, s.pos.as_u32())).collect();
        assert_eq!(located, [(1, 1040), (3, 1120)]);
    }
}
