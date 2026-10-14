//! Fragment likelihood phasing of called SNP genotypes. Full read evidence is
//! captured before seqair recycles its store; chromosome assignments are
//! scored jointly with shared CpG rates and local `PairHMM` realignment.

use crate::metrics::{PileupMetrics, methylation::CpgSide};
#[cfg(test)]
use crate::utils::logging::ThisIsABug as _;
#[cfg(any(feature = "experimental-seqair", test))]
use color_eyre::eyre::OptionExt as _;
use seqair_types::{Base, BaseQuality, Pos0, Pos1, Strand};
use std::num::{NonZeroU8, NonZeroU64};
#[cfg(any(feature = "experimental-seqair", test))]
use tracing::debug;
#[cfg(test)]
use tracing::warn;

mod mask;
#[cfg(test)]
mod problem;
#[cfg(test)]
mod solver;
pub use mask::PhaseMask;

#[cfg(test)]
use problem::{SiteAllele, build_problem};
#[cfg(test)]
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
#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Cost(u8);

#[cfg(test)]
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

#[cfg(test)]
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

#[cfg(test)]
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
#[cfg(test)]
fn legacy_phase_segment(pileups: &mut [PileupMetrics], is_core: impl Fn(u64) -> bool) {
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

#[cfg(feature = "experimental-seqair")]
pub fn phase_segment(
    pileups: &mut [PileupMetrics],
    segment: &crate::sequence::Segment,
    fragments: &[crate::call::haplotypes::Fragment],
    mask: &PhaseMask,
) -> color_eyre::eyre::Result<()> {
    use color_eyre::eyre::OptionExt as _;
    let end = segment
        .start
        .checked_add(u64::try_from(segment.sequence.len())?)
        .ok_or_eyre("Phasing reference coordinate overflow")?;
    let mut calls = Vec::new();
    // Separate inference domains also prevent reads, CpG partners and realignment
    // flanks from importing masked evidence into otherwise unmasked calls.
    for region in mask.unmasked(segment.contig.as_str(), segment.start..end) {
        if pileups
            .iter()
            .filter(|p| region.contains(&u64::from(p.pos)))
            .filter_map(|p| CalledSite::of(p, 0))
            .take(2)
            .count()
            < 2
        {
            continue;
        }
        calls.extend(phase_unmasked(pileups, segment, fragments, region)?);
    }
    for pileup in pileups.iter_mut() {
        pileup.phase = None;
    }
    for (index, call) in calls {
        pileups.get_mut(index).ok_or_eyre("Phase call outside segment")?.phase =
            Some(Box::new(ColumnPhase::Phased(call)));
    }
    Ok(())
}

#[cfg(feature = "experimental-seqair")]
fn phase_unmasked(
    pileups: &[PileupMetrics],
    segment: &crate::sequence::Segment,
    fragments: &[crate::call::haplotypes::Fragment],
    region: std::ops::Range<u64>,
) -> color_eyre::eyre::Result<Vec<(usize, PhaseCall)>> {
    use crate::call::haplotypes::{
        Reference, Site,
        phase::{State, parameters, score},
    };
    use color_eyre::eyre::{OptionExt as _, ensure};
    use std::collections::BTreeMap;
    let mut reference = Reference {
        start: u32::try_from(region.start)?,
        bases: segment
            .sequence
            .get(
                usize::try_from(region.start - segment.start)?
                    ..usize::try_from(region.end - segment.start)?,
            )
            .ok_or_eyre("Unmasked interval outside phasing segment")?
            .iter()
            .copied()
            .map(Base::from)
            .collect(),
    };
    let sites: Vec<_> = pileups
        .iter()
        .enumerate()
        .filter(|(_, pileup)| region.contains(&u64::from(pileup.pos)))
        .filter_map(|(index, pileup)| CalledSite::of(pileup, index))
        .collect();
    // Homozygous SNPs belong to both candidate chromosomes, including CpG gains and losses.
    for pileup in pileups.iter() {
        if let Some(crate::call::variant_calling::EstimatedGenotype {
            genotype: crate::call::variant_calling::GenotypeTag::HomAlt(alt),
            ..
        }) = pileup.pos_metrics.extended.genotype
            && let Some(base) = pileup
                .alts
                .get(usize::from(alt.get()) - 1)
                .filter(|alt| alt.call == crate::metrics::AltCall::RealVariant)
                .map(|alt| alt.base)
            && let Some(target) = pileup
                .pos
                .checked_sub(reference.start)
                .and_then(|offset| reference.bases.get_mut(offset as usize))
        {
            *target = base;
        }
    }
    let linkage = crate::call::haplotypes::phase::links(
        &sites
            .iter()
            .map(|site| {
                (site.pos.as_u32(), site.alleles.map(|allele| allele.base), site.methylable)
            })
            .collect::<Vec<_>>(),
        fragments,
    )?;
    let links: BTreeMap<_, _> = linkage
        .edges
        .into_iter()
        .filter(|((a, b), _)| {
            [a, b].into_iter().all(|&index| {
                linkage
                    .counts
                    .get(index)
                    .is_some_and(|counts| counts.iter().all(|&count| count >= 2))
            })
        })
        .collect();
    let mut components: Vec<usize> = (0..sites.len()).collect();
    for (&(a, b), &count) in &links {
        if count >= 2 {
            unite(&mut components, a, b)?;
        }
    }
    let mut dependencies = std::collections::BTreeSet::new();
    for (a, site) in sites.iter().enumerate() {
        if !site.alleles.iter().any(|allele| allele.base == Base::C) {
            continue;
        }
        let Some(partner) = site.pos.as_u32().checked_add(1) else {
            continue;
        };
        if let Ok(b) = sites.binary_search_by_key(&partner, |site| site.pos.as_u32())
            && sites
                .get(b)
                .is_some_and(|site| site.alleles.iter().any(|allele| allele.base == Base::G))
        {
            dependencies.insert((a, b));
            unite(&mut components, a, b)?;
        }
    }
    let mut groups = BTreeMap::<usize, Vec<usize>>::new();
    for index in 0..sites.len() {
        groups.entry(root(&components, index)?).or_default().push(index);
    }
    let params = parameters()?;
    let mut calls = Vec::new();
    let mut evaluated = 0usize;
    let mut pruned = 0usize;
    for group in groups.into_values().filter(|group| group.len() >= 2) {
        let mut remaining = group.as_slice();
        while let Some(&anchor) = remaining.first() {
            let anchor_pos = sites.get(anchor).ok_or_eyre("Missing phase anchor")?.pos.as_u32();
            let count = remaining
                .iter()
                .take(32)
                .take_while(|&&index| {
                    sites
                        .get(index)
                        .is_some_and(|site| site.pos.as_u32().saturating_sub(anchor_pos) <= 2000)
                })
                .count();
            let candidates = remaining.get(..count).ok_or_eyre("Invalid phase block bounds")?;
            let mut ordered = Vec::with_capacity(count);
            let mut pending = candidates.to_vec();
            while !pending.is_empty() {
                let next = pending
                    .iter()
                    .position(|&index| {
                        ordered.is_empty()
                            || ordered.iter().any(|&active| {
                                let pair = (index.min(active), index.max(active));
                                links.contains_key(&pair) || dependencies.contains(&pair)
                            })
                    })
                    .unwrap_or(0);
                ordered.push(pending.remove(next));
            }
            let block = ordered.as_slice();
            remaining = remaining.get(count..).ok_or_eyre("Invalid phase remainder")?;
            if block.len() < 2 {
                continue;
            }
            let domains = block
                .iter()
                .map(|&index| {
                    let site = sites.get(index).ok_or_eyre("Missing phase site")?;
                    Ok(Site {
                        pos: site.pos.as_u32(),
                        reference: reference
                            .base_at(site.pos.as_u32())
                            .ok_or_eyre("Phase site outside reference")?,
                        alleles: site.alleles.iter().map(|allele| allele.base).collect(),
                    })
                })
                .collect::<color_eyre::eyre::Result<Vec<_>>>()?;
            let included: Vec<_> = fragments
                .iter()
                .filter(|fragment| {
                    domains.iter().any(|site| {
                        fragment.reads.iter().any(|read| {
                            read.observations.binary_search_by_key(&site.pos, |obs| obs.pos).is_ok()
                        })
                    })
                })
                .collect();
            let first = domains
                .first()
                .and_then(|site| site.alleles.first())
                .copied()
                .ok_or_eyre("Missing anchor allele")?;
            let mut states = vec![State { first: vec![first] }];
            let mut held = Vec::new();
            for (step, domain) in domains.iter().enumerate().skip(1) {
                states = states
                    .into_iter()
                    .chain(held.drain(..))
                    .flat_map(|state| {
                        domain.alleles.iter().map(move |&base| {
                            let mut next = state.clone();
                            next.first.push(base);
                            next
                        })
                    })
                    .collect();
                if states.len() > 48 {
                    let prefix = domains.get(..=step).ok_or_eyre("Missing phase prefix")?;
                    let prefix_reads: Vec<_> = included
                        .iter()
                        .copied()
                        .filter(|fragment| {
                            prefix.iter().any(|site| {
                                fragment.reads.iter().any(|read| {
                                    read.observations
                                        .binary_search_by_key(&site.pos, |obs| obs.pos)
                                        .is_ok()
                                })
                            })
                        })
                        .collect();
                    let posterior = score(&reference, prefix, &states, &prefix_reads, &params)?;
                    let mut ranked: Vec<_> =
                        states.into_iter().zip(posterior.probabilities).collect();
                    ranked.sort_by(|a, b| b.1.total_cmp(&a.1));
                    // Keep a competing orientation at every site even when a beam would discard it.
                    let mut protected = std::collections::BTreeSet::new();
                    protected.insert(0usize);
                    for offset in 1..=step {
                        for base in prefix
                            .get(offset)
                            .ok_or_eyre("Missing protected phase site")?
                            .alleles
                            .iter()
                        {
                            if let Some(index) = ranked
                                .iter()
                                .position(|(state, _)| state.first.get(offset) == Some(base))
                            {
                                protected.insert(index);
                            }
                        }
                    }
                    ensure!(
                        protected.len() <= 64,
                        "Phase protection exceeded the chromosome-pair budget"
                    );
                    let keep = 32usize.max(protected.len());
                    let mut next = Vec::new();
                    for (index, (state, _)) in ranked.into_iter().enumerate() {
                        if protected.contains(&index) || next.len() < keep {
                            next.push(state);
                        } else if held.len() < 16 {
                            held.push(state);
                        } else {
                            pruned += 1;
                        }
                    }
                    states = next;
                }
            }
            states.extend(held);
            let posterior = score(&reference, &domains, &states, &included, &params)?;
            evaluated += 1;
            let best = posterior
                .probabilities
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .map(|(index, _)| index)
                .ok_or_eyre("No phase posterior")?;
            let best = states.get(best).ok_or_eyre("Missing best phase state")?;
            let mut confident: Vec<usize> = (0..block.len()).collect();
            for (a, &index_a) in block.iter().enumerate() {
                for (b, &index_b) in block.iter().enumerate().skip(a + 1) {
                    if links
                        .get(&(index_a.min(index_b), index_a.max(index_b)))
                        .copied()
                        .unwrap_or(0)
                        < 2
                    {
                        continue;
                    }
                    let low = domains
                        .get(a)
                        .and_then(|site| site.alleles.first())
                        .ok_or_eyre("Missing first allele")?;
                    let high = domains
                        .get(b)
                        .and_then(|site| site.alleles.first())
                        .ok_or_eyre("Missing second allele")?;
                    let cis = (best.first.get(a) == Some(low)) == (best.first.get(b) == Some(high));
                    let probability: f64 = states
                        .iter()
                        .zip(&posterior.probabilities)
                        .filter(|(state, _)| {
                            ((state.first.get(a) == Some(low))
                                == (state.first.get(b) == Some(high)))
                                == cis
                        })
                        .map(|(_, probability)| probability)
                        .sum();
                    if probability >= 0.99 {
                        unite(&mut confident, a, b)?;
                    }
                }
            }
            let mut emitted = BTreeMap::<usize, Vec<usize>>::new();
            for (local, &index) in block.iter().enumerate() {
                let site = sites.get(index).ok_or_eyre("Missing emitted phase site")?;
                if segment.is_core(site.pos.as_u64()) {
                    emitted.entry(root(&confident, local)?).or_default().push(local);
                }
            }
            for mut members in emitted.into_values().filter(|members| members.len() >= 2) {
                members.sort_by_key(|&local| block.get(local).copied());
                let anchor = *members.first().ok_or_eyre("No emitted phase anchor")?;
                let anchor_site = sites
                    .get(*block.get(anchor).ok_or_eyre("Missing block anchor")?)
                    .ok_or_eyre("Missing emitted anchor")?;
                let set = PhaseSet::at(anchor_site.pos)
                    .ok_or_eyre("Phase anchor has no 1-based coordinate")?;
                let flip = best.first.get(anchor)
                    != Some(
                        &anchor_site.alleles.first().ok_or_eyre("Missing anchor genotype")?.base,
                    );
                for local in members {
                    let site = sites
                        .get(*block.get(local).ok_or_eyre("Missing phase block site")?)
                        .ok_or_eyre("Missing phase call site")?;
                    let selected = best.first.get(local).ok_or_eyre("Missing phased base")?;
                    let allele = site
                        .alleles
                        .iter()
                        .find(|allele| (allele.base == *selected) != flip)
                        .ok_or_eyre("Phased base outside genotype")?;
                    calls.push((site.pileup, PhaseCall { set, first: allele.allele }));
                }
            }
        }
    }
    let phased = calls.len();
    debug!(
        sites = sites.len(),
        fragments = fragments.len(),
        evaluated_blocks = evaluated,
        phased_sites = phased,
        pruned_phase_states = pruned,
        "PairHMM phased segment"
    );
    Ok(calls)
}

#[cfg(feature = "experimental-seqair")]
fn root(parents: &[usize], mut index: usize) -> color_eyre::eyre::Result<usize> {
    loop {
        let parent = *parents.get(index).ok_or_eyre("Phase component index outside segment")?;
        if parent == index {
            return Ok(index);
        }
        index = parent;
    }
}

#[cfg(feature = "experimental-seqair")]
fn unite(parents: &mut [usize], a: usize, b: usize) -> color_eyre::eyre::Result<()> {
    use color_eyre::eyre::OptionExt as _;
    let a = root(parents, a)?;
    let b = root(parents, b)?;
    *parents.get_mut(a.max(b)).ok_or_eyre("Missing phase component root")? = a.min(b);
    Ok(())
}

#[cfg(feature = "experimental-seqair")]
use crate::{
    call::variant_calling::GenotypeTag,
    metrics::{AltCall, methylation::cpg_origin},
};
#[cfg(feature = "experimental-seqair")]
#[derive(Debug, Clone, Copy)]
struct HetAllele {
    allele: GtAllele,
    base: Base,
}

#[cfg(feature = "experimental-seqair")]
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
#[cfg(feature = "experimental-seqair")]
#[derive(Debug, Clone)]
struct CalledSite {
    /// Where the column sits in the slice [`build_problem`] was given.
    pub pileup: usize,
    pub pos: Pos0,
    /// The genotype's two alleles, ascending.
    alleles: [HetAllele; 2],
    /// The CpG cytosine one allele is, when conversion makes it read as the
    /// other: a C/T site's top strand or a G/A site's bottom strand.
    methylable: Option<CpgSide>,
}

#[cfg(feature = "experimental-seqair")]
impl CalledSite {
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
}
