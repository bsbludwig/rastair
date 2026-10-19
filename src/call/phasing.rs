//! Fragment likelihood phasing of called SNP genotypes. Full read evidence is
//! captured before seqair recycles its store; chromosome assignments are
//! scored jointly with shared CpG rates and local `PairHMM` realignment.

use crate::metrics::{PileupMetrics, methylation::CpgSide};
use color_eyre::eyre::OptionExt as _;
use seqair_types::{Base, Pos0, Pos1};
use std::num::{NonZeroU8, NonZeroU64};
use tracing::debug;

mod mask;
pub use mask::PhaseMask;

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

/// Where a heterozygous site sits in its phase block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PhaseCall {
    /// VCF `PS`: the block's first emitted heterozygous site.
    pub set: PhaseSet,
    /// The allele haplotype 1 carries.
    pub first: GtAllele,
}

/// What the molecular models keep on one column: the regional SNV evidence
/// until the genotype is assigned, then the record's place in a phase block,
/// if it has one. Indel records are never phased.
#[derive(Debug, Clone)]
#[expect(
    clippy::large_enum_variant,
    reason = "Columns already box phase data; inline features avoid another allocation per ML candidate"
)]
pub enum ColumnPhase {
    /// A count-genotyped column's phase.
    Phased(PhaseCall),
    Likelihood(crate::call::haplotypes::features::ColumnEvidence),
    /// A molecular genotype's column: its phase, once phased.
    Genotyped(Option<PhaseCall>),
}

impl ColumnPhase {
    /// Record `call` for the column's SNV record.
    fn set_call(slot: &mut Option<Box<Self>>, call: PhaseCall) {
        // A molecular genotype's column keeps saying so: `GTM` reads it.
        let genotyped = matches!(slot.as_deref(), Some(Self::Genotyped(_)));
        *slot = Some(Box::new(if genotyped {
            Self::Genotyped(Some(call))
        } else {
            Self::Phased(call)
        }));
    }
}

/// The phase call a column's SNV record carries, if any.
pub(crate) fn call_of(pileup: &PileupMetrics) -> Option<PhaseCall> {
    match pileup.phase.as_deref() {
        Some(ColumnPhase::Phased(call)) => Some(*call),
        Some(ColumnPhase::Genotyped(call)) => *call,
        Some(ColumnPhase::Likelihood(_)) | None => None,
    }
}

#[cfg(test)]
pub(crate) mod fixtures;

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
    segment: &crate::sequence::Segment,
    fragments: &[crate::call::haplotypes::Fragment],
    mask: &PhaseMask,
    profile: &crate::call::haplotypes::ErrorProfile,
) -> color_eyre::eyre::Result<()> {
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
        calls.extend(phase_unmasked(pileups, segment, fragments, region, profile)?);
    }
    for pileup in pileups.iter_mut() {
        if let Some(ColumnPhase::Genotyped(phase)) = pileup.phase.as_deref_mut() {
            *phase = None;
        } else {
            pileup.phase = None;
        }
    }
    for (index, call) in calls {
        let pileup = pileups.get_mut(index).ok_or_eyre("Phase call outside segment")?;
        ColumnPhase::set_call(&mut pileup.phase, call);
    }
    Ok(())
}

fn phase_unmasked(
    pileups: &[PileupMetrics],
    segment: &crate::sequence::Segment,
    fragments: &[crate::call::haplotypes::Fragment],
    region: std::ops::Range<u64>,
    profile: &crate::call::haplotypes::ErrorProfile,
) -> color_eyre::eyre::Result<Vec<(usize, PhaseCall)>> {
    use crate::call::haplotypes::{
        Reference, Scorer, Site,
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
        .iter()
        .filter(|((a, b), _)| {
            [a, b].into_iter().all(|&index| {
                linkage
                    .counts
                    .get(index)
                    .is_some_and(|counts| counts.iter().all(|&count| count >= 2))
            })
        })
        .map(|(&pair, &count)| (pair, count))
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
    let mut params = parameters()?;
    params.error_profile = profile.clone();
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
            let included: Vec<_> = linkage
                .fragments_covering(block.iter().copied())
                .into_iter()
                .filter_map(|index| fragments.get(index))
                .collect();
            let first = domains
                .first()
                .and_then(|site| site.alleles.first())
                .copied()
                .ok_or_eyre("Missing anchor allele")?;
            // Every beam step rescores the same fragments against the same
            // read-local templates; the memo lives as long as the block does.
            let mut scorer = Scorer::new(&params)?;
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
                    let prefix_reads: Vec<_> = linkage
                        .fragments_covering(
                            block.get(..=step).ok_or_eyre("Missing phase prefix")?.iter().copied(),
                        )
                        .into_iter()
                        .filter_map(|index| fragments.get(index))
                        .collect();
                    let posterior =
                        score(&reference, prefix, &states, &prefix_reads, &params, &mut scorer)?;
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
            let posterior = score(&reference, &domains, &states, &included, &params, &mut scorer)?;
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

fn root(parents: &[usize], mut index: usize) -> color_eyre::eyre::Result<usize> {
    loop {
        let parent = *parents.get(index).ok_or_eyre("Phase component index outside segment")?;
        if parent == index {
            return Ok(index);
        }
        index = parent;
    }
}

fn unite(parents: &mut [usize], a: usize, b: usize) -> color_eyre::eyre::Result<()> {
    let a = root(parents, a)?;
    let b = root(parents, b)?;
    *parents.get_mut(a.max(b)).ok_or_eyre("Missing phase component root")? = a.min(b);
    Ok(())
}

use crate::{
    call::variant_calling::GenotypeTag,
    metrics::{AltCall, methylation::cpg_origin},
};
#[derive(Debug, Clone, Copy)]
struct HetAllele {
    allele: GtAllele,
    base: Base,
}

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

/// A heterozygous column the phaser can place.
#[derive(Debug, Clone)]
struct CalledSite {
    /// Where the column sits in the slice it was built from.
    pub pileup: usize,
    pub pos: Pos0,
    /// The genotype's two alleles, ascending.
    alleles: [HetAllele; 2],
    /// The CpG cytosine one allele is, when conversion makes it read as the
    /// other: a C/T site's top strand or a G/A site's bottom strand.
    methylable: Option<CpgSide>,
}

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_zero_hash_is_no_fragment() {
        assert!(FragmentId::new(0).is_none());
        assert!(FragmentId::new(7).is_some());
    }
}
