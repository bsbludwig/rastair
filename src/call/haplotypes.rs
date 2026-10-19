//! The molecular likelihood model: fragments of captured read pairs scored
//! against candidate chromosomes by a TAPS-aware pair-HMM with shared CpG
//! methylation rates. `regional` searches diploid genotypes over a region's
//! observed alleles for the SNV forests and the molecular genotype, `phase`
//! links called heterozygotes into phase blocks.

use color_eyre::eyre::{OptionExt as _, Result, ensure};
use seqair_types::{Base, Probability, Strand};
use std::collections::BTreeSet;

mod error_model;
mod regional;
pub use error_model::ErrorProfile;
mod likelihood;
pub use likelihood::LikelihoodParams;
pub(crate) mod features;
pub(crate) mod genotype;
mod methylation;
mod pair_hmm;
pub(crate) mod phase;
pub(crate) mod reads;
pub(crate) use pair_hmm::Scorer;
pub use pair_hmm::{Indel, SequenceEvidence, SequenceObservation, SequencePosition};
pub use regional::ExpansionParams;

use likelihood::log_add;

#[derive(Debug, Clone)]
pub struct Reference {
    pub start: u32,
    pub bases: Vec<Base>,
}

impl Reference {
    pub fn base_at(&self, pos: u32) -> Option<Base> {
        self.bases.get(usize::try_from(pos.checked_sub(self.start)?).ok()?).copied()
    }
}

#[derive(Debug, Clone)]
pub struct Site {
    pub pos: u32,
    pub reference: Base,
    pub alleles: Vec<Base>,
}

impl Site {
    pub fn new(pos: u32, reference: Base, alleles: impl IntoIterator<Item = Base>) -> Result<Self> {
        ensure!(reference.known_index().is_some(), "Unknown reference cannot be a SNP candidate");
        let mut bases = vec![reference];
        for base in alleles {
            ensure!(base.known_index().is_some(), "Unknown base cannot be a SNP allele");
            if !bases.contains(&base) {
                bases.push(base);
            }
        }
        ensure!(bases.len() > 1, "A candidate needs a competing allele");
        Ok(Self { pos, reference, alleles: bases })
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Observation {
    pub pos: u32,
    pub base: Base,
    pub quality: Option<u8>,
}

#[derive(Debug, Clone)]
pub struct ReadEvidence {
    pub mapq: Option<u8>,
    pub strand: Strand,
    pub observations: Vec<Observation>,
    pub sequence: Option<SequenceEvidence>,
}

impl ReadEvidence {
    /// Whether the read reaches `pos`: anywhere in its aligned span when its
    /// sequence was captured, else at an observed base.
    fn covers(&self, pos: u32) -> bool {
        if let Some(sequence) = &self.sequence {
            return sequence.bases.first().zip(sequence.bases.last()).is_some_and(
                |(first, last)| first.position.pos <= pos && pos <= last.position.pos,
            );
        }
        self.observations.binary_search_by_key(&pos, |observation| observation.pos).is_ok()
    }
}

#[derive(Debug, Clone)]
pub struct Fragment {
    pub reads: Vec<ReadEvidence>,
}

#[derive(Debug, Clone)]
pub struct Haplotype {
    pub alleles: Vec<(u32, Base)>,
    pub indels: Vec<Indel>,
}

impl Haplotype {
    pub fn allele_at(&self, pos: u32) -> Option<Base> {
        self.alleles.iter().find(|(p, _)| *p == pos).map(|(_, base)| *base)
    }

    pub fn base_at(&self, reference: &Reference, pos: u32) -> Option<Base> {
        if self.indels.iter().any(|indel| indel.deletes(pos)) {
            return None;
        }
        self.allele_at(pos).or_else(|| reference.base_at(pos))
    }

    fn adjacent_base(&self, reference: &Reference, pos: u32, forward: bool) -> Option<(u32, Base)> {
        let mut adjacent =
            if forward { pos.checked_add(1)? } else { pos.checked_sub(1).unwrap_or(pos) };
        for indel in &self.indels {
            if (forward && indel.start == adjacent) || (!forward && indel.start == pos) {
                let inserted = if forward { indel.inserted.first() } else { indel.inserted.last() };
                if let Some(&base) = inserted {
                    return Some((indel.start, base));
                }
            }
            if indel.deletes(adjacent) {
                adjacent = if forward {
                    indel.start.checked_add(u32::try_from(indel.deleted.len()).ok()?)?
                } else {
                    indel.start.checked_sub(1)?
                };
            }
        }
        if adjacent == pos {
            return None;
        }
        self.base_at(reference, adjacent).map(|base| (adjacent, base))
    }
}

#[derive(Debug, Clone)]
pub struct PairEvidence {
    pub first: usize,
    pub second: usize,
    pub log_likelihood: f64,
    pub log_prior: f64,
    pub posterior: Probability,
}

/// Every compared chromosome's log likelihood of every included fragment,
/// row-major by chromosome. Fragments not realigned were scored by the
/// aligned model, because their mates' alignments could not be joined into
/// one profile.
struct LikelihoodMatrix {
    values: Vec<Vec<f64>>,
    fragments: usize,
    realigned_fragments: usize,
}

/// Score `haplotypes` against the fragments reaching the sites' loci or the
/// span of a compared indel: a neutral pilot pass, the shared CpG rates it
/// supports, and the pass under those rates.
fn likelihood_matrix_geometry(
    reference: &Reference,
    sites: &[Site],
    haplotypes: &[Haplotype],
    indel: Option<&Indel>,
    fragments: &[Fragment],
    params: &LikelihoodParams,
    scorer: &mut pair_hmm::Scorer,
) -> Result<LikelihoodMatrix> {
    let loci = scored_loci(reference, sites);
    let mut indel_span: Option<std::ops::Range<u32>> = None;
    let edits = haplotypes
        .iter()
        .flat_map(|haplotype| &haplotype.indels)
        .chain(indel)
        .collect::<std::collections::BTreeSet<_>>();
    for edit in edits {
        let end = edit.equivalent_end(reference)?;
        match &mut indel_span {
            Some(span) => {
                span.start = span.start.min(edit.start);
                span.end = span.end.max(end);
            }
            None => indel_span = Some(edit.start..end),
        }
    }
    let included: Vec<_> = fragments
        .iter()
        .filter(|fragment| {
            fragment.reads.iter().any(|read| {
                loci.iter().any(|&pos| read.covers(pos))
                    || indel_span.as_ref().is_some_and(|span| {
                        read.observations.iter().any(|obs| span.contains(&obs.pos))
                    })
            })
        })
        .collect();
    for fragment in &included {
        ensure!(
            !fragment.reads.is_empty() && fragment.reads.len() <= 2,
            "Fragments require one or two reads"
        );
        for read in &fragment.reads {
            ensure!(read.mapq != Some(255), "Missing mapping quality must be represented by None");
            ensure!(
                read.observations.iter().all(|obs| obs.quality != Some(255)),
                "Missing base quality must be represented by None"
            );
            ensure!(
                read.observations.windows(2).all(|pair| match pair {
                    [a, b] => a.pos < b.pos,
                    _ => true,
                }),
                "A read must have one observation per position in ascending order"
            );
        }
    }
    let mut geometry = sites.to_vec();
    if let Some(span) = &indel_span {
        let pos = span.end.saturating_sub(1);
        if let Some(base) = reference.base_at(pos) {
            // A normalized gap can stand for placements throughout a repeat.
            // Include the right boundary before clipping alignment profiles.
            geometry.push(Site { pos, reference: base, alleles: vec![base] });
        }
    }
    let realignment = pair_hmm::Prepared::new(reference, &geometry, indel, &included, params)?
        .reaching(haplotypes)?;
    let context = methylation::Context::new(reference, haplotypes, &loci, realignment.interval())?;
    let templates = haplotypes
        .iter()
        .map(|haplotype| realignment.template(reference, haplotype, &loci, params))
        .collect::<Result<Vec<_>>>()?;
    let aligned = |haplotype_index: usize,
                   fragment_index: usize,
                   scoring: &methylation::Scoring<'_>|
     -> Result<f64> {
        likelihood::fragment_likelihood_with_model(
            included.get(fragment_index).ok_or_eyre("Missing scoring fragment")?,
            haplotypes.get(haplotype_index).ok_or_eyre("Missing scoring haplotype")?,
            reference,
            &loci,
            params,
            scoring,
            context.haplotypes.get(haplotype_index).ok_or_eyre("Missing CpG context")?,
        )
    };
    // One pass over the matrix: every realigned pair is planned first, so the
    // forward passes run as one batch. Haplotypes that agree inside a
    // fragment's window share its first plan, as they shared a cache.
    let mut pass = |rates: Option<&methylation::Rates>, neutral: bool| -> Result<Vec<Vec<f64>>> {
        enum Slot {
            Planned(usize),
            Value(f64),
        }
        let mut plans = Vec::new();
        let mut distinct: Vec<Slot> = Vec::new();
        let mut cache =
            std::collections::BTreeMap::<(usize, Vec<(u32, u8)>, &[Indel]), usize>::new();
        let mut slots = Vec::with_capacity(haplotypes.len());
        for (haplotype_index, haplotype) in haplotypes.iter().enumerate() {
            let mut row = Vec::with_capacity(included.len());
            for (fragment_index, fragment) in included.iter().enumerate() {
                let scoring = if neutral {
                    methylation::Scoring::Neutral(&context)
                } else if let Some(rates) = rates {
                    methylation::Scoring::Shrinkage(rates, fragment_index)
                } else {
                    methylation::Scoring::Ordinary
                };
                let key = (
                    fragment_index,
                    realignment.signature_at(fragment_index, haplotype, params),
                    haplotype.indels.as_slice(),
                );
                if let Some(&slot) = cache.get(&key) {
                    row.push(slot);
                    continue;
                }
                let template = templates.get(haplotype_index).ok_or_eyre("Missing template")?;
                let slot =
                    match realignment.plan(fragment_index, fragment, template, params, &scoring)? {
                        Some(plan) => {
                            plans.push(plan);
                            Slot::Planned(plans.len() - 1)
                        }
                        None => Slot::Value(aligned(haplotype_index, fragment_index, &scoring)?),
                    };
                distinct.push(slot);
                cache.insert(key, distinct.len() - 1);
                row.push(distinct.len() - 1);
            }
            slots.push(row);
        }
        // The shrinkage pass's rates change with every call, so only the
        // rate-free passes are worth remembering; the memo key includes the
        // rates, so remembering either is exact.
        let scored = scorer.run(&plans, rates.is_none())?;
        let values = slots
            .iter()
            .map(|row| {
                row.iter()
                    .map(|&slot| {
                        match distinct.get(slot).ok_or_eyre("Missing likelihood slot")? {
                            Slot::Planned(plan) => {
                                scored.get(*plan).copied().ok_or_eyre("Missing planned likelihood")
                            }
                            Slot::Value(value) => Ok(*value),
                        }
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let mut values = values;
        error_model::cap_molecules(&mut values, &included, reference, params)?;
        Ok(values)
    };
    let values = if params.non_taps {
        pass(None, false)?
    } else {
        let pilot = pass(None, true)?;
        let rates = methylation::estimate(
            &context,
            haplotypes,
            sites,
            &included,
            &pilot,
            params,
            indel.is_some(),
        )?;
        pass(Some(&rates), false)?
    };
    Ok(LikelihoodMatrix {
        values,
        fragments: included.len(),
        realigned_fragments: realignment.realigned_fragments(),
    })
}

fn posterior_pairs(
    haplotypes: &[Haplotype],
    sites: &[Site],
    matrix: &[Vec<f64>],
    params: &LikelihoodParams,
    has_indel: bool,
) -> Result<Vec<PairEvidence>> {
    let priors = likelihood::retained_pair_priors(haplotypes, sites, params, has_indel)?;
    let mut pairs = Vec::new();
    for (first, _) in haplotypes.iter().enumerate() {
        let a = matrix.get(first).ok_or_eyre("Missing haplotype likelihoods")?;
        for (second, _) in haplotypes.iter().enumerate().skip(first) {
            let b = matrix.get(second).ok_or_eyre("Missing haplotype likelihoods")?;
            ensure!(a.len() == b.len(), "Haplotype likelihood dimensions disagree");
            let value = if first == second {
                a.iter().sum::<f64>()
            } else {
                a.iter().zip(b).map(|(&a, &b)| log_add(a, b) - std::f64::consts::LN_2).sum()
            };
            ensure!(value.is_finite(), "Pair likelihood is not finite");
            pairs.push(PairEvidence {
                first,
                second,
                log_likelihood: value,
                log_prior: *priors.get(pairs.len()).ok_or_eyre("Missing pair prior")?,
                posterior: Probability::ZERO,
            });
        }
    }
    let maximum = pairs
        .iter()
        .map(|pair| pair.log_likelihood + pair.log_prior)
        .reduce(f64::max)
        .ok_or_eyre("No haplotype pairs")?;
    let weights: Vec<_> =
        pairs.iter().map(|pair| (pair.log_likelihood + pair.log_prior - maximum).exp()).collect();
    let normalization: f64 = weights.iter().sum();
    ensure!(
        normalization.is_finite() && normalization > 0.0,
        "Invalid haplotype pair normalization"
    );
    // Normalizing relative weights avoids subtracting a large negative log
    // normalizer, whose rounding can make aggregate probabilities exceed one.
    for (pair, weight) in pairs.iter_mut().zip(weights) {
        pair.posterior = probability(weight / normalization)?;
    }
    Ok(pairs)
}

fn probability(value: f64) -> Result<Probability> {
    ensure!(
        value.is_finite() && (-1e-12..=1.0 + 1e-12).contains(&value),
        "Invalid posterior {value}"
    );
    Probability::new(value.clamp(0.0, 1.0)).map_err(Into::into)
}

fn pair_bases(pair: &PairEvidence, haplotypes: &[Haplotype], pos: u32) -> Result<(Base, Base)> {
    let allele = |index: usize| {
        haplotypes
            .get(index)
            .and_then(|haplotype| haplotype.allele_at(pos))
            .ok_or_eyre("Missing focal haplotype allele")
    };
    Ok((allele(pair.first)?, allele(pair.second)?))
}

fn genotype_index(a: Base, b: Base) -> Option<usize> {
    let (a, b) = (a.known_index()?, b.known_index()?);
    let (low, high) = (a.min(b), a.max(b));
    Some(high * (high + 1) / 2 + low)
}

fn scored_loci(reference: &Reference, sites: &[Site]) -> Vec<u32> {
    let can_be = |pos, base| {
        sites.iter().find(|site| site.pos == pos).map_or_else(
            || reference.base_at(pos) == Some(base),
            |site| site.alleles.contains(&base),
        )
    };
    let mut loci = BTreeSet::new();
    for site in sites {
        loci.insert(site.pos);
        if site.alleles.contains(&Base::C)
            && let Some(pos) = site.pos.checked_add(1)
            && can_be(pos, Base::G)
        {
            loci.insert(pos);
        }
        if site.alleles.contains(&Base::G)
            && let Some(pos) = site.pos.checked_sub(1)
            && can_be(pos, Base::C)
        {
            loci.insert(pos);
        }
    }
    loci.into_iter().collect()
}

#[cfg(test)]
mod tests;
