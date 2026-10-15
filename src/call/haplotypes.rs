use color_eyre::eyre::{OptionExt as _, Result, ensure};
use seqair_types::{Base, Probability, Strand};
use std::collections::BTreeSet;

mod likelihood;
pub use likelihood::{LikelihoodEngine, LikelihoodParams};
mod expansion;
pub use expansion::ExpansionParams;
pub(crate) mod features;
pub(crate) mod genotype;
mod methylation;
mod pair_hmm;
#[cfg(any(feature = "experimental-seqair", test))]
pub(crate) mod phase;
pub(crate) mod reads;
#[cfg(feature = "experimental-seqair")]
pub(crate) use pair_hmm::Scorer;
pub use pair_hmm::{Indel, SequenceEvidence, SequenceObservation, SequencePosition};
mod diagnostics;
pub use diagnostics::{DiagnosticReadParams, HaplotypeDiagnosticParams, diagnose};

#[cfg(test)]
use likelihood::cpg_anchor;
#[cfg(test)]
use likelihood::fragment_likelihood;
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
    pub annotations: Option<Box<reads::ReadAnnotations>>,
}

impl ReadEvidence {
    fn covers(&self, pos: u32, engine: LikelihoodEngine) -> bool {
        if engine == LikelihoodEngine::PairHmm
            && let Some(sequence) = &self.sequence
        {
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
    pub indel: Option<Indel>,
}

impl Haplotype {
    pub fn allele_at(&self, pos: u32) -> Option<Base> {
        self.alleles.iter().find(|(p, _)| *p == pos).map(|(_, base)| *base)
    }

    pub fn base_at(&self, reference: &Reference, pos: u32) -> Option<Base> {
        if self.indel.as_ref().is_some_and(|indel| indel.deletes(pos)) {
            return None;
        }
        self.allele_at(pos).or_else(|| reference.base_at(pos))
    }

    fn adjacent_base(&self, reference: &Reference, pos: u32, forward: bool) -> Option<(u32, Base)> {
        let mut adjacent =
            if forward { pos.checked_add(1)? } else { pos.checked_sub(1).unwrap_or(pos) };
        if let Some(indel) = &self.indel {
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

#[derive(Debug, Clone, Copy)]
pub struct AlleleEvidence {
    pub base: Base,
    pub log_bayes_factor: f64,
    pub posterior: Probability,
    pub heterozygous: Probability,
    pub homozygous: Probability,
    pub expected_fraction: f64,
}

#[derive(Debug, Clone)]
pub struct SnpEvidence {
    pub pos: u32,
    pub reference: Base,
    pub alleles: Vec<AlleleEvidence>,
    pub genotype_entropy: f64,
    /// Unphased A/C/G/T genotypes in VCF triangular order, independent of REF/ALT order.
    pub genotype_posteriors: [Probability; 10],
    pub haplotypes: Vec<Haplotype>,
    pub pairs: Vec<PairEvidence>,
    pub sites: Vec<Site>,
    pub fragments: usize,
    pub focal_fragments: usize,
    pub focal_observations: usize,
    pub focal_read_observations: usize,
    pub focal_overlapping_fragments: usize,
    pub focal_discordant_fragments: usize,
    /// Candidate sites excluded by connectivity or resource limits.
    pub omitted_sites: usize,
    pub disconnected_sites: usize,
    pub resource_limited_sites: usize,
    pub pruned_haplotypes: usize,
    pub rescued_haplotypes: usize,
    /// Maximum discarded pair-posterior mass in a single pruning step, conditional on that step's proposals.
    pub maximum_pruned_mass: f64,
    pub expansion_steps: usize,
    pub offered_linking_fragments: usize,
    pub retained_linking_fragments: usize,
    pub realigned_fragments: usize,
    pub aligned_fallback_fragments: usize,
    pub omitted_indels: usize,
}

impl SnpEvidence {
    pub fn allele(&self, base: Base) -> Option<&AlleleEvidence> {
        self.alleles.iter().find(|allele| allele.base == base)
    }

    /// Conditional on this context's retained candidate haplotypes.
    pub fn genotypes_at(&self, pos: u32) -> Result<[Probability; 10]> {
        ensure!(
            self.sites.iter().any(|site| site.pos == pos),
            "Position is outside the inferred context"
        );
        let mut values = [0.0; 10];
        for pair in &self.pairs {
            let (a, b) = pair_bases(pair, &self.haplotypes, pos)?;
            let index = genotype_index(a, b).ok_or_eyre("Unknown context genotype")?;
            *values.get_mut(index).ok_or_eyre("Invalid context genotype")? += *pair.posterior;
        }
        let mut result = [Probability::ZERO; 10];
        for (result, value) in result.iter_mut().zip(values) {
            *result = probability(value)?;
        }
        Ok(result)
    }
}

/// All focal A/C/G/T hypotheses survive bounded expansion, including latent converted alleles.
#[cfg(test)]
pub fn infer_snp(
    reference: &Reference,
    focal: u32,
    neighbors: &[Site],
    fragments: &[Fragment],
    params: &LikelihoodParams,
    max_sites: usize,
) -> Result<SnpEvidence> {
    infer_snp_with_indels(reference, focal, neighbors, &[], fragments, params, max_sites)
}

#[cfg(test)]
pub fn infer_snp_with_indels(
    reference: &Reference,
    focal: u32,
    neighbors: &[Site],
    indels: &[Indel],
    fragments: &[Fragment],
    params: &LikelihoodParams,
    max_sites: usize,
) -> Result<SnpEvidence> {
    let search =
        ExpansionParams { max_sites: u16::try_from(max_sites)?, ..ExpansionParams::default() };
    infer(reference, focal, neighbors, indels, fragments, params, &search, false)
}

pub fn infer_snp_expanded(
    reference: &Reference,
    focal: u32,
    neighbors: &[Site],
    indels: &[Indel],
    fragments: &[Fragment],
    params: &LikelihoodParams,
    search: &ExpansionParams,
) -> Result<SnpEvidence> {
    infer(reference, focal, neighbors, indels, fragments, params, search, true)
}

#[expect(
    clippy::too_many_arguments,
    reason = "The search and scoring boundary share validated evidence and configuration"
)]
fn infer(
    reference: &Reference,
    focal: u32,
    neighbors: &[Site],
    indels: &[Indel],
    fragments: &[Fragment],
    params: &LikelihoodParams,
    search: &ExpansionParams,
    require_links: bool,
) -> Result<SnpEvidence> {
    params.validate()?;
    search.validate()?;
    let ref_base = reference.base_at(focal).ok_or_eyre("Focal site outside reference")?;
    let focal_site = Site::new(focal, ref_base, Base::KNOWN)?;
    let mut candidates: Vec<_> =
        neighbors.iter().filter(|site| site.pos != focal).cloned().collect();
    if require_links {
        candidates.retain(|site| site.pos.abs_diff(focal) <= search.expansion_radius);
    }
    candidates.sort_by_key(|site| (site.pos.abs_diff(focal), site.pos));
    ensure!(
        candidates.iter().map(|site| site.pos).collect::<BTreeSet<_>>().len() == candidates.len(),
        "Duplicate candidate positions"
    );
    for site in &candidates {
        ensure!(
            reference.base_at(site.pos) == Some(site.reference),
            "Candidate reference mismatch"
        );
        ensure!(
            site.alleles.len() >= 2 && site.alleles.contains(&site.reference),
            "Invalid neighboring alleles"
        );
        ensure!(
            site.alleles.iter().all(|base| base.known_index().is_some()),
            "Unknown neighboring allele"
        );
        ensure!(
            site.alleles.iter().map(|base| *base as u8).collect::<BTreeSet<_>>().len()
                == site.alleles.len(),
            "Duplicate neighboring allele"
        );
    }
    let offered = candidates.len();
    let linking_fragments = |sites: &[Site]| {
        fragments
            .iter()
            .filter(|fragment| {
                fragment.reads.iter().any(|read| read.covers(focal, params.likelihood_engine))
                    && fragment.reads.iter().any(|read| {
                        sites.iter().any(|site| {
                            site.pos != focal && read.covers(site.pos, params.likelihood_engine)
                        })
                    })
            })
            .count()
    };
    let offered_linking_fragments = linking_fragments(&candidates);
    let mut proposed_indels = Vec::new();
    if params.likelihood_engine == LikelihoodEngine::PairHmm {
        for indel in indels {
            let indel = indel.normalized(reference)?;
            if !proposed_indels.contains(&indel) {
                proposed_indels.push(indel);
            }
        }
        proposed_indels.sort_by_key(|indel| (indel.start.abs_diff(focal), indel.start));
    }
    let selected_indel = proposed_indels
        .iter()
        .find(|indel| {
            indel.len() <= usize::from(params.max_indel_length)
                && !indel.deletes(focal)
                && candidates.iter().all(|site| !indel.deletes(site.pos))
                && (!require_links
                    || indel.start.abs_diff(focal) <= search.window_radius
                    || fragments
                        .iter()
                        .filter(|fragment| {
                            fragment
                                .reads
                                .iter()
                                .any(|read| read.covers(focal, params.likelihood_engine))
                                && fragment
                                    .reads
                                    .iter()
                                    .any(|read| read.covers(indel.start, params.likelihood_engine))
                        })
                        .take(usize::from(search.min_linking_fragments))
                        .count()
                        >= usize::from(search.min_linking_fragments))
        })
        .cloned();
    let expanded = expansion::expand(
        reference,
        focal_site,
        &candidates,
        selected_indel.as_ref(),
        fragments,
        params,
        search,
        require_links,
    )?;
    let sites = expanded.sites;
    let haplotypes = expanded.haplotypes;
    let matrix = likelihood_matrix_geometry(
        reference,
        &sites,
        &haplotypes,
        selected_indel.as_ref(),
        fragments,
        params,
        search.read_local,
    )?;
    let pairs =
        posterior_pairs(&haplotypes, &sites, &matrix.values, params, selected_indel.is_some())?;
    let mut genotypes = [0.0; 10];
    for pair in &pairs {
        let (a, b) = pair_bases(pair, &haplotypes, focal)?;
        let index = genotype_index(a, b).ok_or_eyre("Unknown focal genotype")?;
        *genotypes.get_mut(index).ok_or_eyre("Invalid genotype index")? += *pair.posterior;
    }
    let entropy = -genotypes.iter().filter(|&&p| p > 0.0).map(|p| p * p.ln()).sum::<f64>();
    let mut alleles = Vec::with_capacity(3);
    for base in Base::KNOWN.into_iter().filter(|base| *base != ref_base) {
        let (mut present_likelihood, mut absent_likelihood, mut present_prior, mut absent_prior) =
            (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        let (mut present, mut heterozygous, mut homozygous, mut fraction) = (0.0, 0.0, 0.0, 0.0);
        for pair in &pairs {
            let (a, b) = pair_bases(pair, &haplotypes, focal)?;
            let copies = u8::from(a == base) + u8::from(b == base);
            if copies == 0 {
                absent_likelihood =
                    log_add(absent_likelihood, pair.log_likelihood + pair.log_prior);
                absent_prior = log_add(absent_prior, pair.log_prior);
            } else {
                present_likelihood =
                    log_add(present_likelihood, pair.log_likelihood + pair.log_prior);
                present_prior = log_add(present_prior, pair.log_prior);
                present += *pair.posterior;
                fraction += *pair.posterior * f64::from(copies) / 2.0;
                if copies == 1 {
                    heterozygous += *pair.posterior;
                } else {
                    homozygous += *pair.posterior;
                }
            }
        }
        alleles.push(AlleleEvidence {
            base,
            log_bayes_factor: (present_likelihood - present_prior)
                - (absent_likelihood - absent_prior),
            posterior: probability(present)?,
            heterozygous: probability(heterozygous)?,
            homozygous: probability(homozygous)?,
            expected_fraction: fraction,
        });
    }
    let mut focal_observations = 0;
    let mut focal_read_observations = 0;
    let mut focal_overlapping_fragments = 0;
    let mut focal_discordant_fragments = 0;
    for fragment in fragments {
        let base = likelihood::consensus_at(fragment, focal)?;
        focal_observations += usize::from(base.source.is_some());
        focal_read_observations += base.measurements;
        focal_overlapping_fragments += usize::from(base.measurements == 2);
        focal_discordant_fragments += usize::from(base.discordant);
    }
    let focal_fragments = fragments
        .iter()
        .filter(|fragment| {
            fragment.reads.iter().any(|read| read.covers(focal, params.likelihood_engine))
        })
        .count();
    let used = sites.len() - 1;
    let mut genotype_posteriors = [Probability::ZERO; 10];
    for (posterior, value) in genotype_posteriors.iter_mut().zip(genotypes) {
        *posterior = probability(value)?;
    }
    Ok(SnpEvidence {
        pos: focal,
        reference: ref_base,
        alleles,
        genotype_entropy: entropy,
        genotype_posteriors,
        retained_linking_fragments: linking_fragments(&sites),
        haplotypes,
        pairs,
        sites,
        fragments: matrix.fragments,
        focal_fragments,
        focal_observations,
        focal_read_observations,
        focal_overlapping_fragments,
        focal_discordant_fragments,
        omitted_sites: offered - used,
        disconnected_sites: expanded.disconnected_sites,
        resource_limited_sites: expanded.resource_limited_sites,
        pruned_haplotypes: expanded.pruned_haplotypes,
        rescued_haplotypes: expanded.rescued_haplotypes,
        maximum_pruned_mass: expanded.maximum_pruned_mass,
        expansion_steps: expanded.steps,
        offered_linking_fragments,
        realigned_fragments: matrix.realigned_fragments,
        aligned_fallback_fragments: matrix.aligned_fallback_fragments,
        omitted_indels: proposed_indels.len() - usize::from(selected_indel.is_some()),
    })
}

struct LikelihoodMatrix {
    values: Vec<Vec<f64>>,
    fragments: usize,
    realigned_fragments: usize,
    aligned_fallback_fragments: usize,
}

fn likelihood_matrix_geometry(
    reference: &Reference,
    sites: &[Site],
    haplotypes: &[Haplotype],
    indel: Option<&Indel>,
    fragments: &[Fragment],
    params: &LikelihoodParams,
    local: bool,
) -> Result<LikelihoodMatrix> {
    let loci = scored_loci(reference, sites);
    let indel_span = if local {
        let mut span: Option<std::ops::Range<u32>> = None;
        let edits = haplotypes
            .iter()
            .filter_map(|haplotype| haplotype.indel.as_ref())
            .chain(indel)
            .collect::<std::collections::BTreeSet<_>>();
        for edit in edits {
            let end = edit.equivalent_end(reference)?;
            match &mut span {
                Some(span) => {
                    span.start = span.start.min(edit.start);
                    span.end = span.end.max(end);
                }
                None => span = Some(edit.start..end),
            }
        }
        span
    } else {
        None
    };
    let included: Vec<_> = fragments
        .iter()
        .filter(|fragment| {
            fragment.reads.iter().any(|read| {
                loci.iter().any(|&pos| read.covers(pos, params.likelihood_engine))
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
    let realignment = if params.likelihood_engine == LikelihoodEngine::PairHmm {
        let mut geometry = sites.to_vec();
        if let Some(span) = &indel_span {
            let pos = span.end.saturating_sub(1);
            if let Some(base) = reference.base_at(pos) {
                // A normalized gap can stand for placements throughout a repeat.
                // Include the right boundary before clipping alignment profiles.
                geometry.push(Site { pos, reference: base, alleles: vec![base] });
            }
        }
        let prepared = pair_hmm::Prepared::new(reference, &geometry, indel, &included, params)?;
        Some(if local { prepared.localize() } else { prepared })
    } else {
        None
    };
    let context = methylation::Context::new(
        reference,
        haplotypes,
        &loci,
        realignment.as_ref().map(pair_hmm::Prepared::interval),
    )?;
    let templates = haplotypes
        .iter()
        .map(|haplotype| {
            realignment
                .as_ref()
                .map(|prepared| prepared.template(reference, haplotype, &loci, params))
                .transpose()
        })
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
    let mut scorer = pair_hmm::Scorer::new(params);
    // One pass over the matrix: every realigned pair is planned first, so the
    // forward passes run as one batch. In read-local mode haplotypes that agree
    // inside a fragment's window share its first plan, as they shared a cache.
    let mut pass = |rates: Option<&methylation::Rates>, neutral: bool| -> Result<Vec<Vec<f64>>> {
        enum Slot {
            Planned(usize),
            Value(f64),
        }
        let mut plans = Vec::new();
        let mut distinct: Vec<Slot> = Vec::new();
        let mut cache =
            std::collections::BTreeMap::<(usize, Vec<(u32, u8)>, Option<&Indel>), usize>::new();
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
                let key = if local {
                    realignment.as_ref().map(|prepared| {
                        (
                            fragment_index,
                            prepared.signature_at(fragment_index, haplotype, params),
                            haplotype.indel.as_ref(),
                        )
                    })
                } else {
                    None
                };
                if let Some(key) = &key
                    && let Some(&slot) = cache.get(key)
                {
                    row.push(slot);
                    continue;
                }
                let planned = match (&realignment, templates.get(haplotype_index)) {
                    (Some(prepared), Some(Some(template))) => {
                        prepared.plan(fragment_index, fragment, template, params, &scoring)?
                    }
                    _ => None,
                };
                let slot = match planned {
                    Some(plan) => {
                        plans.push(plan);
                        Slot::Planned(plans.len() - 1)
                    }
                    None => Slot::Value(aligned(haplotype_index, fragment_index, &scoring)?),
                };
                distinct.push(slot);
                if let Some(key) = key {
                    cache.insert(key, distinct.len() - 1);
                }
                row.push(distinct.len() - 1);
            }
            slots.push(row);
        }
        let scored = scorer.run(&plans, false)?;
        slots
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
            .collect()
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
        realigned_fragments: realignment
            .as_ref()
            .map_or(0, |prepared| prepared.realigned_fragments()),
        aligned_fallback_fragments: realignment
            .as_ref()
            .map_or(0, |prepared| prepared.fallback_fragments()),
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

#[cfg(test)]
fn enumerate(sites: &[Site]) -> Vec<Haplotype> {
    let mut haplotypes = vec![Haplotype { alleles: Vec::new(), indel: None }];
    for site in sites {
        haplotypes = haplotypes
            .iter()
            .flat_map(|haplotype| {
                site.alleles.iter().map(|&base| {
                    let mut haplotype = haplotype.clone();
                    haplotype.alleles.push((site.pos, base));
                    haplotype
                })
            })
            .collect();
    }
    haplotypes
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
