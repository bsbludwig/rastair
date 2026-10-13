use color_eyre::eyre::{OptionExt as _, Result, ensure};
use seqair_types::{Base, Probability, Strand};
use std::collections::BTreeSet;

mod likelihood;
pub use likelihood::{LikelihoodEngine, LikelihoodParams, MethylationModel};
mod methylation;
mod pair_hmm;
pub use pair_hmm::{Indel, SequenceEvidence, SequenceObservation, SequencePosition};
mod diagnostics;
pub use diagnostics::{DiagnosticReadParams, HaplotypeDiagnosticParams, diagnose};

#[cfg(test)]
use likelihood::fragment_likelihood;
use likelihood::{Node, cpg_anchor, log_add, log_sum, pair_prior, quadrature};

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
    /// Other candidates were omitted to keep the Cartesian search bounded.
    pub omitted_sites: usize,
    /// The bounded prototype integrates at most one CpG rate per haplotype.
    /// It reduces context rather than tying unrelated CpG rates together.
    pub cpg_context_limited: bool,
    pub cpg_omitted_sites: usize,
    pub cpg_omitted_indels: usize,
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
}

/// SNP inference under the selected read likelihood engine. Every focal A/C/G/T hypothesis is
/// retained, even when no read directly shows that genetic base. Methylation
/// is marginalized under independent uniform priors for the CpGs on the two
/// haplotypes; the same rate is shared by all fragments observing that CpG.
///
/// Context is reduced until each haplotype has at most one scored CpG. This
/// makes the prototype's numerical integral at most two-dimensional without
/// silently sharing a rate between different CpGs. `PairHMM` flank CpGs use a
/// fixed prior-predictive methylation fraction rather than additional rate fits.
/// The neutral and shrinkage alternatives retain multiple CpGs; shrinkage uses
/// shared rates learned with a neutral pilot and excludes each scored molecule.
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

pub fn infer_snp_with_indels(
    reference: &Reference,
    focal: u32,
    neighbors: &[Site],
    indels: &[Indel],
    fragments: &[Fragment],
    params: &LikelihoodParams,
    max_sites: usize,
) -> Result<SnpEvidence> {
    params.validate()?;
    ensure!(
        (1..=3).contains(&max_sites),
        "Bounded SNP search supports one to three candidate sites"
    );
    let ref_base = reference.base_at(focal).ok_or_eyre("Focal site outside reference")?;
    let focal_site = Site::new(focal, ref_base, Base::KNOWN)?;
    let mut candidates: Vec<_> =
        neighbors.iter().filter(|site| site.pos != focal).cloned().collect();
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
    candidates.truncate(max_sites - 1);
    let bounded_sites = candidates.len();
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
    let mut limited = false;
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
    let mut selected_indel = proposed_indels
        .iter()
        .find(|indel| {
            indel.len() <= usize::from(params.max_indel_length)
                && !indel.deletes(focal)
                && candidates.iter().all(|site| !indel.deletes(site.pos))
        })
        .cloned();
    let offered_indel = selected_indel.is_some();
    let (sites, haplotypes, loci) = loop {
        let mut sites = candidates.clone();
        sites.push(focal_site.clone());
        sites.sort_by_key(|site| site.pos);
        let mut haplotypes = enumerate(&sites);
        if let Some(indel) = &selected_indel {
            let with_indel: Vec<_> = haplotypes
                .iter()
                .cloned()
                .map(|mut haplotype| {
                    haplotype.indel = Some(indel.clone());
                    haplotype
                })
                .collect();
            haplotypes.extend(with_indel);
        }
        let loci = scored_loci(reference, &sites);
        let multiple_cpgs = !params.non_taps
            && params.methylation_model == MethylationModel::Integrated
            && haplotypes
                .iter()
                .map(|haplotype| {
                    if params.likelihood_engine == LikelihoodEngine::PairHmm {
                        pair_hmm::tracked_cpgs(reference, haplotype, &loci)
                            .map(|anchors| anchors.len() > 1)
                    } else {
                        Ok(loci
                            .iter()
                            .filter_map(|&pos| cpg_anchor(reference, haplotype, pos))
                            .collect::<BTreeSet<_>>()
                            .len()
                            > 1)
                    }
                })
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .any(|multiple| multiple);
        if !multiple_cpgs {
            break (sites, haplotypes, loci);
        }
        if candidates.pop().is_none() {
            ensure!(
                selected_indel.take().is_some(),
                "A single SNP has multiple independent observed CpGs"
            );
        }
        limited = true;
    };
    let included: Vec<_> = fragments
        .iter()
        .filter(|fragment| {
            fragment
                .reads
                .iter()
                .any(|read| loci.iter().any(|&pos| read.covers(pos, params.likelihood_engine)))
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
        Some(pair_hmm::Prepared::new(
            reference,
            &sites,
            selected_indel.as_ref(),
            &included,
            params,
        )?)
    } else {
        None
    };
    let context = methylation::Context::new(
        reference,
        &haplotypes,
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
    let score = |haplotype_index: usize,
                 fragment_index: usize,
                 scoring: &methylation::Scoring<'_>|
     -> Result<f64> {
        let haplotype = haplotypes.get(haplotype_index).ok_or_eyre("Missing scoring haplotype")?;
        let fragment = included.get(fragment_index).ok_or_eyre("Missing scoring fragment")?;
        if let Some(prepared) = &realignment
            && let Some(template) = templates.get(haplotype_index).and_then(Option::as_ref)
            && let Some(value) = prepared.fragment_likelihood_with_model(
                fragment_index,
                fragment,
                template,
                params,
                scoring,
            )?
        {
            return Ok(value);
        }
        likelihood::fragment_likelihood_with_model(
            fragment,
            haplotype,
            reference,
            &loci,
            params,
            scoring,
            context.haplotypes.get(haplotype_index).ok_or_eyre("Missing CpG context")?,
        )
    };
    let rates = if !params.non_taps && params.methylation_model == MethylationModel::Shrinkage {
        let pilot = (0..haplotypes.len())
            .map(|haplotype| {
                (0..included.len())
                    .map(|fragment| {
                        score(haplotype, fragment, &methylation::Scoring::Neutral(&context))
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        Some(methylation::estimate(
            &context,
            &haplotypes,
            &sites,
            &included,
            &pilot,
            params,
            selected_indel.is_some(),
        )?)
    } else {
        None
    };
    let nodes = quadrature(params.quadrature_order)?;
    let mut matrix = Vec::with_capacity(haplotypes.len());
    for (haplotype_index, haplotype) in haplotypes.iter().enumerate() {
        let methylable = !params.non_taps
            && params.methylation_model == MethylationModel::Integrated
            && if realignment.is_some() {
                !pair_hmm::tracked_cpgs(reference, haplotype, &loci)?.is_empty()
            } else {
                loci.iter().any(|&pos| cpg_anchor(reference, haplotype, pos).is_some())
            };
        let nodes: Vec<_> = if methylable {
            nodes.clone()
        } else {
            vec![Node { methylation: 0.0, log_weight: 0.0 }]
        };
        let values = nodes
            .iter()
            .map(|node| {
                included
                    .iter()
                    .enumerate()
                    .map(|(index, _)| {
                        let scoring = match params.methylation_model {
                            MethylationModel::Neutral => methylation::Scoring::Neutral(&context),
                            MethylationModel::Shrinkage if rates.is_some() => {
                                methylation::Scoring::Shrinkage(
                                    rates.as_ref().ok_or_eyre("Missing estimated rates")?,
                                    index,
                                )
                            }
                            _ => methylation::Scoring::Integrated(node.methylation),
                        };
                        score(haplotype_index, index, &scoring)
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        matrix.push((nodes, values));
    }
    let mut pairs = Vec::new();
    for (first, a) in haplotypes.iter().enumerate() {
        let (nodes_a, values_a) = matrix.get(first).ok_or_eyre("Missing haplotype likelihoods")?;
        for (second, b) in haplotypes.iter().enumerate().skip(first) {
            let (nodes_b, values_b) =
                matrix.get(second).ok_or_eyre("Missing haplotype likelihoods")?;
            let mut integral = f64::NEG_INFINITY;
            for (node_a, likelihoods_a) in nodes_a.iter().zip(values_a) {
                if first == second {
                    integral =
                        log_add(integral, node_a.log_weight + likelihoods_a.iter().sum::<f64>());
                } else {
                    for (node_b, likelihoods_b) in nodes_b.iter().zip(values_b) {
                        let likelihood = likelihoods_a
                            .iter()
                            .zip(likelihoods_b)
                            .map(|(&a, &b)| log_add(a, b) - std::f64::consts::LN_2)
                            .sum::<f64>();
                        integral =
                            log_add(integral, node_a.log_weight + node_b.log_weight + likelihood);
                    }
                }
            }
            ensure!(integral.is_finite(), "Integrated pair likelihood is not finite");
            pairs.push(PairEvidence {
                first,
                second,
                log_likelihood: integral,
                log_prior: pair_prior(
                    a,
                    b,
                    &sites,
                    *params.snp_prior,
                    *params.indel_prior,
                    selected_indel.is_some(),
                )?,
                posterior: Probability::ZERO,
            });
        }
    }
    let normalization = log_sum(pairs.iter().map(|pair| pair.log_likelihood + pair.log_prior));
    let mut genotypes = [0.0; 10];
    for pair in &mut pairs {
        pair.posterior = probability((pair.log_likelihood + pair.log_prior - normalization).exp())?;
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
        fragments: included.len(),
        focal_fragments,
        focal_observations,
        focal_read_observations,
        focal_overlapping_fragments,
        focal_discordant_fragments,
        omitted_sites: offered - used,
        cpg_context_limited: limited,
        cpg_omitted_sites: bounded_sites - used,
        cpg_omitted_indels: usize::from(offered_indel && selected_indel.is_none()),
        offered_linking_fragments,
        realigned_fragments: realignment
            .as_ref()
            .map_or(0, |prepared| prepared.realigned_fragments()),
        aligned_fallback_fragments: realignment
            .as_ref()
            .map_or(0, |prepared| prepared.fallback_fragments()),
        omitted_indels: proposed_indels.len() - usize::from(selected_indel.is_some()),
    })
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
