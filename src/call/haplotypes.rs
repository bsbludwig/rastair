use color_eyre::eyre::{OptionExt as _, Result, ensure};
use seqair_types::{Base, Probability, Strand};
use std::collections::BTreeSet;

mod likelihood;
pub use likelihood::LikelihoodParams;
mod diagnostics;
pub use diagnostics::{DiagnosticReadParams, HaplotypeDiagnosticParams, diagnose};

use likelihood::{Node, cpg_anchor, fragment_likelihood, log_add, log_sum, pair_prior, quadrature};

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
}

#[derive(Debug, Clone)]
pub struct Fragment {
    pub reads: Vec<ReadEvidence>,
}

#[derive(Debug, Clone)]
pub struct Haplotype {
    pub alleles: Vec<(u32, Base)>,
}

impl Haplotype {
    pub fn allele_at(&self, pos: u32) -> Option<Base> {
        self.alleles.iter().find(|(p, _)| *p == pos).map(|(_, base)| *base)
    }

    pub fn base_at(&self, reference: &Reference, pos: u32) -> Option<Base> {
        self.allele_at(pos).or_else(|| reference.base_at(pos))
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
    pub haplotypes: Vec<Haplotype>,
    pub pairs: Vec<PairEvidence>,
    pub sites: Vec<Site>,
    pub fragments: usize,
    pub focal_fragments: usize,
    pub focal_observations: usize,
    /// Other candidates were omitted to keep the Cartesian search bounded.
    pub omitted_sites: usize,
    /// This aligned prototype integrates at most one CpG rate per haplotype.
    /// It reduces context rather than tying unrelated CpG rates together.
    pub cpg_context_limited: bool,
}

impl SnpEvidence {
    pub fn allele(&self, base: Base) -> Option<&AlleleEvidence> {
        self.alleles.iter().find(|allele| allele.base == base)
    }
}

/// Alignment-conditioned SNP inference. Every focal A/C/G/T hypothesis is
/// retained, even when no read directly shows that genetic base. Methylation
/// is marginalized under independent uniform priors for the CpGs on the two
/// haplotypes; the same rate is shared by all fragments observing that CpG.
///
/// Context is reduced until each haplotype has at most one observed CpG. This
/// makes the prototype's numerical integral at most two-dimensional without
/// silently sharing a rate between different CpGs. It is not realignment.
pub fn infer_snp(
    reference: &Reference,
    focal: u32,
    neighbors: &[Site],
    fragments: &[Fragment],
    params: &LikelihoodParams,
    max_sites: usize,
) -> Result<SnpEvidence> {
    params.validate()?;
    ensure!(
        (1..=3).contains(&max_sites),
        "Aligned prototype supports one to three candidate sites"
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
    let mut limited = false;
    let (sites, haplotypes, loci) = loop {
        let mut sites = candidates.clone();
        sites.push(focal_site.clone());
        sites.sort_by_key(|site| site.pos);
        let haplotypes = enumerate(&sites);
        let loci = scored_loci(reference, &sites);
        let multiple_cpgs = !params.non_taps
            && haplotypes.iter().any(|haplotype| {
                loci.iter()
                    .filter_map(|&pos| cpg_anchor(reference, haplotype, pos))
                    .collect::<BTreeSet<_>>()
                    .len()
                    > 1
            });
        if !multiple_cpgs {
            break (sites, haplotypes, loci);
        }
        ensure!(candidates.pop().is_some(), "A single SNP has multiple independent observed CpGs");
        limited = true;
    };
    let included: Vec<_> = fragments
        .iter()
        .filter(|fragment| {
            fragment.reads.iter().any(|read| {
                read.observations.iter().any(|obs| loci.binary_search(&obs.pos).is_ok())
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
    let nodes = quadrature(params.quadrature_order)?;
    let mut matrix = Vec::with_capacity(haplotypes.len());
    for haplotype in &haplotypes {
        let methylable = !params.non_taps
            && loci.iter().any(|&pos| cpg_anchor(reference, haplotype, pos).is_some());
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
                    .map(|fragment| {
                        fragment_likelihood(
                            fragment,
                            haplotype,
                            reference,
                            &loci,
                            node.methylation,
                            params,
                        )
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
                log_prior: pair_prior(a, b, &sites, *params.snp_prior)?,
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
    let focal_observations = fragments
        .iter()
        .flat_map(|fragment| &fragment.reads)
        .flat_map(|read| &read.observations)
        .filter(|obs| obs.pos == focal)
        .count();
    let focal_fragments = fragments
        .iter()
        .filter(|fragment| {
            fragment.reads.iter().any(|read| read.observations.iter().any(|obs| obs.pos == focal))
        })
        .count();
    let used = sites.len() - 1;
    Ok(SnpEvidence {
        pos: focal,
        reference: ref_base,
        alleles,
        genotype_entropy: entropy,
        haplotypes,
        pairs,
        sites,
        fragments: included.len(),
        focal_fragments,
        focal_observations,
        omitted_sites: offered - used,
        cpg_context_limited: limited,
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
    let mut haplotypes = vec![Haplotype { alleles: Vec::new() }];
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
