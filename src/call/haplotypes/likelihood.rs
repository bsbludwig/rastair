use super::{Fragment, Haplotype, Reference, Site};
use crate::metrics::methylation::CpgSide;
use color_eyre::eyre::{Result, ensure};
use seqair_types::{Base, Probability, Strand};

#[derive(Debug, Clone, clap::Args, serde::Serialize)]
pub struct LikelihoodParams {
    /// Conversion probability for methylated cytosines in TAPS libraries
    #[arg(long, default_value = "0.95")]
    pub conversion_rate: Probability,
    /// Conversion probability for unmethylated cytosines in TAPS libraries
    #[arg(long, default_value = "0.003")]
    pub false_conversion_rate: Probability,
    /// Score ordinary sequencing reads without TAPS conversion
    #[arg(long)]
    pub non_taps: bool,
    /// Prior probability of a non-reference genotype at each candidate site
    #[arg(long, default_value = "0.001")]
    pub snp_prior: Probability,
    /// Gauss-Legendre nodes for integrating a uniform CpG methylation prior
    #[arg(long, default_value_t = 16, value_parser = clap::value_parser!(u16).range(2..=64))]
    pub quadrature_order: u16,
}

impl LikelihoodParams {
    pub fn validate(&self) -> Result<()> {
        ensure!(
            self.non_taps || self.conversion_rate > self.false_conversion_rate,
            "Methylated conversion must exceed unmethylated conversion"
        );
        ensure!(
            *self.snp_prior > 0.0 && *self.snp_prior < 1.0,
            "SNP prior must be strictly between zero and one"
        );
        ensure!((2..=64).contains(&self.quadrature_order), "Quadrature order must be 2..=64");
        Ok(())
    }

    fn conversion(&self, methylation: Option<f64>) -> f64 {
        if self.non_taps {
            0.0
        } else {
            let methylation = methylation.unwrap_or(0.0);
            methylation * *self.conversion_rate + (1.0 - methylation) * *self.false_conversion_rate
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Node {
    pub methylation: f64,
    pub log_weight: f64,
}

pub(super) fn quadrature(order: u16) -> Result<Vec<Node>> {
    ensure!((2..=64).contains(&order), "Quadrature order must be 2..=64");
    let n = usize::from(order);
    let nf = f64::from(order);
    let mut nodes = Vec::with_capacity(n);
    for i in 0..n.div_ceil(2) {
        let index = f64::from(u16::try_from(i)?);
        let mut root = (std::f64::consts::PI * (index + 0.75) / (nf + 0.5)).cos();
        let mut derivative = 0.0;
        let mut converged = false;
        for _ in 0..64 {
            let (p, previous) = legendre(order, root);
            derivative = nf * (root * p - previous) / (root * root - 1.0);
            let next = root - p / derivative;
            if (next - root).abs() < 1e-14 {
                root = next;
                let (p, previous) = legendre(order, root);
                derivative = nf * (root * p - previous) / (root * root - 1.0);
                converged = true;
                break;
            }
            root = next;
        }
        ensure!(converged, "Gauss-Legendre root did not converge");
        let weight = 1.0 / ((1.0 - root * root) * derivative * derivative);
        nodes.push(Node { methylation: (1.0 - root) / 2.0, log_weight: weight.ln() });
        if nodes.len() < n && !(n % 2 == 1 && i == n / 2) {
            nodes.push(Node { methylation: (1.0 + root) / 2.0, log_weight: weight.ln() });
        }
    }
    nodes.sort_by(|a, b| a.methylation.total_cmp(&b.methylation));
    ensure!(nodes.len() == n, "Wrong quadrature node count");
    Ok(nodes)
}

fn legendre(order: u16, x: f64) -> (f64, f64) {
    let (mut previous, mut current) = (1.0, x);
    for k in 2..=order {
        let k = f64::from(k);
        let next = ((2.0 * k - 1.0) * x * current - (k - 1.0) * previous) / k;
        (previous, current) = (current, next);
    }
    (current, previous)
}

pub(super) fn log_add(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let upper = a.max(b);
    upper + (-(a - b).abs()).exp().ln_1p()
}

pub(super) fn log_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values.into_iter().fold(f64::NEG_INFINITY, log_add)
}

fn sequencing_emission(observed: Base, template: Base, quality: Option<u8>) -> f64 {
    if observed.known_index().is_none() {
        return 1.0;
    }
    // Q0 and missing qualities have no discriminating power; treating Q0 as
    // certainty of a mismatch would turn uninformative reads into alt evidence.
    let error = quality.map_or(0.75, |q| 10.0_f64.powf(-f64::from(q) / 10.0).min(0.75));
    if observed == template { 1.0 - error } else { error / 3.0 }
}

pub(super) fn cpg_anchor(reference: &Reference, haplotype: &Haplotype, pos: u32) -> Option<u32> {
    let base = haplotype.base_at(reference, pos)?;
    match base {
        Base::C if haplotype.base_at(reference, pos.checked_add(1)?) == Some(Base::G) => Some(pos),
        Base::G if haplotype.base_at(reference, pos.checked_sub(1)?) == Some(Base::C) => {
            pos.checked_sub(1)
        }
        _ => None,
    }
}

pub(super) fn fragment_likelihood(
    fragment: &Fragment,
    haplotype: &Haplotype,
    reference: &Reference,
    loci: &[u32],
    methylation: f64,
    params: &LikelihoodParams,
) -> Result<f64> {
    ensure!(fragment.reads.len() <= 2, "A likelihood fragment may contain at most two reads");
    let known = fragment
        .reads
        .iter()
        .map(|read| read.strand)
        .find(|strand| *strand != Strand::Unknown)
        .filter(|strand| {
            fragment
                .reads
                .iter()
                .all(|read| read.strand == Strand::Unknown || read.strand == *strand)
        });
    // Mates have one molecular orientation. A known mate resolves an unknown
    // one; disagreeing assignments instead leave that orientation uncertain.
    let orientations: &[Strand] = match known {
        Some(Strand::OT) => &[Strand::OT],
        Some(Strand::OB) => &[Strand::OB],
        _ => &[Strand::OT, Strand::OB],
    };
    let mut result = f64::NEG_INFINITY;
    for &orientation in orientations {
        // Mapping failure is one event per alignment, rather than an
        // independent chance for every observed base of that alignment.
        for mask in 0..(1usize << fragment.reads.len()) {
            let mut probability = if orientations.len() == 2 { 0.5 } else { 1.0 };
            for (read_index, read) in fragment.reads.iter().enumerate() {
                let bad = mask & (1 << read_index) != 0;
                let mismapped = read.mapq.map_or(1.0, |q| 10.0_f64.powf(-f64::from(q) / 10.0));
                probability *= if bad { mismapped } else { 1.0 - mismapped };
            }
            if probability == 0.0 {
                continue;
            }
            let mut likelihood = probability.ln();
            for &pos in loci {
                let Some(base) = haplotype.base_at(reference, pos) else { continue };
                if base.known_index().is_none() {
                    continue;
                }
                let side =
                    [CpgSide::C, CpgSide::G].into_iter().find(|side| side.unmod_base() == base);
                let rate =
                    params.conversion(cpg_anchor(reference, haplotype, pos).map(|_| methylation));
                let (mut unconverted, mut converted) = (0.0, 0.0);
                for (read_index, read) in fragment.reads.iter().enumerate() {
                    let bad = mask & (1 << read_index) != 0;
                    let strand = orientation;
                    if let Ok(index) = read.observations.binary_search_by_key(&pos, |obs| obs.pos)
                        && let Some(obs) = read.observations.get(index)
                    {
                        let ordinary = if bad {
                            if obs.base.known_index().is_some() { 0.25 } else { 1.0 }
                        } else {
                            sequencing_emission(obs.base, base, obs.quality)
                        };
                        let modified =
                            if let Some(side) = side.filter(|side| side.strand() == strand) {
                                if bad {
                                    ordinary
                                } else {
                                    sequencing_emission(obs.base, side.mod_base(), obs.quality)
                                }
                            } else {
                                ordinary
                            };
                        unconverted += ordinary.ln();
                        converted += modified.ln();
                    }
                }
                // Overlapping mates share a conversion event at this position.
                likelihood += log_add((1.0 - rate).ln() + unconverted, rate.ln() + converted);
            }
            result = log_add(result, likelihood);
        }
    }
    ensure!(result.is_finite(), "Fragment likelihood is not finite");
    Ok(result)
}

pub(super) fn pair_prior(
    first: &Haplotype,
    second: &Haplotype,
    sites: &[Site],
    prior: f64,
) -> Result<f64> {
    let mut result = 0.0;
    let mut heterozygotes = 0u32;
    for site in sites {
        let a = first
            .allele_at(site.pos)
            .ok_or_else(|| color_eyre::eyre::eyre!("Missing haplotype allele"))?;
        let b = second
            .allele_at(site.pos)
            .ok_or_else(|| color_eyre::eyre::eyre!("Missing haplotype allele"))?;
        let genotypes = site.alleles.len() * (site.alleles.len() + 1) / 2;
        let probability = if a == site.reference && b == site.reference {
            1.0 - prior
        } else {
            prior / f64::from(u16::try_from(genotypes - 1)?)
        };
        result += probability.ln();
        heterozygotes += u32::from(a != b);
    }
    // The genotype prior is shared equally by its possible phase assignments.
    // Otherwise adding neighboring hets would change the focal SNP's prior.
    result -= f64::from(heterozygotes.saturating_sub(1)) * std::f64::consts::LN_2;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quadrature_integrates_uniform_polynomials() -> Result<()> {
        for order in [2, 3, 8, 16, 32, 64] {
            let nodes = quadrature(order)?;
            for power in 0..u32::from(order).min(12) {
                let integral: f64 = nodes
                    .iter()
                    .map(|node| node.log_weight.exp() * node.methylation.powf(f64::from(power)))
                    .sum();
                assert!((integral - 1.0 / f64::from(power + 1)).abs() < 1e-12);
            }
        }
        Ok(())
    }

    #[test]
    fn log_sum_survives_underflow_and_empty_input() {
        assert_eq!(log_sum([]), f64::NEG_INFINITY);
        assert!((log_add(-10_000.0, -10_000.0) + 10_000.0 - std::f64::consts::LN_2).abs() < 1e-12);
    }
}
