use super::{Fragment, Haplotype, Reference, Site};
use crate::metrics::methylation::CpgSide;
use color_eyre::eyre::{Result, ensure};
use seqair_types::{Base, Probability, Strand};

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum, serde::Serialize)]
pub enum LikelihoodEngine {
    Aligned,
    PairHmm,
}

#[derive(Debug, Clone, clap::Args, serde::Serialize)]
pub struct LikelihoodParams {
    /// Read likelihood engine; `PairHMM` sums over local alignments
    #[arg(long, value_enum, default_value = "aligned")]
    pub likelihood_engine: LikelihoodEngine,
    /// Strength of the symmetric Beta prior for estimated CpG methylation rates
    #[arg(long, default_value_t = 4.0)]
    pub methylation_prior_strength: f64,
    /// Minimum pilot probability that a molecule carries a CpG before estimating its rate
    #[arg(long, default_value = "0.9")]
    pub methylation_min_cpg_probability: Probability,
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
    /// Gauss-Legendre nodes for estimating each shared CpG methylation rate
    #[arg(long, default_value_t = 16, value_parser = clap::value_parser!(u16).range(2..=64))]
    pub quadrature_order: u16,
    /// Probability of opening each insertion or deletion in the `PairHMM`
    #[arg(long, default_value = "0.00003162277660168379")]
    pub gap_open: Probability,
    /// Probability of extending an insertion or deletion in the `PairHMM`
    #[arg(long, default_value = "0.1")]
    pub gap_extend: Probability,
    /// Prior probability of a non-reference genotype at the nuisance indel site
    #[arg(long, default_value = "0.0001")]
    pub indel_prior: Probability,
    /// Local realignment flank around selected SNPs and the nuisance indel
    #[arg(long, default_value_t = 12, value_parser = clap::value_parser!(u16).range(2..=100))]
    pub realignment_flank: u16,
    /// Maximum length of a CIGAR insertion/deletion proposed as a nuisance allele
    #[arg(long, default_value_t = 20, value_parser = clap::value_parser!(u16).range(1..=100))]
    pub max_indel_length: u16,
}

impl LikelihoodParams {
    pub fn validate(&self) -> Result<()> {
        ensure!(
            self.methylation_prior_strength.is_finite() && self.methylation_prior_strength >= 2.0,
            "Methylation prior strength must be finite and at least two"
        );
        ensure!(
            self.non_taps || self.conversion_rate > self.false_conversion_rate,
            "Methylated conversion must exceed unmethylated conversion"
        );
        ensure!(
            *self.snp_prior > 0.0 && *self.snp_prior < 1.0,
            "SNP prior must be strictly between zero and one"
        );
        ensure!((2..=64).contains(&self.quadrature_order), "Quadrature order must be 2..=64");
        ensure!(
            *self.gap_open > 0.0 && *self.gap_open < 0.5,
            "Gap open must be between zero and 0.5"
        );
        ensure!(
            *self.gap_extend > 0.0 && *self.gap_extend < 1.0,
            "Gap extension must be between zero and one"
        );
        ensure!(
            *self.indel_prior > 0.0 && *self.indel_prior < 1.0,
            "Indel prior must be between zero and one"
        );
        ensure!((2..=100).contains(&self.realignment_flank), "Realignment flank must be 2..=100");
        ensure!((1..=100).contains(&self.max_indel_length), "Maximum indel length must be 1..=100");
        Ok(())
    }

    pub(super) fn conversion(&self, methylation: Option<f64>) -> f64 {
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

pub(super) fn sequencing_emission(observed: Base, template: Base, quality: Option<u8>) -> f64 {
    if observed.known_index().is_none() {
        return 1.0;
    }
    // Q0 and missing qualities have no discriminating power; treating Q0 as
    // certainty of a mismatch would turn uninformative reads into alt evidence.
    let error = quality.map_or(0.75, |q| 10.0_f64.powf(-f64::from(q) / 10.0).min(0.75));
    if observed == template { 1.0 - error } else { error / 3.0 }
}

#[derive(Debug, Clone, Copy)]
pub(super) struct BaseMeasurement {
    pub base: Base,
    pub quality: Option<u8>,
    pub mapq: Option<u8>,
    pub read_index: usize,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct MolecularBase {
    pub base: Base,
    pub quality: Option<u8>,
    pub source: Option<usize>,
    pub measurements: usize,
    pub discordant: bool,
}

pub(super) fn consensus(
    measurements: impl IntoIterator<Item = BaseMeasurement>,
) -> Result<MolecularBase> {
    let mut known = [None; 2];
    let mut known_count = 0;
    let mut count = 0;
    for measurement in measurements {
        count += 1;
        ensure!(
            count <= 2 && measurement.read_index < 2,
            "A molecular base requires at most two mate observations"
        );
        ensure!(measurement.quality != Some(255), "Missing base quality must be None");
        if measurement.base.known_index().is_some() {
            *known
                .get_mut(known_count)
                .ok_or_else(|| color_eyre::eyre::eyre!("Too many known mate observations"))? =
                Some(measurement);
            known_count += 1;
        }
    }
    let mut result = MolecularBase {
        base: Base::Unknown,
        quality: None,
        source: None,
        measurements: known_count,
        discordant: false,
    };
    let rank = |m: &BaseMeasurement| {
        (m.quality.unwrap_or(0), m.mapq.unwrap_or(0), std::cmp::Reverse(m.read_index))
    };
    if let Some(best) = known.iter().flatten().max_by_key(|measurement| rank(measurement)) {
        result.base = best.base;
        result.quality = best.quality;
        result.source = Some(best.read_index);
        if let Some(other) =
            known.iter().flatten().find(|m| m.read_index != best.read_index && m.base != best.base)
        {
            result.discordant = true;
            // Both calls observe one nucleotide. Agreement cannot add Phred
            // confidence, and a competing call spends that confidence instead.
            result.quality = best.quality.map(|q| q.saturating_sub(other.quality.unwrap_or(0)));
        }
    }
    Ok(result)
}

pub(super) fn consensus_at(fragment: &Fragment, pos: u32) -> Result<MolecularBase> {
    consensus(fragment.reads.iter().enumerate().filter_map(|(read_index, read)| {
        let index = read.observations.binary_search_by_key(&pos, |obs| obs.pos).ok()?;
        let obs = read.observations.get(index)?;
        Some(BaseMeasurement { base: obs.base, quality: obs.quality, mapq: read.mapq, read_index })
    }))
}

pub(super) fn cpg_anchor(reference: &Reference, haplotype: &Haplotype, pos: u32) -> Option<u32> {
    let base = haplotype.base_at(reference, pos)?;
    match base {
        Base::C
            if haplotype
                .adjacent_base(reference, pos, true)
                .is_some_and(|(_, base)| base == Base::G) =>
        {
            Some(pos)
        }
        Base::G => haplotype
            .adjacent_base(reference, pos, false)
            .filter(|(_, base)| *base == Base::C)
            .map(|(pos, _)| pos),
        _ => None,
    }
}

#[cfg(test)]
pub(super) fn fragment_likelihood(
    fragment: &Fragment,
    haplotype: &Haplotype,
    reference: &Reference,
    loci: &[u32],
    methylation: f64,
    params: &LikelihoodParams,
) -> Result<f64> {
    fragment_likelihood_with_model(
        fragment,
        haplotype,
        reference,
        loci,
        params,
        &super::methylation::Scoring::Fixed(methylation),
        &BTreeSet::new(),
    )
}

use std::collections::BTreeSet;

pub(super) fn fragment_likelihood_with_model(
    fragment: &Fragment,
    haplotype: &Haplotype,
    reference: &Reference,
    loci: &[u32],
    params: &LikelihoodParams,
    scoring: &super::methylation::Scoring<'_>,
    cpgs: &BTreeSet<super::methylation::Cpg>,
) -> Result<f64> {
    ensure!(fragment.reads.len() <= 2, "A likelihood fragment may contain at most two reads");
    let bases = loci
        .iter()
        .map(|&pos| Ok((pos, consensus_at(fragment, pos)?)))
        .collect::<Result<Vec<_>>>()?;
    let owns = |read_index| bases.iter().any(|(_, base)| base.source == Some(read_index));
    let orientations = orientations_from(
        fragment
            .reads
            .iter()
            .filter(|read| {
                read.observations.iter().any(|obs| {
                    obs.base.known_index().is_some() && loci.binary_search(&obs.pos).is_ok()
                })
            })
            .map(|read| read.strand),
    );
    let mut result = f64::NEG_INFINITY;
    for &orientation in orientations {
        // Mapping failure is one event per alignment, rather than an
        // independent chance for every observed base of that alignment.
        for mask in 0..(1usize << fragment.reads.len()) {
            if fragment
                .reads
                .iter()
                .enumerate()
                .any(|(index, _)| !owns(index) && mask & (1 << index) != 0)
            {
                continue;
            }
            let mut probability = if orientations.len() == 2 { 0.5 } else { 1.0 };
            for (read_index, read) in fragment.reads.iter().enumerate() {
                if !owns(read_index) {
                    continue;
                }
                let bad = mask & (1 << read_index) != 0;
                let mismapped = read.mapq.map_or(1.0, |q| 10.0_f64.powf(-f64::from(q) / 10.0));
                probability *= if bad { mismapped } else { 1.0 - mismapped };
            }
            if probability == 0.0 {
                continue;
            }
            let mut likelihood = probability.ln();
            for &(pos, observed) in &bases {
                let Some(source) = observed.source else { continue };
                let Some(base) = haplotype.base_at(reference, pos) else { continue };
                if base.known_index().is_none() {
                    continue;
                }
                let position = super::SequencePosition::reference(pos);
                if !params.non_taps
                    && let Some((ordinary, modified)) =
                        scoring.neutral_bases(position, orientation, base)
                {
                    likelihood += if mask & (1 << source) != 0 {
                        0.25_f64.ln()
                    } else {
                        (0.5 * sequencing_emission(observed.base, ordinary, observed.quality)
                            + 0.5 * sequencing_emission(observed.base, modified, observed.quality))
                        .ln()
                    };
                    continue;
                }
                let side =
                    [CpgSide::C, CpgSide::G].into_iter().find(|side| side.unmod_base() == base);
                let cpg = cpgs.iter().find(|cpg| cpg.contains(position));
                let methylation = if cpg_anchor(reference, haplotype, pos).is_some() {
                    Some(scoring.rate(cpg))
                } else {
                    None
                };
                let rate = params.conversion(methylation);
                let bad = mask & (1 << source) != 0;
                let unconverted = if bad {
                    0.25
                } else {
                    sequencing_emission(observed.base, base, observed.quality)
                }
                .ln();
                let converted = if let Some(side) = side.filter(|side| side.strand() == orientation)
                {
                    if bad {
                        0.25
                    } else {
                        sequencing_emission(observed.base, side.mod_base(), observed.quality)
                    }
                    .ln()
                } else {
                    unconverted
                };
                likelihood += log_add((1.0 - rate).ln() + unconverted, rate.ln() + converted);
            }
            result = log_add(result, likelihood);
        }
    }
    ensure!(result.is_finite(), "Fragment likelihood is not finite");
    Ok(result)
}

pub(super) fn orientations_from(strands: impl IntoIterator<Item = Strand>) -> &'static [Strand] {
    let mut known = None;
    for strand in strands.into_iter().filter(|strand| *strand != Strand::Unknown) {
        if known.is_some_and(|known| known != strand) {
            return &[Strand::OT, Strand::OB];
        }
        known = Some(strand);
    }
    // Mates have one molecular orientation. A known mate resolves an unknown
    // one; disagreeing assignments instead leave that orientation uncertain.
    match known {
        Some(Strand::OT) => &[Strand::OT],
        Some(Strand::OB) => &[Strand::OB],
        _ => &[Strand::OT, Strand::OB],
    }
}

pub(super) fn pair_prior(
    first: &Haplotype,
    second: &Haplotype,
    sites: &[Site],
    prior: f64,
    indel_prior: f64,
    indel_alternatives: usize,
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
    if indel_alternatives > 0 {
        let alleles = indel_alternatives + 1;
        let alternate_genotypes = alleles * (alleles + 1) / 2 - 1;
        result += if first.indel.is_none() && second.indel.is_none() {
            (1.0 - indel_prior).ln()
        } else {
            (indel_prior / f64::from(u32::try_from(alternate_genotypes)?)).ln()
        };
        heterozygotes += u32::from(first.indel != second.indel);
    }
    // The genotype prior is shared equally by its possible phase assignments.
    // Otherwise adding neighboring hets would change the focal SNP's prior.
    result -= f64::from(heterozygotes.saturating_sub(1)) * std::f64::consts::LN_2;
    Ok(result)
}

pub(super) fn retained_pair_priors(
    haplotypes: &[Haplotype],
    sites: &[Site],
    params: &LikelihoodParams,
    has_indel: bool,
) -> Result<Vec<f64>> {
    let indels = haplotypes
        .iter()
        .filter_map(|haplotype| haplotype.indel.as_ref())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let mut counts = std::collections::BTreeMap::<Vec<(u8, u8)>, usize>::new();
    let mut pairs = Vec::new();
    for (index, first) in haplotypes.iter().enumerate() {
        for second in haplotypes.iter().skip(index) {
            let mut key = Vec::with_capacity(sites.len() + usize::from(has_indel));
            let mut heterozygotes = 0u32;
            for site in sites {
                let a = first
                    .allele_at(site.pos)
                    .ok_or_else(|| color_eyre::eyre::eyre!("Missing retained allele"))?
                    as u8;
                let b = second
                    .allele_at(site.pos)
                    .ok_or_else(|| color_eyre::eyre::eyre!("Missing retained allele"))?
                    as u8;
                key.push((a.min(b), a.max(b)));
                heterozygotes += u32::from(a != b);
            }
            if has_indel {
                let allele = |indel: Option<&super::Indel>| -> Result<u8> {
                    match indel {
                        None => Ok(0),
                        Some(indel) => Ok(u8::try_from(
                            indels.binary_search(&indel).map_err(|_| {
                                color_eyre::eyre::eyre!("Missing candidate indel prior")
                            })? + 1,
                        )?),
                    }
                };
                let a = allele(first.indel.as_ref())?;
                let b = allele(second.indel.as_ref())?;
                key.push((a.min(b), a.max(b)));
                heterozygotes += u32::from(a != b);
            }
            *counts.entry(key.clone()).or_default() += 1;
            pairs.push((
                key,
                heterozygotes.saturating_sub(1),
                pair_prior(
                    first,
                    second,
                    sites,
                    *params.snp_prior,
                    *params.indel_prior,
                    if has_indel { indels.len() } else { 0 },
                )?,
            ));
        }
    }
    pairs
        .into_iter()
        .map(|(key, phase_bits, prior)| {
            let count = *counts
                .get(&key)
                .ok_or_else(|| color_eyre::eyre::eyre!("Missing retained phase count"))?;
            if 1usize.checked_shl(phase_bits) == Some(count) {
                return Ok(prior);
            }
            // Pruning removes phase assignments. Divide a genotype's prior over its
            // surviving assignments, rather than penalizing it for the search budget.
            Ok(prior + f64::from(phase_bits) * std::f64::consts::LN_2 - (count as f64).ln())
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mate_consensus_never_adds_quality_and_penalizes_disagreement() -> Result<()> {
        for (
            first_base,
            second_base,
            first_quality,
            second_quality,
            expected_base,
            expected_quality,
        ) in [
            (Base::A, Base::A, 40, 40, Base::A, 40),
            (Base::A, Base::A, 20, 40, Base::A, 40),
            (Base::A, Base::G, 40, 30, Base::A, 10),
            (Base::A, Base::G, 20, 40, Base::G, 20),
            (Base::C, Base::T, 40, 40, Base::C, 0),
        ] {
            let measurements = [
                BaseMeasurement {
                    base: first_base,
                    quality: Some(first_quality),
                    mapq: Some(60),
                    read_index: 0,
                },
                BaseMeasurement {
                    base: second_base,
                    quality: Some(second_quality),
                    mapq: Some(60),
                    read_index: 1,
                },
            ];
            let base = consensus(measurements)?;
            assert_eq!(base.base, expected_base);
            assert_eq!(base.quality, Some(expected_quality));
            assert_eq!(base.measurements, 2);
            assert_eq!(base.discordant, first_base != second_base);
            let reversed = consensus(measurements.into_iter().rev())?;
            assert_eq!(base.base, reversed.base);
            assert_eq!(base.quality, reversed.quality);
            let emission_sum = Base::KNOWN
                .into_iter()
                .map(|template| sequencing_emission(base.base, template, base.quality))
                .sum::<f64>();
            assert!((emission_sum - 1.0).abs() < 1e-12);
            if expected_quality == 0 {
                assert!(Base::KNOWN.into_iter().all(|template| sequencing_emission(
                    base.base,
                    template,
                    base.quality
                ) == 0.25));
            }
        }
        Ok(())
    }

    #[test]
    fn consensus_ignores_unknown_bases_and_missing_qualities_are_neutral() -> Result<()> {
        let known =
            BaseMeasurement { base: Base::G, quality: Some(30), mapq: Some(60), read_index: 1 };
        let unknown = BaseMeasurement {
            base: Base::Unknown,
            quality: Some(40),
            mapq: Some(60),
            read_index: 0,
        };
        let base = consensus([unknown, known])?;
        assert_eq!(base.base, Base::G);
        assert_eq!(base.quality, Some(30));
        assert_eq!(base.measurements, 1);
        let base = consensus([BaseMeasurement { quality: None, ..known }])?;
        assert!(Base::KNOWN.into_iter().all(|template| sequencing_emission(
            base.base,
            template,
            base.quality
        ) == 0.25));
        let base = consensus([unknown])?;
        assert_eq!(base.source, None);
        assert_eq!(base.measurements, 0);
        assert!(consensus([BaseMeasurement { quality: Some(255), ..known }]).is_err());
        assert!(consensus([BaseMeasurement { read_index: 2, ..known }]).is_err());
        assert!(consensus([known; 3]).is_err());
        Ok(())
    }

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
