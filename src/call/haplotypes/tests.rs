//! Tests of the likelihood model through the regional search: a focal site
//! proposed with every base, its neighbours with their alleles, and what the
//! diploid posterior says about each allele at the focal site.

use super::*;
use Base::*;
use likelihood::{cpg_anchor, fragment_likelihood};

pub(super) fn params(non_taps: bool) -> LikelihoodParams {
    LikelihoodParams {
        error_profile: super::ErrorProfile { repeats: Vec::new(), ..Default::default() },
        conversion_rate: Probability::new(0.95).expect("probability"),
        false_conversion_rate: Probability::new(0.003).expect("probability"),
        snp_prior: Probability::new(0.001).expect("probability"),
        non_taps,
        quadrature_order: 16,
        methylation_prior_strength: 4.0,
        methylation_min_cpg_probability: Probability::new(0.9).expect("probability"),
        gap_open: Probability::new(10.0_f64.powf(-4.5)).expect("probability"),
        gap_extend: Probability::new(0.1).expect("probability"),
        indel_prior: Probability::new(0.0001).expect("probability"),
        realignment_flank: 12,
        max_indel_length: 20,
    }
}

/// The window the aligned model calibrates on: the loci themselves.
pub(super) fn loci_window(loci: &[u32]) -> std::ops::Range<u32> {
    let first = loci.first().copied().expect("loci");
    let last = loci.last().copied().expect("loci");
    first..last + 1
}

/// Every combination of the sites' alleles.
pub(super) fn enumerate(sites: &[Site]) -> Vec<Haplotype> {
    let mut haplotypes = vec![Haplotype { alleles: Vec::new(), indels: Vec::new() }];
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

#[derive(Debug, Clone, Copy)]
pub(super) struct AlleleEvidence {
    pub base: Base,
    pub log_bayes_factor: f64,
    pub posterior: Probability,
    pub heterozygous: Probability,
    pub homozygous: Probability,
}

/// What the regional model says about one focal site.
#[derive(Debug)]
pub(super) struct SnpEvidence {
    pub alleles: Vec<AlleleEvidence>,
    pub haplotypes: Vec<Haplotype>,
    pub pairs: Vec<PairEvidence>,
    pub sites: Vec<Site>,
    /// Unphased A/C/G/T genotypes at the focal site in VCF triangular order.
    pub genotype_posteriors: [Probability; 10],
    pub realigned_fragments: usize,
    pub aligned_fallback_fragments: usize,
    pub focal_fragments: usize,
    pub focal_observations: usize,
    pub focal_read_observations: usize,
    pub focal_overlapping_fragments: usize,
    pub focal_discordant_fragments: usize,
    /// Fragments covering the focal site and another site of the region.
    pub linking_fragments: usize,
}

impl SnpEvidence {
    pub fn allele(&self, base: Base) -> Option<&AlleleEvidence> {
        self.alleles.iter().find(|allele| allele.base == base)
    }

    /// Conditional on the retained chromosomes.
    pub fn genotypes_at(&self, pos: u32) -> Result<[Probability; 10]> {
        ensure!(self.sites.iter().any(|site| site.pos == pos), "Position is outside the region");
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

pub(super) fn infer_snp(
    reference: &Reference,
    focal: u32,
    neighbors: &[Site],
    fragments: &[Fragment],
    params: &LikelihoodParams,
    max_sites: usize,
) -> Result<SnpEvidence> {
    infer_snp_with_indels(reference, focal, neighbors, &[], fragments, params, max_sites)
}

/// The regional model over the focal site, proposed with every base, and the
/// `max_sites - 1` nearest of `neighbors`; each of `indels` is an alternate
/// at the candidate anchoring it.
pub(super) fn infer_snp_with_indels(
    reference: &Reference,
    focal: u32,
    neighbors: &[Site],
    indels: &[Indel],
    fragments: &[Fragment],
    params: &LikelihoodParams,
    max_sites: usize,
) -> Result<SnpEvidence> {
    params.validate()?;
    let ref_base = reference.base_at(focal).ok_or_eyre("Focal site outside reference")?;
    let mut sites = vec![Site::new(focal, ref_base, Base::KNOWN)?];
    let mut nearest: Vec<_> = neighbors.iter().filter(|site| site.pos != focal).cloned().collect();
    nearest.sort_by_key(|site| (site.pos.abs_diff(focal), site.pos));
    for site in nearest.into_iter().take(max_sites.saturating_sub(1)) {
        ensure!(
            reference.base_at(site.pos) == Some(site.reference),
            "Candidate reference mismatch"
        );
        sites.push(site);
    }
    let mut candidates: Vec<_> =
        sites.into_iter().map(|site| regional::Candidate { site, indels: Vec::new() }).collect();
    for indel in indels {
        let indel = indel.normalized(reference)?;
        let anchor = indel.start.checked_sub(1).ok_or_eyre("Indel at the reference start")?;
        match candidates.iter_mut().find(|candidate| candidate.site.pos == anchor) {
            Some(candidate) => candidate.indels.push(indel),
            None => {
                let base =
                    reference.base_at(anchor).ok_or_eyre("Indel anchor outside reference")?;
                candidates.push(regional::Candidate {
                    site: Site { pos: anchor, reference: base, alleles: vec![base] },
                    indels: vec![indel],
                });
            }
        }
    }
    candidates.sort_by_key(|candidate| candidate.site.pos);
    let search = ExpansionParams { max_sites: 64, ..ExpansionParams::default() };
    let mut scorer = Scorer::new(params)?;
    let region = regional::infer(reference, &candidates, fragments, params, &search, &mut scorer)?;
    ensure!(region.used == candidates.len(), "The region did not reach every candidate");

    let mut genotypes = [0.0; 10];
    for pair in &region.pairs {
        let (a, b) = pair_bases(pair, &region.haplotypes, focal)?;
        let index = genotype_index(a, b).ok_or_eyre("Unknown focal genotype")?;
        *genotypes.get_mut(index).ok_or_eyre("Invalid genotype index")? += *pair.posterior;
    }
    let mut genotype_posteriors = [Probability::ZERO; 10];
    for (posterior, value) in genotype_posteriors.iter_mut().zip(genotypes) {
        *posterior = probability(value)?;
    }
    let mut alleles = Vec::with_capacity(3);
    for base in Base::KNOWN.into_iter().filter(|base| *base != ref_base) {
        let (mut present_likelihood, mut absent_likelihood, mut present_prior, mut absent_prior) =
            (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        let (mut present, mut heterozygous, mut homozygous) = (0.0, 0.0, 0.0);
        for pair in &region.pairs {
            let (a, b) = pair_bases(pair, &region.haplotypes, focal)?;
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
        });
    }
    let covers = |fragment: &Fragment, pos| fragment.reads.iter().any(|read| read.covers(pos));
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
    let focal_fragments = fragments.iter().filter(|fragment| covers(fragment, focal)).count();
    let linking_fragments = fragments
        .iter()
        .filter(|fragment| {
            covers(fragment, focal)
                && region.sites.iter().any(|site| site.pos != focal && covers(fragment, site.pos))
        })
        .count();
    Ok(SnpEvidence {
        alleles,
        aligned_fallback_fragments: region.fragments - region.realigned,
        haplotypes: region.haplotypes,
        pairs: region.pairs,
        sites: region.sites,
        genotype_posteriors,
        realigned_fragments: region.realigned,
        focal_fragments,
        focal_observations,
        focal_read_observations,
        focal_overlapping_fragments,
        focal_discordant_fragments,
        linking_fragments,
    })
}

#[test]
fn posterior_normalization_survives_large_common_log_offsets() -> Result<()> {
    let sites = [Site::new(1001, A, [G])?];
    let haplotypes = enumerate(&sites);
    let params = params(true);
    let baseline =
        posterior_pairs(&haplotypes, &sites, &[vec![-1.0; 100], vec![-1.05; 100]], &params, false)?;
    let shifted = posterior_pairs(
        &haplotypes,
        &sites,
        &[vec![-10_000_001.0; 100], vec![-10_000_001.05; 100]],
        &params,
        false,
    )?;
    assert!((shifted.iter().map(|pair| *pair.posterior).sum::<f64>() - 1.0).abs() < 1e-12);
    for (a, b) in baseline.iter().zip(&shifted) {
        assert!((*a.posterior - *b.posterior).abs() < 1e-6);
    }
    Ok(())
}

pub(super) fn reference(bases: &[u8]) -> Reference {
    Reference { start: 1000, bases: bases.iter().copied().map(Base::from).collect() }
}

fn read(strand: Strand, bases: &[(u32, Base)]) -> ReadEvidence {
    ReadEvidence {
        strand,
        mapq: Some(60),
        sequence: None,
        observations: bases
            .iter()
            .map(|&(pos, base)| Observation { pos, base, quality: Some(40) })
            .collect(),
    }
}

fn fragment(strand: Strand, bases: &[(u32, Base)]) -> Fragment {
    Fragment { reads: vec![read(strand, bases)] }
}

fn copies(count: usize, strand: Strand, bases: &[(u32, Base)]) -> Vec<Fragment> {
    std::iter::repeat_n(fragment(strand, bases), count).collect()
}

fn allele(evidence: &SnpEvidence, base: Base) -> &AlleleEvidence {
    evidence.allele(base).expect("alternate evidence")
}

fn probability_sum(evidence: &SnpEvidence) -> f64 {
    evidence.pairs.iter().map(|pair| *pair.posterior).sum()
}

/// Five `A/G` hets linked pairwise by fragments are phased as one region:
/// the best pair is all-`A` against all-`G`, and every site's genotype is
/// the het.
#[test]
fn a_chain_of_linked_sites_is_phased_jointly() -> Result<()> {
    let reference = Reference { start: 1000, bases: vec![A; 1000] };
    let positions = [1001, 1121, 1241, 1361, 1481];
    let neighbors =
        positions.iter().skip(1).map(|&pos| Site::new(pos, A, [G])).collect::<Result<Vec<_>>>()?;
    let mut reads = Vec::new();
    for pair in positions.windows(2) {
        let [a, b] = pair else { continue };
        reads.extend(copies(10, Strand::OT, &[(*a, A), (*b, A)]));
        reads.extend(copies(10, Strand::OB, &[(*a, G), (*b, G)]));
    }
    let result = infer_snp(&reference, 1001, &neighbors, &reads, &params(true), 5)?;
    assert_eq!(result.sites.len(), 5);
    assert!(*allele(&result, G).heterozygous > 0.999);
    assert!((probability_sum(&result) - 1.0).abs() < 1e-12);
    let best = result
        .pairs
        .iter()
        .max_by(|a, b| a.posterior.total_cmp(&b.posterior))
        .ok_or_eyre("best pair")?;
    let a = result.haplotypes.get(best.first).ok_or_eyre("first")?;
    let b = result.haplotypes.get(best.second).ok_or_eyre("second")?;
    assert!(
        positions.iter().all(|&pos| a.allele_at(pos) == Some(A))
            || positions.iter().all(|&pos| a.allele_at(pos) == Some(G))
    );
    assert!(positions.iter().all(|&pos| a.allele_at(pos) != b.allele_at(pos)));
    for &pos in &positions {
        let genotypes = result.genotypes_at(pos)?;
        assert!(
            *genotypes.get(genotype_index(A, G).ok_or_eyre("index")?).ok_or_eyre("A/G")?
                > Probability::new(0.999)?
        );
    }
    assert_eq!(result.genotypes_at(1001)?, result.genotype_posteriors);
    assert!(result.genotypes_at(1999).is_err());
    Ok(())
}

#[test]
fn retained_phase_assignments_share_one_genotype_prior_after_pruning() -> Result<()> {
    let sites = [Site::new(1001, A, [G])?, Site::new(1002, A, [G])?];
    let haplotypes = [
        Haplotype { alleles: vec![(1001, A), (1002, A)], indels: Vec::new() },
        Haplotype { alleles: vec![(1001, G), (1002, G)], indels: Vec::new() },
    ];
    let priors = likelihood::retained_pair_priors(&haplotypes, &sites, &params(true), false)?;
    assert!(
        (priors.get(1).ok_or_eyre("het prior")?.exp() - (0.001_f64 / 2.0).powi(2)).abs() < 1e-18
    );
    Ok(())
}

#[test]
fn cpg_partner_hypotheses_are_joint_when_a_variant_removes_conversion_context() -> Result<()> {
    let reference = reference(b"ACGAAAAA");
    let neighbors = [Site::new(1001, C, [T])?];
    let mut reads = copies(6, Strand::OT, &[(1001, T), (1002, G)]);
    reads.extend(copies(2, Strand::OT, &[(1001, C), (1002, G)]));
    reads.extend(copies(4, Strand::OT, &[(1001, T), (1002, A)]));
    let result = infer_snp(&reference, 1002, &neighbors, &reads, &params(false), 2)?;
    assert_eq!(result.sites.len(), 2, "the conversion-ambiguous partner is part of the context");
    let genotypes = result.genotypes_at(1001)?;
    assert!(
        *genotypes.get(genotype_index(C, T).ok_or_eyre("C/T index")?).ok_or_eyre("C/T")?
            > Probability::new(0.99)?
    );
    assert!(*allele(&result, A).heterozygous > 0.99);
    Ok(())
}

#[test]
fn a_partner_variant_does_not_create_a_false_cpg_at_an_uncovered_focal_base() -> Result<()> {
    let reference = reference(b"ATGAAACA");
    let neighbors = [Site::new(1002, G, [A])?];
    let reads = copies(20, Strand::OB, &[(1002, A)]);
    let result = infer_snp(&reference, 1001, &neighbors, &reads, &params(false), 2)?;
    assert_eq!(result.focal_observations, 0);
    assert_eq!(result.sites.len(), 2);
    assert!(
        *allele(&result, C).posterior < 0.01,
        "the partner can be a genetic A, rather than conversion implying a de-novo C"
    );
    Ok(())
}

#[test]
fn an_isolated_snp_is_scored_without_a_called_genotype() -> Result<()> {
    let mut reads = copies(12, Strand::OT, &[(1001, A)]);
    reads.extend(copies(12, Strand::OB, &[(1001, G)]));
    let evidence = infer_snp(&reference(b"AAA"), 1001, &[], &reads, &params(true), 3)?;
    assert!(*allele(&evidence, G).heterozygous > 0.999);
    assert!(*allele(&evidence, G).homozygous < 1e-5);
    assert_eq!(evidence.haplotypes.len(), 4);
    assert_eq!(evidence.pairs.len(), 10);
    assert_eq!(evidence.focal_fragments, 24);
    assert!((probability_sum(&evidence) - 1.0).abs() < 1e-12);
    Ok(())
}

#[test]
fn zero_depth_returns_the_prior_and_zero_bayes_factors() -> Result<()> {
    let evidence = infer_snp(&reference(b"AAA"), 1001, &[], &[], &params(true), 3)?;
    for alt in &evidence.alleles {
        assert!((*alt.posterior - 0.001 * 4.0 / 9.0).abs() < 1e-12);
        assert!(alt.log_bayes_factor.abs() < 1e-12);
    }
    assert_eq!(evidence.focal_fragments, 0);
    Ok(())
}

#[test]
fn neutral_conversion_preserves_other_alleles_and_the_opposite_strand() -> Result<()> {
    for (sequence, pos, ordinary, converted, confounded, opposite, other) in [
        (b"ACGA".as_slice(), 1001, C, T, Strand::OT, Strand::OB, A),
        (b"ACGA".as_slice(), 1002, G, A, Strand::OB, Strand::OT, T),
        (b"ATGA".as_slice(), 1001, C, T, Strand::OT, Strand::OB, A),
    ] {
        let reference = reference(sequence);
        let sites = [Site::new(pos, reference.base_at(pos).ok_or_eyre("base")?, Base::KNOWN)?];
        let haplotypes = enumerate(&sites);
        let loci = scored_loci(&reference, &sites);
        let params = params(false);
        let context =
            methylation::Context::new(&reference, &haplotypes, &loci, loci_window(&loci))?;
        let score = |template, observed, strand| -> Result<f64> {
            let (index, haplotype) = haplotypes
                .iter()
                .enumerate()
                .find(|(_, haplotype)| haplotype.allele_at(pos) == Some(template))
                .ok_or_eyre("haplotype")?;
            likelihood::fragment_likelihood_with_model(
                &fragment(strand, &[(pos, observed)]),
                haplotype,
                &reference,
                &loci,
                &params,
                &methylation::Scoring::Neutral(&context),
                context.haplotypes.get(index).ok_or_eyre("context")?,
            )
        };
        for observed in Base::KNOWN {
            assert!(
                (score(ordinary, observed, confounded)? - score(converted, observed, confounded)?)
                    .abs()
                    < 1e-12
            );
        }
        assert!(score(other, other, confounded)? > score(ordinary, other, confounded)? + 8.0);
        assert!(
            score(converted, converted, opposite)? > score(ordinary, converted, opposite)? + 8.0
        );
    }
    Ok(())
}

#[test]
fn simplified_models_keep_multiple_cpg_neighbors_and_non_taps_likelihoods() -> Result<()> {
    let reference = reference(b"ACGAAACGA");
    let neighbors = [Site::new(1006, C, [T])?];
    let mut reads = copies(10, Strand::OT, &[(1001, T), (1006, C)]);
    reads.extend(copies(10, Strand::OB, &[(1001, C), (1006, T)]));
    let ordinary = infer_snp(&reference, 1001, &neighbors, &reads, &params(true), 3)?;
    let mut params = params(false);
    let result = infer_snp(&reference, 1001, &neighbors, &reads, &params, 3)?;
    assert_eq!(result.sites.len(), 2);
    assert_eq!(result.linking_fragments, 20);
    assert!((probability_sum(&result) - 1.0).abs() < 1e-12);
    params.non_taps = true;
    let control = infer_snp(&reference, 1001, &neighbors, &reads, &params, 3)?;
    assert_eq!(ordinary.genotype_posteriors, control.genotype_posteriors);
    Ok(())
}

#[test]
fn adding_a_neighbor_does_not_change_the_focal_prior() -> Result<()> {
    let neighbor = Site::new(1003, A, [G])?;
    let reads = copies(8, Strand::OT, &[(1003, G)]);
    let evidence = infer_snp(&reference(b"AAAAA"), 1001, &[neighbor], &reads, &params(true), 3)?;
    for alt in &evidence.alleles {
        assert!((*alt.posterior - 0.001 * 4.0 / 9.0).abs() < 1e-12);
        assert!(alt.log_bayes_factor.abs() < 1e-12);
    }
    assert_eq!(evidence.sites.len(), 2);
    Ok(())
}

#[test]
fn converted_cpg_reads_are_explained_without_a_genetic_t() -> Result<()> {
    let mut reads = copies(19, Strand::OT, &[(1001, T)]);
    reads.extend(copies(1, Strand::OT, &[(1001, C)]));
    reads.extend(copies(20, Strand::OB, &[(1001, C)]));
    let taps = infer_snp(&reference(b"ACGA"), 1001, &[], &reads, &params(false), 3)?;
    let ordinary = infer_snp(&reference(b"ACGA"), 1001, &[], &reads, &params(true), 3)?;
    assert!(*allele(&taps, T).posterior < 0.01);
    assert!(*allele(&ordinary, T).heterozygous > 0.999);
    assert!(allele(&taps, T).log_bayes_factor < allele(&ordinary, T).log_bayes_factor);
    Ok(())
}

#[test]
fn opposite_strand_evidence_identifies_a_true_c_t_snp() -> Result<()> {
    let mut reads = copies(20, Strand::OT, &[(1001, T)]);
    reads.extend(copies(10, Strand::OB, &[(1001, C)]));
    reads.extend(copies(10, Strand::OB, &[(1001, T)]));
    let evidence = infer_snp(&reference(b"ACGA"), 1001, &[], &reads, &params(false), 3)?;
    assert!(*allele(&evidence, T).heterozygous > 0.999);
    Ok(())
}

#[test]
fn a_converted_denovo_c_allele_can_look_like_reference_t() -> Result<()> {
    let mut reads = copies(19, Strand::OT, &[(1001, T)]);
    reads.extend(copies(1, Strand::OT, &[(1001, C)]));
    reads.extend(copies(20, Strand::OB, &[(1001, C)]));
    let evidence = infer_snp(&reference(b"ATGA"), 1001, &[], &reads, &params(false), 3)?;
    assert!(*allele(&evidence, C).homozygous > 0.99);
    Ok(())
}

#[test]
fn bottom_strand_conversion_has_the_symmetric_model() -> Result<()> {
    let mut reads = copies(19, Strand::OB, &[(1002, A)]);
    reads.extend(copies(1, Strand::OB, &[(1002, G)]));
    reads.extend(copies(20, Strand::OT, &[(1002, G)]));
    let evidence = infer_snp(&reference(b"ACGA"), 1002, &[], &reads, &params(false), 3)?;
    assert!(*allele(&evidence, A).posterior < 0.01);
    Ok(())
}

#[test]
fn compound_heterozygotes_retain_both_nonreference_alleles() -> Result<()> {
    let mut reads = copies(12, Strand::OB, &[(1001, C)]);
    reads.extend(copies(12, Strand::OT, &[(1001, G)]));
    let evidence = infer_snp(&reference(b"CAG"), 1001, &[], &reads, &params(false), 3)?;
    assert!(*allele(&evidence, C).heterozygous > 0.999);
    assert!(*allele(&evidence, G).heterozygous > 0.999);
    Ok(())
}

#[test]
fn haplotype_sequence_decides_cpg_creation_and_destruction() -> Result<()> {
    let reference = reference(b"ACGA");
    let intact = Haplotype { alleles: vec![(1002, G)], indels: Vec::new() };
    let disrupted = Haplotype { alleles: vec![(1002, T)], indels: Vec::new() };
    assert_eq!(cpg_anchor(&reference, &intact, 1001), Some(1001));
    assert_eq!(cpg_anchor(&reference, &disrupted, 1001), None);
    let reference = super::tests::reference(b"ATGA");
    let created = Haplotype { alleles: vec![(1001, C)], indels: Vec::new() };
    assert_eq!(cpg_anchor(&reference, &created, 1001), Some(1001));
    Ok(())
}

#[test]
fn overlapping_mates_share_the_conversion_event() -> Result<()> {
    let reference = reference(b"ACGA");
    let haplotype = Haplotype { alleles: vec![(1001, C)], indels: Vec::new() };
    let reads =
        Fragment { reads: vec![read(Strand::OT, &[(1001, T)]), read(Strand::OT, &[(1001, T)])] };
    let likelihood =
        fragment_likelihood(&reads, &haplotype, &reference, &[1001], 0.5, &params(false))?.exp();
    let rate = 0.5 * 0.95 + 0.5 * 0.003;
    let expected = rate * (1.0_f64 - 1e-4) + (1.0 - rate) * (1e-4_f64 / 3.0);
    assert!((likelihood - expected).abs() < 1e-6);
    assert!(likelihood > 0.45, "two converted mates are one methylation event, not two");
    let evidence = infer_snp(&reference, 1001, &[], &[reads], &params(false), 1)?;
    assert_eq!(evidence.focal_fragments, 1);
    assert_eq!(evidence.focal_observations, 1);
    assert_eq!(evidence.focal_read_observations, 2);
    assert_eq!(evidence.focal_overlapping_fragments, 1);
    assert_eq!(evidence.focal_discordant_fragments, 0);
    Ok(())
}

#[test]
fn discordant_overlap_is_one_low_confidence_observation_on_both_taps_strands() -> Result<()> {
    for (pos, original, converted, strand) in [(1001, C, T, Strand::OT), (1002, G, A, Strand::OB)] {
        let fragment = Fragment {
            reads: vec![read(strand, &[(pos, original)]), read(strand, &[(pos, converted)])],
        };
        let evidence = infer_snp(&reference(b"ACGA"), pos, &[], &[fragment], &params(false), 1)?;
        assert_eq!(evidence.focal_observations, 1);
        assert_eq!(evidence.focal_read_observations, 2);
        assert_eq!(evidence.focal_discordant_fragments, 1);
        for allele in &evidence.alleles {
            assert!(allele.log_bayes_factor.abs() < 1e-12);
        }
    }
    Ok(())
}

#[test]
fn missing_mapping_quality_and_q0_do_not_create_evidence() -> Result<()> {
    for missing_mapq in [false, true] {
        let mut reads = copies(10, Strand::OT, &[(1001, G)]);
        for fragment in &mut reads {
            for read in &mut fragment.reads {
                if missing_mapq {
                    read.mapq = None;
                } else {
                    for obs in &mut read.observations {
                        obs.quality = Some(0);
                    }
                }
            }
        }
        let evidence = infer_snp(&reference(b"AAA"), 1001, &[], &reads, &params(true), 1)?;
        for alt in &evidence.alleles {
            assert!(alt.log_bayes_factor.abs() < 1e-12);
        }
    }
    Ok(())
}

#[test]
fn unknown_orientation_is_marginalized_once_per_fragment() -> Result<()> {
    let reference = reference(b"ACGA");
    let haplotype = Haplotype { alleles: vec![(1001, C)], indels: Vec::new() };
    let read = fragment(Strand::Unknown, &[(1001, T), (1002, A)]);
    let likelihood =
        fragment_likelihood(&read, &haplotype, &reference, &[1001, 1002], 0.5, &params(false))?
            .exp();
    assert!(likelihood < 0.001, "one molecule cannot independently choose OT at C and OB at G");
    Ok(())
}

#[test]
fn fragment_linkage_distinguishes_cis_and_trans_pairs() -> Result<()> {
    let mut reads = copies(10, Strand::OT, &[(1001, A), (1003, G)]);
    reads.extend(copies(10, Strand::OB, &[(1001, G), (1003, A)]));
    let evidence = infer_snp(
        &reference(b"AAAAA"),
        1001,
        &[Site::new(1003, A, [G])?],
        &reads,
        &params(true),
        3,
    )?;
    let best = evidence
        .pairs
        .iter()
        .max_by(|a, b| a.posterior.total_cmp(&b.posterior))
        .expect("best pair");
    let haplotype = |index: usize| {
        evidence
            .haplotypes
            .get(index)
            .expect("haplotype")
            .alleles
            .iter()
            .map(|(_, base)| *base)
            .collect::<Vec<_>>()
    };
    let a = haplotype(best.first);
    let b = haplotype(best.second);
    assert!((a == [A, G] && b == [G, A]) || (a == [G, A] && b == [A, G]));
    assert!(*best.posterior > 0.999);
    Ok(())
}

#[test]
fn unrelated_cpg_rates_retain_their_neighbor_links() -> Result<()> {
    let reference = reference(b"ACGAACGA");
    let reads = copies(4, Strand::OT, &[(1001, T), (1005, T)]);
    let evidence =
        infer_snp(&reference, 1001, &[Site::new(1005, C, [T])?], &reads, &params(false), 3)?;
    assert_eq!(evidence.sites.len(), 2);
    assert_eq!(evidence.linking_fragments, 4);
    Ok(())
}

#[test]
fn overlapping_singleton_alternates_cannot_amplify_the_three_smoke_patterns() -> Result<()> {
    for (context, ref_ot, alt_ot, ref_ob, alt, non_taps) in [
        (b"ACG".as_slice(), 12, 2, 10, T, false),
        (b"CAG".as_slice(), 5, 0, 3, G, false),
        (b"CAT".as_slice(), 5, 0, 7, T, true),
    ] {
        let reference = reference(context);
        let reference_base = reference.base_at(1001).ok_or_eyre("reference")?;
        let mut reads = copies(ref_ot, Strand::OT, &[(1001, reference_base)]);
        reads.extend(copies(alt_ot, Strand::OT, &[(1001, alt)]));
        reads.extend(copies(ref_ob, Strand::OB, &[(1001, reference_base)]));
        let no_alt = infer_snp(&reference, 1001, &[], &reads, &params(non_taps), 1)?;
        let alt_read = read(Strand::OB, &[(1001, alt)]);
        reads.push(Fragment { reads: vec![alt_read.clone()] });
        let single = infer_snp(&reference, 1001, &[], &reads, &params(non_taps), 1)?;
        reads.last_mut().ok_or_eyre("alt fragment")?.reads.push(alt_read);
        let paired = infer_snp(&reference, 1001, &[], &reads, &params(non_taps), 1)?;
        assert_eq!(paired.focal_fragments, ref_ot + alt_ot + ref_ob + 1);
        assert_eq!(paired.focal_observations, paired.focal_fragments);
        assert_eq!(paired.focal_read_observations, paired.focal_fragments + 1);
        assert!(*allele(&no_alt, alt).posterior < 0.01);
        assert!(*allele(&paired, alt).posterior < 0.05, "{context:?}: {:?}", allele(&paired, alt));
        assert!(
            (allele(&paired, alt).log_bayes_factor - allele(&single, alt).log_bayes_factor).abs()
                < 1e-12
        );
        for obs in reads
            .last_mut()
            .ok_or_eyre("alt fragment")?
            .reads
            .iter_mut()
            .flat_map(|read| &mut read.observations)
        {
            obs.quality = Some(20);
        }
        let lower_quality = infer_snp(&reference, 1001, &[], &reads, &params(non_taps), 1)?;
        assert!(
            allele(&lower_quality, alt).log_bayes_factor
                < allele(&paired, alt).log_bayes_factor - 4.0
        );
    }
    Ok(())
}

#[test]
fn the_site_cap_keeps_the_nearest_neighbours() -> Result<()> {
    let reference = reference(b"ACGAAACGAACA");
    let reads = copies(6, Strand::OT, &[(1001, T), (1006, C), (1010, G)]);
    let neighbors = [Site::new(1006, C, [T])?, Site::new(1010, C, [G])?];
    let capped = infer_snp(&reference, 1001, &neighbors, &reads, &params(false), 1)?;
    assert_eq!(capped.sites.len(), 1);
    let expanded = infer_snp(&reference, 1001, &neighbors, &reads, &params(false), 3)?;
    assert_eq!(expanded.sites.len(), 3);
    assert_eq!(expanded.linking_fragments, 6);
    assert!((probability_sum(&expanded) - 1.0).abs() < 1e-12);
    Ok(())
}

#[test]
fn malformed_candidates_and_unknown_reference_are_errors() -> Result<()> {
    assert!(infer_snp(&reference(b"ANA"), 1001, &[], &[], &params(false), 3).is_err());
    let site = Site::new(1003, C, [T])?;
    assert!(infer_snp(&reference(b"AAAAA"), 1001, &[site], &[], &params(false), 3).is_err());
    Ok(())
}

#[test]
fn adjacent_variants_form_a_cpg_only_on_the_cis_haplotype() -> Result<()> {
    let mut reads = copies(8, Strand::OT, &[(1001, T), (1002, G)]);
    reads.extend(copies(8, Strand::OB, &[(1001, C), (1002, A)]));
    reads.extend(copies(8, Strand::OT, &[(1001, T), (1002, T)]));
    reads.extend(copies(8, Strand::OB, &[(1001, T), (1002, T)]));
    let evidence = infer_snp(
        &reference(b"ATTA"),
        1001,
        &[Site::new(1002, T, [G])?],
        &reads,
        &params(false),
        3,
    )?;
    let best = evidence
        .pairs
        .iter()
        .max_by(|a, b| a.posterior.total_cmp(&b.posterior))
        .expect("best pair");
    let haplotype = |index: usize| {
        evidence
            .haplotypes
            .get(index)
            .expect("haplotype")
            .alleles
            .iter()
            .map(|(_, base)| *base)
            .collect::<Vec<_>>()
    };
    let a = haplotype(best.first);
    let b = haplotype(best.second);
    assert!((a == [C, G] && b == [T, T]) || (a == [T, T] && b == [C, G]));
    assert!(*best.posterior > 0.999);
    Ok(())
}

#[test]
fn fully_converted_strand_alone_cannot_distinguish_c_from_t() -> Result<()> {
    let reference = reference(b"ACGA");
    let c = Haplotype { alleles: vec![(1001, C)], indels: Vec::new() };
    let t = Haplotype { alleles: vec![(1001, T)], indels: Vec::new() };
    let fragment = fragment(Strand::OT, &[(1001, T)]);
    let mut params = params(false);
    params.conversion_rate = Probability::ONE;
    let c_likelihood = fragment_likelihood(&fragment, &c, &reference, &[1001], 1.0, &params)?;
    let t_likelihood = fragment_likelihood(&fragment, &t, &reference, &[1001], 1.0, &params)?;
    assert!((c_likelihood - t_likelihood).abs() < 1e-12);
    Ok(())
}

#[test]
fn unknown_observed_bases_do_not_change_mapping_evidence() -> Result<()> {
    let reference = reference(b"AAA");
    let haplotype = Haplotype { alleles: vec![(1001, A)], indels: Vec::new() };
    let plain = fragment(Strand::OT, &[(1001, A)]);
    let unknown = fragment(Strand::OT, &[(1001, A), (1002, Unknown)]);
    let a = fragment_likelihood(&plain, &haplotype, &reference, &[1001, 1002], 0.0, &params(true))?;
    let b =
        fragment_likelihood(&unknown, &haplotype, &reference, &[1001, 1002], 0.0, &params(true))?;
    assert!((a - b).abs() < 1e-12);
    Ok(())
}

#[test]
fn a_known_mate_resolves_an_unknown_molecular_orientation() -> Result<()> {
    let reference = reference(b"ACGA");
    let haplotype = Haplotype { alleles: vec![(1001, C)], indels: Vec::new() };
    let known =
        Fragment { reads: vec![read(Strand::OB, &[(1001, T)]), read(Strand::OB, &[(1001, T)])] };
    let unknown = Fragment {
        reads: vec![read(Strand::Unknown, &[(1001, T)]), read(Strand::OB, &[(1001, T)])],
    };
    let a = fragment_likelihood(&known, &haplotype, &reference, &[1001], 0.5, &params(false))?;
    let b = fragment_likelihood(&unknown, &haplotype, &reference, &[1001], 0.5, &params(false))?;
    assert!((a - b).abs() < 1e-12);
    Ok(())
}

#[test]
fn disagreeing_mate_assignments_are_one_uncertain_orientation() -> Result<()> {
    let reference = reference(b"ACGA");
    let haplotype = Haplotype { alleles: vec![(1001, C)], indels: Vec::new() };
    let disagreeing =
        Fragment { reads: vec![read(Strand::OT, &[(1001, T)]), read(Strand::OB, &[(1001, T)])] };
    let unknown = Fragment {
        reads: vec![read(Strand::Unknown, &[(1001, T)]), read(Strand::Unknown, &[(1001, T)])],
    };
    let a =
        fragment_likelihood(&disagreeing, &haplotype, &reference, &[1001], 0.5, &params(false))?;
    let b = fragment_likelihood(&unknown, &haplotype, &reference, &[1001], 0.5, &params(false))?;
    assert!((a - b).abs() < 1e-12);
    Ok(())
}
