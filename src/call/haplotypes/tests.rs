use super::*;
use Base::*;

pub(super) fn params(non_taps: bool) -> LikelihoodParams {
    LikelihoodParams {
        conversion_rate: Probability::new(0.95).expect("probability"),
        false_conversion_rate: Probability::new(0.003).expect("probability"),
        snp_prior: Probability::new(0.001).expect("probability"),
        non_taps,
        quadrature_order: 16,
        likelihood_engine: LikelihoodEngine::Aligned,
        methylation_model: MethylationModel::Integrated,
        methylation_prior_strength: 4.0,
        methylation_min_cpg_probability: Probability::new(0.9).expect("probability"),
        gap_open: Probability::new(10.0_f64.powf(-4.5)).expect("probability"),
        gap_extend: Probability::new(0.1).expect("probability"),
        indel_prior: Probability::new(0.0001).expect("probability"),
        realignment_flank: 12,
        max_indel_length: 20,
    }
}

fn reference(bases: &[u8]) -> Reference {
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
        let mut params = params(false);
        params.methylation_model = MethylationModel::Neutral;
        let context = methylation::Context::new(&reference, &haplotypes, &loci, None)?;
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
    for model in [MethylationModel::Neutral, MethylationModel::Shrinkage] {
        let mut params = params(false);
        params.methylation_model = model;
        let result = infer_snp(&reference, 1001, &neighbors, &reads, &params, 3)?;
        assert_eq!(result.sites.len(), 2);
        assert_eq!(result.cpg_omitted_sites, 0);
        assert_eq!(result.retained_linking_fragments, 20);
        assert!((probability_sum(&result) - 1.0).abs() < 1e-12);
        params.non_taps = true;
        let control = infer_snp(&reference, 1001, &neighbors, &reads, &params, 3)?;
        assert_eq!(ordinary.genotype_posteriors, control.genotype_posteriors);
    }
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
    let intact = Haplotype { alleles: vec![(1002, G)], indel: None };
    let disrupted = Haplotype { alleles: vec![(1002, T)], indel: None };
    assert_eq!(cpg_anchor(&reference, &intact, 1001), Some(1001));
    assert_eq!(cpg_anchor(&reference, &disrupted, 1001), None);
    let reference = super::tests::reference(b"ATGA");
    let created = Haplotype { alleles: vec![(1001, C)], indel: None };
    assert_eq!(cpg_anchor(&reference, &created, 1001), Some(1001));
    Ok(())
}

#[test]
fn overlapping_mates_share_the_conversion_event() -> Result<()> {
    let reference = reference(b"ACGA");
    let haplotype = Haplotype { alleles: vec![(1001, C)], indel: None };
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
    let haplotype = Haplotype { alleles: vec![(1001, C)], indel: None };
    let read = fragment(Strand::Unknown, &[(1001, T), (1002, A)]);
    let likelihood =
        fragment_likelihood(&read, &haplotype, &reference, &[1001, 1002], 0.5, &params(false))?
            .exp();
    assert!(likelihood < 0.001, "one molecule cannot independently choose OT at C and OB at G");
    Ok(())
}

#[test]
fn methylation_integral_matches_an_analytic_small_problem() -> Result<()> {
    let reads = copies(10, Strand::OT, &[(1001, T)]);
    let evidence = infer_snp(&reference(b"ACGA"), 1001, &[], &reads, &params(false), 1)?;
    let pair = evidence
        .pairs
        .iter()
        .find(|pair| {
            pair_bases(pair, &evidence.haplotypes, 1001).is_ok_and(|bases| bases == (C, C))
        })
        .expect("C/C hypothesis");
    let emission =
        |rate: f64| (1.0 - 1e-6) * (rate * (1.0 - 1e-4) + (1.0 - rate) * 1e-4 / 3.0) + 1e-6 * 0.25;
    let (a, b) = (emission(0.003), emission(0.95));
    let expected = (b.powi(11) - a.powi(11)) / (11.0 * (b - a));
    assert!((pair.log_likelihood - expected.ln()).abs() < 1e-10);
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
fn unrelated_cpg_rates_are_not_silently_tied_together() -> Result<()> {
    let reference = reference(b"ACGAACGA");
    let reads = copies(4, Strand::OT, &[(1001, T), (1005, T)]);
    let evidence =
        infer_snp(&reference, 1001, &[Site::new(1005, C, [T])?], &reads, &params(false), 3)?;
    assert_eq!(evidence.sites.len(), 1);
    assert_eq!(evidence.omitted_sites, 1);
    assert!(evidence.cpg_context_limited);
    assert_eq!(evidence.cpg_omitted_sites, 1);
    assert_eq!(evidence.offered_linking_fragments, 4);
    assert_eq!(evidence.retained_linking_fragments, 0);
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
fn cpg_reduction_reports_lost_links_separately_from_the_site_cap() -> Result<()> {
    let reference = reference(b"ACGAAACGAACA");
    let reads = copies(6, Strand::OT, &[(1001, T), (1006, C), (1010, G)]);
    let neighbors = [Site::new(1006, C, [T])?, Site::new(1010, C, [G])?];
    let capped = infer_snp(&reference, 1001, &neighbors, &reads, &params(false), 1)?;
    assert_eq!(capped.omitted_sites, 2);
    assert_eq!(capped.cpg_omitted_sites, 0);
    assert_eq!(capped.offered_linking_fragments, 0);
    let reduced = infer_snp(&reference, 1001, &neighbors, &reads, &params(false), 3)?;
    assert_eq!(reduced.cpg_omitted_sites, 2);
    assert_eq!(reduced.offered_linking_fragments, 6);
    assert_eq!(reduced.retained_linking_fragments, 0);
    assert!((reduced.genotype_posteriors.iter().map(|p| **p).sum::<f64>() - 1.0).abs() < 1e-12);
    Ok(())
}

#[test]
fn adjacent_background_cpg_does_not_break_an_isolated_snp() -> Result<()> {
    let evidence = infer_snp(&reference(b"CGCG"), 1001, &[], &[], &params(false), 3)?;
    assert!(!evidence.cpg_context_limited);
    Ok(())
}

#[test]
fn independent_two_cpg_oracle_recovers_linkage_against_a_singleton_error() -> Result<()> {
    use super::likelihood::{log_add, log_sum, pair_prior, quadrature, sequencing_emission};
    let reference = reference(b"CAGAAACGAA");
    let mut reads = copies(2, Strand::OT, &[(1001, A), (1007, G)]);
    reads.extend(copies(3, Strand::OT, &[(1001, A), (1007, T)]));
    reads.extend(copies(3, Strand::OB, &[(1001, A), (1007, G)]));
    let alternate = read(Strand::OB, &[(1001, G), (1007, T)]);
    reads.push(Fragment { reads: vec![alternate.clone(), alternate] });
    for read in reads.iter_mut().flat_map(|fragment| &mut fragment.reads) {
        read.mapq = Some(254);
    }
    let params = params(false);
    let neighbor = Site::new(1007, G, [T])?;
    let bounded = infer_snp(&reference, 1001, std::slice::from_ref(&neighbor), &reads, &params, 3)?;
    assert_eq!(bounded.cpg_omitted_sites, 1);
    assert_eq!(bounded.offered_linking_fragments, 9);
    assert_eq!(bounded.retained_linking_fragments, 0);
    let sites = [Site::new(1001, A, Base::KNOWN)?, neighbor];
    let haplotypes = enumerate(&sites);
    let loci = scored_loci(&reference, &sites);
    // Nine molecules give a polynomial of degree at most nine per CpG rate.
    // Five Gauss-Legendre nodes integrate it exactly under the uniform prior.
    let nodes = quadrature(5)?;
    let mut matrix = Vec::new();
    for haplotype in &haplotypes {
        let anchors: Vec<_> = loci
            .iter()
            .filter_map(|&pos| cpg_anchor(&reference, haplotype, pos))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        assert!(anchors.len() <= 2);
        let mut grid = vec![(Vec::new(), 0.0)];
        for _ in &anchors {
            grid = grid
                .into_iter()
                .flat_map(|(rates, weight)| {
                    nodes.iter().map(move |node| {
                        let mut rates = rates.clone();
                        rates.push(node.methylation);
                        (rates, weight + node.log_weight)
                    })
                })
                .collect();
        }
        let mut values = Vec::new();
        for (rates, weight) in grid {
            let mut likelihoods = Vec::new();
            for fragment in &reads {
                let strand = fragment.reads.first().ok_or_eyre("read")?.strand;
                let mut likelihood = 0.0;
                for &pos in &loci {
                    let base = haplotype.base_at(&reference, pos).ok_or_eyre("base")?;
                    let methylation = cpg_anchor(&reference, haplotype, pos)
                        .and_then(|anchor| anchors.iter().position(|&a| a == anchor))
                        .and_then(|i| rates.get(i))
                        .copied();
                    let rate = params.conversion(methylation);
                    let modified = [
                        crate::metrics::methylation::CpgSide::C,
                        crate::metrics::methylation::CpgSide::G,
                    ]
                    .into_iter()
                    .find(|side| side.unmod_base() == base && side.strand() == strand)
                    .map_or(base, |side| side.mod_base());
                    let observed = super::likelihood::consensus_at(fragment, pos)?;
                    let unconverted =
                        sequencing_emission(observed.base, base, observed.quality).ln();
                    let converted =
                        sequencing_emission(observed.base, modified, observed.quality).ln();
                    likelihood += log_add((1.0 - rate).ln() + unconverted, rate.ln() + converted);
                }
                likelihoods.push(likelihood);
            }
            values.push((weight, likelihoods));
        }
        matrix.push(values);
    }
    let mut total = Vec::new();
    let mut present = Vec::new();
    for (first, a) in haplotypes.iter().enumerate() {
        for (second, b) in haplotypes.iter().enumerate().skip(first) {
            let a_values = matrix.get(first).ok_or_eyre("first haplotype")?;
            let b_values = matrix.get(second).ok_or_eyre("second haplotype")?;
            let integrated = if first == second {
                log_sum(a_values.iter().map(|(weight, values)| weight + values.iter().sum::<f64>()))
            } else {
                log_sum(a_values.iter().flat_map(|(weight_a, values_a)| {
                    b_values.iter().map(move |(weight_b, values_b)| {
                        weight_a
                            + weight_b
                            + values_a
                                .iter()
                                .zip(values_b)
                                .map(|(&a, &b)| log_add(a, b) - std::f64::consts::LN_2)
                                .sum::<f64>()
                    })
                }))
            };
            let joint = integrated
                + pair_prior(a, b, &sites, *params.snp_prior, *params.indel_prior, false)?;
            total.push(joint);
            if a.allele_at(1001) == Some(G) || b.allele_at(1001) == Some(G) {
                present.push(joint);
            }
        }
    }
    let oracle = (log_sum(present) - log_sum(total)).exp();
    assert!(*allele(&bounded, G).posterior > oracle * 100.0);
    assert!(oracle < 1e-6, "independent CpG posterior={oracle}, bounded={:?}", allele(&bounded, G));
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
    let c = Haplotype { alleles: vec![(1001, C)], indel: None };
    let t = Haplotype { alleles: vec![(1001, T)], indel: None };
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
    let haplotype = Haplotype { alleles: vec![(1001, A)], indel: None };
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
    let haplotype = Haplotype { alleles: vec![(1001, C)], indel: None };
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
    let haplotype = Haplotype { alleles: vec![(1001, C)], indel: None };
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
