use super::*;
use Base::*;

fn params(non_taps: bool) -> LikelihoodParams {
    LikelihoodParams {
        conversion_rate: Probability::new(0.95).expect("probability"),
        false_conversion_rate: Probability::new(0.003).expect("probability"),
        snp_prior: Probability::new(0.001).expect("probability"),
        non_taps,
        quadrature_order: 16,
    }
}

fn reference(bases: &[u8]) -> Reference {
    Reference { start: 1000, bases: bases.iter().copied().map(Base::from).collect() }
}

fn read(strand: Strand, bases: &[(u32, Base)]) -> ReadEvidence {
    ReadEvidence {
        strand,
        mapq: Some(60),
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
    let intact = Haplotype { alleles: vec![(1002, G)] };
    let disrupted = Haplotype { alleles: vec![(1002, T)] };
    assert_eq!(cpg_anchor(&reference, &intact, 1001), Some(1001));
    assert_eq!(cpg_anchor(&reference, &disrupted, 1001), None);
    let reference = super::tests::reference(b"ATGA");
    let created = Haplotype { alleles: vec![(1001, C)] };
    assert_eq!(cpg_anchor(&reference, &created, 1001), Some(1001));
    Ok(())
}

#[test]
fn overlapping_mates_share_the_conversion_event() -> Result<()> {
    let reference = reference(b"ACGA");
    let haplotype = Haplotype { alleles: vec![(1001, C)] };
    let reads =
        Fragment { reads: vec![read(Strand::OT, &[(1001, T)]), read(Strand::OT, &[(1001, T)])] };
    let likelihood =
        fragment_likelihood(&reads, &haplotype, &reference, &[1001], 0.5, &params(false))?.exp();
    let rate = 0.5 * 0.95 + 0.5 * 0.003;
    let expected = rate * (1.0_f64 - 1e-4).powi(2) + (1.0 - rate) * (1e-4_f64 / 3.0).powi(2);
    assert!((likelihood - expected).abs() < 1e-6);
    assert!(likelihood > 0.45, "two converted mates are one methylation event, not two");
    let evidence = infer_snp(&reference, 1001, &[], &[reads], &params(false), 1)?;
    assert_eq!(evidence.focal_fragments, 1);
    assert_eq!(evidence.focal_observations, 2);
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
    let haplotype = Haplotype { alleles: vec![(1001, C)] };
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
    Ok(())
}

#[test]
fn adjacent_background_cpg_does_not_break_an_isolated_snp() -> Result<()> {
    let evidence = infer_snp(&reference(b"CGCG"), 1001, &[], &[], &params(false), 3)?;
    assert!(!evidence.cpg_context_limited);
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
    let c = Haplotype { alleles: vec![(1001, C)] };
    let t = Haplotype { alleles: vec![(1001, T)] };
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
    let haplotype = Haplotype { alleles: vec![(1001, A)] };
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
    let haplotype = Haplotype { alleles: vec![(1001, C)] };
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
    let haplotype = Haplotype { alleles: vec![(1001, C)] };
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
