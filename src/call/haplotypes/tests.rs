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
fn bounded_expansion_follows_transitive_fragment_links_beyond_three_sites() -> Result<()> {
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
    let exact = infer_snp_expanded(
        &reference,
        1001,
        &neighbors,
        &[],
        &reads,
        &params(true),
        &ExpansionParams::default(),
    )?;
    let bounded = infer_snp_expanded(
        &reference,
        1001,
        &neighbors,
        &[],
        &reads,
        &params(true),
        &ExpansionParams { max_haplotypes: 16, ..ExpansionParams::default() },
    )?;
    for result in [&exact, &bounded] {
        assert_eq!(result.sites.len(), 5);
        assert_eq!(result.expansion_steps, 4);
        assert!(*allele(result, G).heterozygous > 0.999);
        assert!((probability_sum(result) - 1.0).abs() < 1e-12);
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
    }
    assert!(bounded.haplotypes.len() <= 16);
    assert!(bounded.pruned_haplotypes > 0);
    assert!(bounded.maximum_pruned_mass < 1e-6);
    assert!((*allele(&exact, G).posterior - *allele(&bounded, G).posterior).abs() < 1e-9);
    Ok(())
}

#[test]
fn disconnected_and_low_confidence_fragments_do_not_extend_a_region() -> Result<()> {
    let reference = reference(b"AAAAAAAAA");
    let neighbors = [Site::new(1004, A, [G])?];
    let read = read(Strand::OT, &[(1001, A), (1004, G)]);
    let duplicated = vec![Fragment { reads: vec![read.clone(), read] }];
    let result = infer_snp_expanded(
        &reference,
        1001,
        &neighbors,
        &[],
        &duplicated,
        &params(true),
        &ExpansionParams::default(),
    )?;
    assert_eq!(result.sites.len(), 1);
    assert_eq!(result.disconnected_sites, 1, "two mates are one linking molecule");
    let mut conflict = duplicated.clone();
    conflict.extend(duplicated);
    for fragment in &mut conflict {
        let other = fragment.reads.get_mut(1).ok_or_eyre("mate")?;
        let observation = other.observations.get_mut(1).ok_or_eyre("observation")?;
        observation.base = A;
    }
    let result = infer_snp_expanded(
        &reference,
        1001,
        &neighbors,
        &[],
        &conflict,
        &params(true),
        &ExpansionParams::default(),
    )?;
    assert_eq!(result.sites.len(), 1, "Q0 consensus cannot provide a connection");
    let mut independent = copies(10, Strand::OT, &[(1001, A)]);
    independent.extend(copies(10, Strand::OT, &[(1004, G)]));
    let result = infer_snp_expanded(
        &reference,
        1001,
        &neighbors,
        &[],
        &independent,
        &params(true),
        &ExpansionParams::default(),
    )?;
    assert_eq!(result.disconnected_sites, 1);
    Ok(())
}

#[test]
fn expansion_limits_preserve_focal_genotypes_and_report_omissions() -> Result<()> {
    let reference = reference(b"AAAAAAAAAAAAAAAA");
    let neighbors =
        (1002..1010).map(|pos| Site::new(pos, A, [C, G, T])).collect::<Result<Vec<_>>>()?;
    let bases = (1001..1010).map(|pos| (pos, A)).collect::<Vec<_>>();
    let mut reads = copies(8, Strand::OT, &bases);
    let alternate = (1001..1010).map(|pos| (pos, G)).collect::<Vec<_>>();
    reads.extend(copies(8, Strand::OB, &alternate));
    let search = ExpansionParams { max_sites: 4, max_haplotypes: 32, ..ExpansionParams::default() };
    let result =
        infer_snp_expanded(&reference, 1001, &neighbors, &[], &reads, &params(true), &search)?;
    assert_eq!(result.sites.len(), 4);
    assert_eq!(result.resource_limited_sites, 5);
    assert_eq!(result.omitted_sites, 5);
    assert!(result.haplotypes.len() <= 32);
    let genotypes = result
        .pairs
        .iter()
        .map(|pair| {
            let (a, b) = pair_bases(pair, &result.haplotypes, 1001)?;
            genotype_index(a, b).ok_or_eyre("genotype")
        })
        .collect::<Result<BTreeSet<_>>>()?;
    assert_eq!(genotypes.len(), 10);
    assert!(result.pairs.iter().all(|pair| pair.log_prior.is_finite()));
    Ok(())
}

#[test]
fn retained_phase_assignments_share_one_genotype_prior_after_pruning() -> Result<()> {
    let sites = [Site::new(1001, A, [G])?, Site::new(1002, A, [G])?];
    let haplotypes = [
        Haplotype { alleles: vec![(1001, A), (1002, A)], indel: None },
        Haplotype { alleles: vec![(1001, G), (1002, G)], indel: None },
    ];
    let priors = likelihood::retained_pair_priors(&haplotypes, &sites, &params(true), false)?;
    assert!(
        (priors.get(1).ok_or_eyre("het prior")?.exp() - (0.001_f64 / 2.0).powi(2)).abs() < 1e-18
    );
    Ok(())
}

#[test]
fn buffered_branches_can_return_when_downstream_evidence_changes_the_ranking() -> Result<()> {
    let reference = reference(b"AAAAAAAAAAAAAAAAAAAA");
    let positions = [1001, 1004, 1007, 1010, 1013, 1016];
    let neighbors =
        positions.iter().skip(1).map(|&pos| Site::new(pos, A, [G])).collect::<Result<Vec<_>>>()?;
    for seed in 1..=32u64 {
        let mut state = seed;
        let mut reads = Vec::new();
        for _ in 0..32 {
            let mut observations = Vec::new();
            for &pos in &positions {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                if state >> 61 != 0 {
                    observations.push((pos, if state & (1 << 40) == 0 { A } else { G }));
                }
            }
            reads.push(fragment(Strand::OT, &observations));
        }
        let search = ExpansionParams { max_haplotypes: 16, ..ExpansionParams::default() };
        let result =
            infer_snp_expanded(&reference, 1001, &neighbors, &[], &reads, &params(true), &search)?;
        if result.rescued_haplotypes > 0 {
            let no_holdouts = infer_snp_expanded(
                &reference,
                1001,
                &neighbors,
                &[],
                &reads,
                &params(true),
                &ExpansionParams { haplotype_holdout: 0, ..search },
            )?;
            assert_eq!(no_holdouts.rescued_haplotypes, 0);
            assert_eq!(result.sites.len(), positions.len());
            assert!(result.haplotypes.len() <= 16);
            return Ok(());
        }
    }
    color_eyre::eyre::bail!("No buffered branch was reconsidered across the deterministic fixtures")
}

#[test]
fn conversion_only_cpg_candidates_do_not_bridge_genetic_regions() -> Result<()> {
    let reference = reference(b"ACGAACGA");
    let neighbors = [Site::new(1005, C, [T])?];
    let mut reads = copies(8, Strand::OT, &[(1001, T), (1005, T)]);
    reads.extend(copies(8, Strand::OB, &[(1001, C), (1005, C)]));
    let result = infer_snp_expanded(
        &reference,
        1001,
        &neighbors,
        &[],
        &reads,
        &params(false),
        &ExpansionParams::default(),
    )?;
    assert_eq!(result.sites.len(), 1);
    reads.extend(copies(4, Strand::OB, &[(1001, C), (1005, T)]));
    let result = infer_snp_expanded(
        &reference,
        1001,
        &neighbors,
        &[],
        &reads,
        &params(false),
        &ExpansionParams::default(),
    )?;
    assert_eq!(result.sites.len(), 2);
    assert_eq!(result.disconnected_sites, 0);
    Ok(())
}

#[test]
fn discovery_reaches_both_150bp_edges_and_respects_the_expansion_bound() -> Result<()> {
    let reference = Reference { start: 800, bases: vec![A; 500] };
    let neighbors = [Site::new(850, A, [G])?, Site::new(1150, A, [G])?, Site::new(1151, A, [G])?];
    let observations = [(850, G), (1000, G), (1150, G), (1151, G)];
    let reads = copies(4, Strand::OT, &observations);
    let result = infer_snp_expanded(
        &reference,
        1000,
        &neighbors,
        &[],
        &reads,
        &params(true),
        &ExpansionParams { expansion_radius: 150, ..ExpansionParams::default() },
    )?;
    assert_eq!(result.sites.iter().map(|site| site.pos).collect::<Vec<_>>(), [850, 1000, 1150]);
    Ok(())
}

#[test]
fn unconverted_denovo_cpg_alleles_link_on_the_conversion_prone_strand() -> Result<()> {
    for (sequence, unconverted, strand) in
        [(b"ATGAAACA".as_slice(), C, Strand::OT), (b"CAAAAACA".as_slice(), G, Strand::OB)]
    {
        let reference = reference(sequence);
        let site =
            Site::new(1001, reference.base_at(1001).ok_or_eyre("reference")?, [unconverted])?;
        let reads = copies(4, strand, &[(1001, unconverted), (1006, C)]);
        let result = infer_snp_expanded(
            &reference,
            1006,
            &[site],
            &[],
            &reads,
            &params(false),
            &ExpansionParams::default(),
        )?;
        assert_eq!(result.sites.len(), 2);
        assert_eq!(result.disconnected_sites, 0);
    }
    Ok(())
}

#[test]
fn distant_nuisance_indels_require_fragment_connections_to_the_focal_region() -> Result<()> {
    let reference = Reference { start: 1000, bases: vec![A; 1000] };
    let indel = Indel { start: 1501, deleted: Vec::new(), inserted: vec![C] };
    let sequence_read = |pos| {
        let mut read = read(Strand::OT, &[(pos, A)]);
        read.sequence = Some(SequenceEvidence {
            bases: vec![SequenceObservation {
                position: SequencePosition::reference(pos),
                base: A,
                quality: Some(40),
            }],
            reference_skip: false,
            indels: Vec::new(),
        });
        read
    };
    let mut reads =
        (0..4).map(|_| Fragment { reads: vec![sequence_read(1001)] }).collect::<Vec<_>>();
    let mut params = params(true);
    params.likelihood_engine = LikelihoodEngine::PairHmm;
    let isolated = infer_snp_expanded(
        &reference,
        1001,
        &[],
        std::slice::from_ref(&indel),
        &reads,
        &params,
        &ExpansionParams::default(),
    )?;
    assert!(isolated.haplotypes.iter().all(|haplotype| haplotype.indel.is_none()));
    assert_eq!(isolated.omitted_indels, 1);
    for fragment in &mut reads {
        fragment.reads.push(sequence_read(1501));
    }
    let connected = infer_snp_expanded(
        &reference,
        1001,
        &[],
        &[indel],
        &reads,
        &params,
        &ExpansionParams::default(),
    )?;
    assert!(connected.haplotypes.iter().any(|haplotype| haplotype.indel.is_some()));
    Ok(())
}

#[test]
fn cpg_partner_hypotheses_are_joint_when_a_variant_removes_conversion_context() -> Result<()> {
    let reference = reference(b"ACGAAAAA");
    let neighbors = [Site::new(1001, C, [T])?];
    let mut reads = copies(6, Strand::OT, &[(1001, T), (1002, G)]);
    reads.extend(copies(2, Strand::OT, &[(1001, C), (1002, G)]));
    reads.extend(copies(4, Strand::OT, &[(1001, T), (1002, A)]));
    let result = infer_snp_expanded(
        &reference,
        1002,
        &neighbors,
        &[],
        &reads,
        &params(false),
        &ExpansionParams::default(),
    )?;
    assert_eq!(
        result.sites.len(),
        2,
        "conversion-ambiguous partner still changes candidate sequence context"
    );
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
    let result = infer_snp_expanded(
        &reference,
        1001,
        &neighbors,
        &[],
        &reads,
        &params(false),
        &ExpansionParams::default(),
    )?;
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
    let mut params = params(false);
    let result = infer_snp(&reference, 1001, &neighbors, &reads, &params, 3)?;
    assert_eq!(result.sites.len(), 2);
    assert_eq!(result.retained_linking_fragments, 20);
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
    assert_eq!(evidence.omitted_sites, 0);
    assert_eq!(evidence.offered_linking_fragments, 4);
    assert_eq!(evidence.retained_linking_fragments, 4);
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
fn site_cap_reports_omissions_without_cpg_link_loss() -> Result<()> {
    let reference = reference(b"ACGAAACGAACA");
    let reads = copies(6, Strand::OT, &[(1001, T), (1006, C), (1010, G)]);
    let neighbors = [Site::new(1006, C, [T])?, Site::new(1010, C, [G])?];
    let capped = infer_snp(&reference, 1001, &neighbors, &reads, &params(false), 1)?;
    assert_eq!(capped.omitted_sites, 2);
    assert_eq!(capped.resource_limited_sites, 2);
    let expanded = infer_snp(&reference, 1001, &neighbors, &reads, &params(false), 3)?;
    assert_eq!(expanded.sites.len(), 3);
    assert_eq!(expanded.retained_linking_fragments, 6);
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
