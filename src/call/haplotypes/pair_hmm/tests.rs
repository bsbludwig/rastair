use super::super::likelihood::sequencing_emission;
use super::*;
use crate::call::haplotypes::{
    LikelihoodEngine, Observation, ReadEvidence, infer_snp, infer_snp_with_indels,
};
use seqair_types::Probability;

fn params(non_taps: bool) -> LikelihoodParams {
    let mut params = crate::call::haplotypes::tests::params(non_taps);
    params.likelihood_engine = LikelihoodEngine::PairHmm;
    params
}

fn reference(sequence: &[u8]) -> Reference {
    Reference { start: 1000, bases: sequence.iter().copied().map(Base::from).collect() }
}

fn read(sequence: &[u8], positions: &[SequencePosition], strand: Strand) -> ReadEvidence {
    assert_eq!(sequence.len(), positions.len());
    let bases: Vec<_> = sequence
        .iter()
        .zip(positions)
        .map(|(&base, &position)| SequenceObservation {
            position,
            base: Base::from(base),
            quality: Some(40),
        })
        .collect();
    ReadEvidence {
        annotations: None,
        mapq: Some(60),
        strand,
        observations: bases
            .iter()
            .filter(|base| base.position.is_reference())
            .map(|base| Observation {
                pos: base.position.pos,
                base: base.base,
                quality: base.quality,
            })
            .collect(),
        sequence: Some(SequenceEvidence { bases, reference_skip: false, indels: Vec::new() }),
    }
}

fn aligned_read(sequence: &[u8], strand: Strand) -> ReadEvidence {
    read(
        sequence,
        &(0..sequence.len())
            .map(|offset| {
                SequencePosition::reference(1000 + u32::try_from(offset).expect("short fixture"))
            })
            .collect::<Vec<_>>(),
        strand,
    )
}

#[test]
fn forward_agrees_with_exhaustive_alignment_paths() -> Result<()> {
    fn enumerate(
        i: usize,
        j: usize,
        state: u8,
        rows: usize,
        columns: usize,
        params: &LikelihoodParams,
    ) -> f64 {
        if i == rows {
            return if state == 2 { 0.0 } else { 1.0 };
        }
        let mut result = 0.0;
        if j < columns {
            let transition =
                if state == 0 { 1.0 - 2.0 * *params.gap_open } else { 1.0 - *params.gap_extend };
            let emission = 0.13 + ((i * 3 + j * 2) % 7) as f64 * 0.09;
            result += transition * emission * enumerate(i + 1, j + 1, 0, rows, columns, params);
        }
        if state != 2 && j > 0 {
            let transition = if state == 0 { *params.gap_open } else { *params.gap_extend };
            result += transition * 0.25 * enumerate(i + 1, j, 1, rows, columns, params);
        }
        if state != 1 && j < columns && i > 0 {
            let transition = if state == 0 { *params.gap_open } else { *params.gap_extend };
            result += transition * enumerate(i, j + 1, 2, rows, columns, params);
        }
        result
    }
    let mut params = params(true);
    params.gap_open = Probability::new(0.07)?;
    params.gap_extend = Probability::new(0.31)?;
    for rows in 1..=4 {
        for columns in 1..=5 {
            let expected: f64 = (0..columns)
                .map(|j| enumerate(0, j, 2, rows, columns, &params) / columns as f64)
                .sum();
            let actual = forward(
                rows,
                columns,
                &params,
                |_| 0.25_f64.ln(),
                |i, j| (0.13 + ((i * 3 + j * 2) % 7) as f64 * 0.09).ln(),
            )?
            .exp();
            assert!((actual - expected).abs() < 1e-13, "{rows}/{columns}: {actual} != {expected}");
        }
    }
    Ok(())
}

#[test]
fn ordinary_emissions_form_a_normalized_observation_model() -> Result<()> {
    let params = params(true);
    for length in 1..=3u32 {
        let mass = forward(length as usize, 4, &params, |_| 0.0, |_, _| 0.0)?;
        let mut probability = 0.0;
        for mut code in 0..4usize.pow(length) {
            let mut observed = Vec::new();
            for _ in 0..length {
                observed.push(*Base::KNOWN.get(code % 4).expect("base"));
                code /= 4;
            }
            probability += (forward(
                length as usize,
                4,
                &params,
                |_| 0.25_f64.ln(),
                |i, j| {
                    sequencing_emission(
                        *observed.get(i).expect("read"),
                        *Base::KNOWN.get(j).expect("haplotype"),
                        Some(30),
                    )
                    .ln()
                },
            )? - mass)
                .exp();
        }
        assert!((probability - 1.0).abs() < 1e-12);
    }
    Ok(())
}

#[test]
fn unmapped_profile_shortcut_matches_the_forward_sum_with_unknown_bases() -> Result<()> {
    let params = params(false);
    let reference = reference(b"ACGTTAGCATACGT");
    let fragment = Fragment { reads: vec![aligned_read(b"ACNTTAG", Strand::OT)] };
    let sites = [Site::new(1002, Base::G, Base::KNOWN)?];
    for indel in
        [None, Some(Indel { start: 1004, deleted: Vec::new(), inserted: vec![Base::C, Base::G] })]
    {
        let haplotype = Haplotype { alleles: vec![(1002, Base::G)], indel: indel.clone() };
        let prepared = Prepared::new(&reference, &sites, indel.as_ref(), &[&fragment], &params)?;
        let template = prepared.template(&reference, &haplotype, &[1002], &params)?;
        let profile = prepared
            .fragments
            .first()
            .and_then(Option::as_ref)
            .and_then(|molecule| molecule.joint.first())
            .ok_or_eyre("profile")?;
        let uniform = |index: usize| {
            profile.get(index).map_or(f64::NEG_INFINITY, |base| {
                if base.base.known_index().is_some() { 0.25_f64.ln() } else { 0.0 }
            })
        };
        let expected =
            forward(profile.len(), template.bases.len(), &params, uniform, |i, _| uniform(i))?
                - forward(profile.len(), template.bases.len(), &params, |_| 0.0, |_, _| 0.0)?;
        let actual = score_with_mapping(profile, &template, Strand::OT, 0.9, &params, 1)?;
        assert!((actual - expected).abs() < 1e-12);
    }
    Ok(())
}

#[test]
fn long_read_likelihoods_do_not_underflow() -> Result<()> {
    let params = params(true);
    let length = 550;
    let value = forward(length, length, &params, |_| 0.25_f64.ln(), |_, _| 0.25_f64.ln())?;
    let mass = forward(length, length, &params, |_| 0.0, |_, _| 0.0)?;
    assert!((value - mass + length as f64 * 4.0_f64.ln()).abs() < 1e-9);
    Ok(())
}

#[test]
fn overlapping_mates_share_one_conversion_draw_in_the_match_emission() -> Result<()> {
    let params = params(false);
    let template = Template {
        bases: vec![TemplateBase {
            position: SequencePosition::reference(1000),
            base: Base::C,
            cpg: true,
            cpg_key: None,
        }],
        local: BTreeMap::new(),
    };
    let profile = vec![consensus((0..2).map(|read_index| BaseMeasurement {
        base: Base::T,
        quality: Some(40),
        mapq: Some(60),
        read_index,
    }))?];
    let rate = params.conversion(Some(0.6));
    let expected = (1.0 - rate) * (0.0001_f64 / 3.0) + rate * 0.9999;
    let actual = score(&profile, &template, Strand::OT, 0.6, &params)?.exp();
    assert!((actual - expected).abs() < 1e-12);
    let independent = ((1.0 - rate) * (0.0001 / 3.0) + rate * 0.9999).powi(2);
    assert!((actual - independent).abs() > 0.1);
    let opposite = score(&profile, &template, Strand::OB, 0.6, &params)?.exp();
    assert!((opposite - (0.0001_f64 / 3.0)).abs() < 1e-12);
    Ok(())
}

#[test]
fn identical_mates_are_exactly_one_profile_for_snp_and_indel_likelihoods() -> Result<()> {
    let reference = reference(b"ACTAGTCATG");
    let read = aligned_read(b"ACTGGTCATG", Strand::OT);
    let indel = Indel { start: 1004, inserted: vec![Base::G], deleted: Vec::new() };
    let single = infer_snp_with_indels(
        &reference,
        1003,
        &[],
        std::slice::from_ref(&indel),
        &[Fragment { reads: vec![read.clone()] }],
        &params(true),
        1,
    )?;
    let paired = infer_snp_with_indels(
        &reference,
        1003,
        &[],
        &[indel],
        &[Fragment { reads: vec![read.clone(), read] }],
        &params(true),
        1,
    )?;
    assert_eq!(single.pairs.len(), paired.pairs.len());
    for (a, b) in single.pairs.iter().zip(&paired.pairs) {
        assert!((a.log_likelihood - b.log_likelihood).abs() < 1e-12);
        assert!((*a.posterior - *b.posterior).abs() < 1e-12);
    }
    assert_eq!(paired.focal_observations, 1);
    assert_eq!(paired.focal_read_observations, 2);
    Ok(())
}

#[test]
fn disjoint_mates_keep_their_separate_sequence_evidence() -> Result<()> {
    let reference = reference(b"ACTAGTCATG");
    let first = read(
        b"ACTAG",
        &(1000..1005).map(SequencePosition::reference).collect::<Vec<_>>(),
        Strand::OT,
    );
    let second = read(
        b"TCATG",
        &(1005..1010).map(SequencePosition::reference).collect::<Vec<_>>(),
        Strand::OT,
    );
    let fragments = [
        Fragment { reads: vec![first.clone()] },
        Fragment { reads: vec![second.clone()] },
        Fragment { reads: vec![first, second] },
    ];
    let params = params(true);
    let sites = [Site::new(1002, Base::T, [Base::G])?];
    let prepared =
        Prepared::new(&reference, &sites, None, &fragments.iter().collect::<Vec<_>>(), &params)?;
    let template = prepared.template(
        &reference,
        &Haplotype { alleles: vec![(1002, Base::T)], indel: None },
        &[1002],
        &params,
    )?;
    let values = fragments
        .iter()
        .enumerate()
        .map(|(index, fragment)| {
            prepared
                .fragment_likelihood(index, fragment, &template, 0.0, &params)?
                .ok_or_eyre("Realignment profile")
        })
        .collect::<Result<Vec<_>>>()?;
    assert!(
        (values.first().ok_or_eyre("first")? + values.get(1).ok_or_eyre("second")?
            - values.get(2).ok_or_eyre("paired")?)
        .abs()
            < 1e-12
    );
    Ok(())
}

#[test]
fn conflict_consensus_is_neutral_and_unequal_qualities_reduce_match_confidence() -> Result<()> {
    let params = params(true);
    let template = Template {
        bases: vec![TemplateBase {
            position: SequencePosition::reference(1000),
            base: Base::A,
            cpg: false,
            cpg_key: None,
        }],
        local: BTreeMap::new(),
    };
    for competing_quality in [40, 30] {
        let profile = vec![consensus([
            BaseMeasurement { base: Base::A, quality: Some(40), mapq: Some(60), read_index: 0 },
            BaseMeasurement {
                base: Base::G,
                quality: Some(competing_quality),
                mapq: Some(60),
                read_index: 1,
            },
        ])?];
        let actual = score(&profile, &template, Strand::OT, 0.0, &params)?.exp();
        let expected = if competing_quality == 40 { 0.25 } else { 0.9 };
        assert!((actual - expected).abs() < 1e-12);
        assert!(actual < 0.9999);
    }
    Ok(())
}

#[test]
fn candidate_insertions_and_deletions_change_cpg_context() -> Result<()> {
    let reference = reference(b"ACTGA");
    let insertion = Haplotype {
        alleles: vec![(1001, Base::C)],
        indel: Some(Indel { start: 1002, deleted: Vec::new(), inserted: vec![Base::G] }),
    };
    let deletion = Haplotype {
        alleles: vec![(1001, Base::C)],
        indel: Some(Indel { start: 1002, deleted: vec![Base::T], inserted: Vec::new() }),
    };
    let unedited = Haplotype { alleles: insertion.alleles.clone(), indel: None };
    assert!(cpgs(&reference, &unedited, 1000, 1005)?.is_empty());
    assert_eq!(
        cpgs(&reference, &insertion, 1000, 1005)?.iter().map(|key| key.c).collect::<BTreeSet<_>>(),
        BTreeSet::from([SequencePosition::reference(1001)])
    );
    assert_eq!(
        cpgs(&reference, &deletion, 1000, 1005)?.iter().map(|key| key.c).collect::<BTreeSet<_>>(),
        BTreeSet::from([SequencePosition::reference(1001)])
    );
    Ok(())
}

#[test]
fn indels_are_left_normalized_and_reference_checked() -> Result<()> {
    let reference = reference(b"GATATACT");
    let indel = Indel { start: 1005, deleted: Vec::new(), inserted: vec![Base::A, Base::T] }
        .normalized(&reference)?;
    assert_eq!(indel.start, 1001);
    assert_eq!(indel.inserted, vec![Base::A, Base::T]);
    assert!(
        Indel { start: 1002, deleted: vec![Base::C], inserted: Vec::new() }
            .normalized(&reference)
            .is_err()
    );
    Ok(())
}

#[test]
fn a_nuisance_indel_does_not_change_the_focal_snp_prior_without_reads() -> Result<()> {
    let evidence = infer_snp_with_indels(
        &reference(b"ACTAGTCATG"),
        1003,
        &[],
        &[Indel { start: 1004, deleted: Vec::new(), inserted: vec![Base::C] }],
        &[],
        &params(true),
        1,
    )?;
    assert_eq!(evidence.haplotypes.len(), 8);
    for allele in &evidence.alleles {
        assert!((*allele.posterior - 0.001 * 4.0 / 9.0).abs() < 1e-12);
        assert!(allele.log_bayes_factor.abs() < 1e-12);
    }
    assert!(
        (evidence.pairs.iter().map(|pair| pair.log_prior.exp()).sum::<f64>() - 1.0).abs() < 1e-12
    );
    Ok(())
}

#[test]
fn realignment_removes_a_false_snp_beside_a_misplaced_insertion() -> Result<()> {
    let reference = reference(b"ACTAGTCATG");
    let positions: Vec<_> = (1000..1004)
        .map(SequencePosition::reference)
        .chain([SequencePosition::insertion(1004, 0)])
        .chain((1004..1010).map(SequencePosition::reference))
        .collect();
    let fragments: Vec<_> = std::iter::repeat_n(
        Fragment { reads: vec![read(b"ACTGAGTCATG", &positions, Strand::OT)] },
        18,
    )
    .collect();
    let aligned = infer_snp(
        &reference,
        1003,
        &[],
        &fragments,
        &crate::call::haplotypes::tests::params(true),
        1,
    )?;
    let hmm = infer_snp_with_indels(
        &reference,
        1003,
        &[],
        &[Indel { start: 1003, deleted: Vec::new(), inserted: vec![Base::G] }],
        &fragments,
        &params(true),
        1,
    )?;
    assert!(*aligned.allele(Base::G).expect("allele").posterior > 0.99);
    assert!(*hmm.allele(Base::G).expect("allele").posterior < 0.01, "{hmm:?}");
    let indel_posterior: f64 = hmm
        .pairs
        .iter()
        .filter(|pair| {
            hmm.haplotypes.get(pair.first).expect("hap").indel.is_some()
                || hmm.haplotypes.get(pair.second).expect("hap").indel.is_some()
        })
        .map(|pair| *pair.posterior)
        .sum();
    assert!(indel_posterior > 0.99);
    assert_eq!(hmm.realigned_fragments, 18);
    Ok(())
}

#[test]
fn missing_mapq_and_base_quality_do_not_favor_an_indel_or_snp() -> Result<()> {
    for missing_mapq in [false, true] {
        let mut read = aligned_read(b"ACTAGTCATG", Strand::OT);
        if missing_mapq {
            read.mapq = None;
        } else {
            for observation in &mut read.observations {
                observation.quality = None;
            }
            for base in &mut read.sequence.as_mut().expect("sequence").bases {
                base.quality = None;
            }
        }
        let evidence = infer_snp_with_indels(
            &reference(b"ACTAGTCATG"),
            1003,
            &[],
            &[Indel { start: 1004, deleted: Vec::new(), inserted: vec![Base::C] }],
            &[Fragment { reads: vec![read] }],
            &params(true),
            1,
        )?;
        for allele in &evidence.alleles {
            assert!(allele.log_bayes_factor.abs() < 1e-11);
        }
    }
    Ok(())
}

#[test]
fn discordant_overlapping_cigars_use_an_explicit_aligned_fallback() -> Result<()> {
    let a = aligned_read(b"ACTAGTCATG", Strand::OT);
    let b = read(
        b"ACTGAGTCATG",
        &(1000..1004)
            .map(SequencePosition::reference)
            .chain([SequencePosition::insertion(1004, 0)])
            .chain((1004..1010).map(SequencePosition::reference))
            .collect::<Vec<_>>(),
        Strand::OT,
    );
    let fragments = [Fragment { reads: vec![a, b] }];
    let aligned = infer_snp(
        &reference(b"ACTAGTCATG"),
        1003,
        &[],
        &fragments,
        &crate::call::haplotypes::tests::params(true),
        1,
    )?;
    let hmm = infer_snp(&reference(b"ACTAGTCATG"), 1003, &[], &fragments, &params(true), 1)?;
    assert_eq!(hmm.aligned_fallback_fragments, 1);
    assert_eq!(hmm.realigned_fragments, 0);
    assert_eq!(hmm.focal_observations, 1);
    assert_eq!(hmm.focal_discordant_fragments, 1);
    for allele in &hmm.alleles {
        assert!(
            (allele.log_bayes_factor
                - aligned.allele(allele.base).expect("allele").log_bayes_factor)
                .abs()
                < 1e-12
        );
    }
    Ok(())
}

#[test]
fn taps_and_non_taps_realign_the_same_converted_reads_differently() -> Result<()> {
    let reference = reference(b"ACTCGTATG");
    let fragments: Vec<_> = (0..24)
        .map(|index| Fragment {
            reads: vec![aligned_read(
                if index < 12 { b"ACTTGTATG" } else { b"ACTCGTATG" },
                if index < 12 { Strand::OT } else { Strand::OB },
            )],
        })
        .collect();
    let taps = infer_snp(&reference, 1003, &[], &fragments, &params(false), 1)?;
    let ordinary = infer_snp(&reference, 1003, &[], &fragments, &params(true), 1)?;
    assert!(*taps.allele(Base::T).expect("allele").posterior < 0.01);
    assert!(*ordinary.allele(Base::T).expect("allele").heterozygous > 0.99);
    Ok(())
}

#[test]
fn realignment_removes_a_false_snp_beside_a_misplaced_deletion() -> Result<()> {
    let reference = reference(b"ACTAGTCATG");
    let positions: Vec<_> =
        (1000..1004).chain(1005..1010).map(SequencePosition::reference).collect();
    let fragments: Vec<_> = std::iter::repeat_n(
        Fragment { reads: vec![read(b"ACAGTCATG", &positions, Strand::OT)] },
        18,
    )
    .collect();
    let aligned = infer_snp(
        &reference,
        1003,
        &[],
        &fragments,
        &crate::call::haplotypes::tests::params(true),
        1,
    )?;
    let hmm = infer_snp_with_indels(
        &reference,
        1003,
        &[],
        &[Indel { start: 1002, deleted: vec![Base::T], inserted: Vec::new() }],
        &fragments,
        &params(true),
        1,
    )?;
    assert!(*aligned.allele(Base::G).expect("allele").posterior > 0.99);
    assert!(*hmm.allele(Base::G).expect("allele").posterior < 0.01);
    Ok(())
}

#[test]
fn a_true_snp_is_retained_on_the_haplotype_carrying_an_indel() -> Result<()> {
    let reference = reference(b"ACTAGTCATG");
    let positions: Vec<_> = (1000..1003)
        .map(SequencePosition::reference)
        .chain([SequencePosition::insertion(1003, 0)])
        .chain((1003..1010).map(SequencePosition::reference))
        .collect();
    let mut fragments: Vec<_> = std::iter::repeat_n(
        Fragment { reads: vec![read(b"ACTGAGCCATG", &positions, Strand::OT)] },
        12,
    )
    .collect();
    fragments.extend(std::iter::repeat_n(
        Fragment { reads: vec![aligned_read(b"ACTAGTCATG", Strand::OB)] },
        12,
    ));
    let hmm = infer_snp_with_indels(
        &reference,
        1005,
        &[],
        &[Indel { start: 1003, deleted: Vec::new(), inserted: vec![Base::G] }],
        &fragments,
        &params(true),
        1,
    )?;
    assert!(*hmm.allele(Base::C).expect("allele").heterozygous > 0.999);
    let best = hmm
        .pairs
        .iter()
        .max_by(|a, b| a.posterior.partial_cmp(&b.posterior).expect("finite probability"))
        .expect("pair");
    let a = hmm.haplotypes.get(best.first).expect("haplotype");
    let b = hmm.haplotypes.get(best.second).expect("haplotype");
    assert!(*best.posterior > 0.99);
    for haplotype in [a, b] {
        assert_eq!(haplotype.indel.is_some(), haplotype.allele_at(1005) == Some(Base::C));
    }
    Ok(())
}

#[test]
fn reads_spanning_a_deletion_at_the_target_are_still_scored() -> Result<()> {
    let positions: Vec<_> =
        (1000..1003).chain(1004..1010).map(SequencePosition::reference).collect();
    let fragments = [Fragment { reads: vec![read(b"ACTGTCATG", &positions, Strand::OT)] }];
    let hmm = infer_snp(&reference(b"ACTAGTCATG"), 1003, &[], &fragments, &params(true), 1)?;
    assert_eq!(hmm.focal_observations, 0);
    assert_eq!(hmm.focal_fragments, 1);
    assert_eq!(hmm.realigned_fragments, 1);
    Ok(())
}

#[test]
fn cpg_annotation_uses_context_outside_the_realignment_interval() -> Result<()> {
    let mut params = params(false);
    params.realignment_flank = 2;
    for (bases, focal, edge) in [(b"ACGTAAA", 1004, 0usize), (b"ATCCGTA", 1001, 3)] {
        let reference = reference(bases);
        let site = Site::new(focal, reference.base_at(focal).expect("reference"), Base::KNOWN)?;
        let haplotype = Haplotype { alleles: vec![(focal, site.reference)], indel: None };
        let prepared = Prepared::new(&reference, &[site], None, &[], &params)?;
        let template = prepared.template(&reference, &haplotype, &[focal], &params)?;
        assert!(template.bases.get(edge).expect("edge").cpg);
        assert!(template.bases.get(edge).expect("edge").cpg_key.is_some());
    }
    Ok(())
}

#[test]
fn internal_inserted_cpgs_have_distinct_rate_keys() -> Result<()> {
    let reference = reference(b"ACTAGTCATG");
    let haplotype = Haplotype {
        alleles: vec![(1003, Base::A)],
        indel: Some(Indel {
            start: 1003,
            deleted: Vec::new(),
            inserted: vec![Base::C, Base::G, Base::C, Base::G],
        }),
    };
    let keys = cpgs(&reference, &haplotype, 1000, 1010)?;
    assert_eq!(keys.len(), 2);
    assert!(keys.iter().all(|key| !key.c.is_reference() && !key.g.is_reference()));
    Ok(())
}

#[test]
fn aligned_fallback_cpg_context_uses_candidate_indel_adjacency() -> Result<()> {
    let context = reference(b"ACTGA");
    for indel in [
        Indel { start: 1002, deleted: vec![Base::T], inserted: Vec::new() },
        Indel { start: 1002, deleted: Vec::new(), inserted: vec![Base::G] },
    ] {
        let haplotype = Haplotype { alleles: vec![(1001, Base::C)], indel: Some(indel) };
        assert_eq!(super::super::likelihood::cpg_anchor(&context, &haplotype, 1001), Some(1001));
    }
    let disrupted = Haplotype {
        alleles: vec![(1001, Base::C)],
        indel: Some(Indel { start: 1002, deleted: Vec::new(), inserted: vec![Base::T] }),
    };
    assert_eq!(super::super::likelihood::cpg_anchor(&reference(b"ACGA"), &disrupted, 1001), None);
    Ok(())
}

#[test]
fn inconsistent_sequence_and_aligned_observations_are_rejected() -> Result<()> {
    let mut read = aligned_read(b"ACTAGTCATG", Strand::OT);
    read.observations.get_mut(3).expect("observation").base = Base::G;
    assert!(
        infer_snp(
            &reference(b"ACTAGTCATG"),
            1003,
            &[],
            &[Fragment { reads: vec![read] }],
            &params(true),
            1
        )
        .is_err()
    );
    Ok(())
}

#[test]
fn an_unused_mate_cannot_change_orientation_or_force_aligned_fallback() -> Result<()> {
    let reference = reference(b"ACTCGTATG");
    let active = aligned_read(b"ACTTGTATG", Strand::OT);
    let single: Vec<_> =
        std::iter::repeat_n(Fragment { reads: vec![active.clone()] }, 18).collect();
    let expected = infer_snp(&reference, 1003, &[], &single, &params(false), 1)?;
    for reference_skip in [false, true] {
        let mut unused = read(
            b"ACTAGTCATG",
            &(1020..1030).map(SequencePosition::reference).collect::<Vec<_>>(),
            Strand::OB,
        );
        unused.sequence.as_mut().expect("sequence").reference_skip = reference_skip;
        let paired: Vec<_> =
            std::iter::repeat_n(Fragment { reads: vec![active.clone(), unused] }, 18).collect();
        let actual = infer_snp(&reference, 1003, &[], &paired, &params(false), 1)?;
        assert_eq!(actual.aligned_fallback_fragments, 0);
        for allele in &actual.alleles {
            assert!(
                (allele.log_bayes_factor
                    - expected.allele(allele.base).expect("allele").log_bayes_factor)
                    .abs()
                    < 1e-12
            );
        }
    }
    Ok(())
}

#[test]
fn read_local_windows_keep_the_molecular_observation_model_normalized() -> Result<()> {
    let reference = Reference { start: 1000, bases: vec![Base::A; 500] };
    let sites = vec![Site::new(1100, Base::A, [Base::C])?, Site::new(1400, Base::A, [Base::G])?];
    let haplotype = Haplotype { alleles: vec![(1100, Base::A), (1400, Base::G)], indel: None };
    let params = params(true);
    let mut total = 0.0;
    for a in Base::KNOWN {
        for b in Base::KNOWN {
            let reads = [(1100, a), (1400, b)].map(|(pos, base)| {
                let observation = Observation { pos, base, quality: Some(40) };
                ReadEvidence {
                    annotations: None,
                    mapq: Some(60),
                    strand: Strand::OT,
                    observations: vec![observation],
                    sequence: Some(SequenceEvidence {
                        bases: vec![SequenceObservation {
                            position: SequencePosition::reference(pos),
                            base,
                            quality: Some(40),
                        }],
                        reference_skip: false,
                        indels: Vec::new(),
                    }),
                }
            });
            let fragment = Fragment { reads: reads.to_vec() };
            let prepared =
                Prepared::new(&reference, &sites, None, &[&fragment], &params)?.localize();
            let template = prepared.template(&reference, &haplotype, &[1100, 1400], &params)?;
            assert_eq!(template.local.len(), 2);
            assert!(template.local.values().all(|template| template.bases.len() <= 25));
            total += prepared
                .fragment_likelihood(0, &fragment, &template, 0.5, &params)?
                .expect("realigned")
                .exp();
        }
    }
    assert!((total - 1.0).abs() < 1e-10, "{total}");
    Ok(())
}

#[test]
fn scaled_forward_matches_the_log_domain_oracle_on_short_and_long_profiles() -> Result<()> {
    let params = params(true);
    for rows in [1, 7, 40, 150, 500] {
        for columns in [1, 9, 63, 180] {
            let insertion = |row| if row % 9 == 0 { 0.0 } else { 0.25_f64.ln() };
            let emission =
                |row, column| (0.001 + ((row * 7 + column * 11) % 19) as f64 * 0.045).ln();
            let original = forward_log(rows, columns, &params, insertion, emission)?;
            let scaled = forward(rows, columns, &params, insertion, emission)?;
            assert!(
                (original - scaled).abs() < 1e-9,
                "{rows}/{columns}: {original} versus {scaled}"
            );
        }
    }
    Ok(())
}

#[test]
fn phase_likelihood_reuse_ignores_unsequenced_mate_gaps_but_preserves_cpg_partners() -> Result<()> {
    let reference = Reference { start: 1000, bases: vec![Base::A; 500] };
    let sites = vec![Site::new(1100, Base::A, [Base::C])?, Site::new(1300, Base::A, [Base::G])?];
    let reads = [1100, 1300].map(|pos| ReadEvidence {
        annotations: None,
        mapq: Some(60),
        strand: Strand::OT,
        observations: vec![Observation { pos, base: Base::A, quality: Some(40) }],
        sequence: Some(SequenceEvidence {
            bases: vec![SequenceObservation {
                position: SequencePosition::reference(pos),
                base: Base::A,
                quality: Some(40),
            }],
            reference_skip: false,
            indels: Vec::new(),
        }),
    });
    let fragment = Fragment { reads: reads.to_vec() };
    let params = params(true);
    let prepared = Prepared::new(&reference, &sites, None, &[&fragment], &params)?.localize();
    let first =
        Haplotype { alleles: vec![(1100, Base::A), (1200, Base::A), (1300, Base::G)], indel: None };
    let mut unobserved = first.clone();
    unobserved.alleles.get_mut(1).expect("gap SNP").1 = Base::G;
    assert_eq!(prepared.signature(&first, &params), prepared.signature(&unobserved, &params));
    let a = prepared.template(&reference, &first, &[1100, 1300], &params)?;
    let b = prepared.template(&reference, &unobserved, &[1100, 1300], &params)?;
    assert_eq!(
        prepared.fragment_likelihood(0, &fragment, &a, 0.5, &params)?,
        prepared.fragment_likelihood(0, &fragment, &b, 0.5, &params)?
    );
    let mut partner = first.clone();
    partner.alleles.push((1113, Base::G));
    partner.alleles.sort_by_key(|(pos, _)| *pos);
    assert_ne!(prepared.signature(&first, &params), prepared.signature(&partner, &params));
    Ok(())
}

/// The mixed emission as production scored it before the dot-product weights
/// and the lane-batched kernel: every branch of the closure, cell by cell.
fn reference_score_with_model(
    profile: &Profile,
    template: &Template,
    orientation: Strand,
    scoring: &super::super::methylation::Scoring<'_>,
    params: &LikelihoodParams,
    bad_reads: usize,
) -> Result<f64> {
    if profile.is_empty() {
        return Ok(0.0);
    }
    if profile.iter().all(|base| base.source.is_none_or(|index| bad_reads & (1 << index) != 0)) {
        return Ok(profile.iter().filter(|base| base.base.known_index().is_some()).count() as f64
            * 0.25_f64.ln());
    }
    let insertions: Vec<_> = profile
        .iter()
        .map(|observed| if observed.base.known_index().is_some() { 0.25 } else { 1.0 })
        .collect();
    let emissions: Vec<[f64; 4]> = profile
        .iter()
        .map(|observed| {
            Base::KNOWN.map(|base| {
                if observed.source.is_some_and(|index| bad_reads & (1 << index) != 0) {
                    0.25
                } else {
                    sequencing_emission(observed.base, base, observed.quality)
                }
            })
        })
        .collect();
    let conversions: Vec<_> = template
        .bases
        .iter()
        .map(|base| {
            let side = [CpgSide::C, CpgSide::G].into_iter().find(|side| {
                side.unmod_base() == base.base && side.strand() == orientation && !params.non_taps
            });
            let rate = params.conversion(base.cpg.then_some(scoring.rate(base.cpg_key.as_ref())));
            (
                base.base.known_index(),
                side.and_then(|side| side.mod_base().known_index()),
                1.0 - rate,
                rate,
                if params.non_taps {
                    None
                } else {
                    scoring.neutral_bases(base.position, orientation, base.base)
                },
            )
        })
        .collect();
    let insertion = |i: usize| insertions.get(i).copied().unwrap_or(0.0);
    let emission = |i: usize, j: usize| {
        let Some(&(ordinary, modified, unconverted_rate, converted_rate, masked)) =
            conversions.get(j)
        else {
            return 0.0;
        };
        let Some(ordinary) = ordinary else { return insertion(i) };
        let Some(values) = emissions.get(i) else { return 0.0 };
        if let Some((a, b)) = masked {
            let Some(a) = a.known_index().and_then(|index| values.get(index)) else {
                return 0.0;
            };
            let Some(b) = b.known_index().and_then(|index| values.get(index)) else {
                return 0.0;
            };
            return (*a + *b) * 0.5;
        }
        let ordinary = values.get(ordinary).copied().unwrap_or(0.0);
        let Some(modified) = modified.and_then(|index| values.get(index).copied()) else {
            return ordinary;
        };
        unconverted_rate * ordinary + converted_rate * modified
    };
    let value =
        forward_scaled_reference(profile.len(), template.bases.len(), params, insertion, emission)?;
    let mass = forward(profile.len(), template.bases.len(), params, |_| 0.0, |_, _| 0.0)?;
    Ok(value - mass)
}

#[test]
#[expect(clippy::cast_possible_truncation, reason = "Test code")]
fn dot_product_weights_and_the_batched_kernel_reproduce_the_mixed_emission_bit_for_bit()
-> Result<()> {
    use super::super::methylation::{Context, Scoring};
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> f64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
        fn below(&mut self, probability: f64) -> bool {
            self.next() < probability
        }
        fn base(&mut self, unknown: f64) -> Base {
            if self.below(unknown) {
                Base::Unknown
            } else {
                Base::KNOWN.get(((self.next() * 4.0) as usize).min(3)).copied().unwrap_or(Base::A)
            }
        }
    }
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    let mut bases: Vec<_> = (0..64).map(|_| rng.base(0.0)).collect();
    for start in [10, 30, 50] {
        *bases.get_mut(start).expect("CpG C") = Base::C;
        *bases.get_mut(start + 1).expect("CpG G") = Base::G;
    }
    let reference = Reference { start: 1000, bases };
    let haplotypes = [Haplotype { alleles: Vec::new(), indel: None }];
    let context = Context::new(&reference, &haplotypes, &[1005, 1040], Some(1000..1064))?;
    let mut compared = 0;
    for trial in 0..160 {
        let params = params(trial % 5 == 0);
        let rows = 1 + (trial * 7) % 13;
        let columns = 1 + (trial * 11) % 17;
        let profile: Profile = (0..rows)
            .map(|_| {
                let base = rng.base(0.15);
                MolecularBase {
                    base,
                    quality: if rng.below(0.2) { None } else { Some((rng.next() * 45.0) as u8) },
                    source: base.known_index().map(|_| usize::from(rng.below(0.5))),
                    measurements: 1,
                    discordant: false,
                }
            })
            .collect();
        let start = 1000 + (rng.next() * 40.0) as u32;
        let template = Template {
            bases: (0..columns)
                .map(|offset| {
                    let pos = start + u32::try_from(offset).expect("short template");
                    let base = if rng.below(0.5) {
                        reference.base_at(pos).unwrap_or(Base::A)
                    } else {
                        rng.base(0.1)
                    };
                    TemplateBase {
                        position: SequencePosition::reference(pos),
                        base,
                        cpg: rng.below(0.3),
                        cpg_key: None,
                    }
                })
                .collect(),
            local: BTreeMap::new(),
        };
        for orientation in [Strand::OT, Strand::OB] {
            for bad_reads in 0..4 {
                let fixed = Scoring::Fixed(rng.next());
                for scoring in [&fixed, &Scoring::Neutral(&context), &Scoring::Ordinary] {
                    let expected = reference_score_with_model(
                        &profile,
                        &template,
                        orientation,
                        scoring,
                        &params,
                        bad_reads,
                    )?;
                    let actual = score_with_model(
                        &profile,
                        &template,
                        orientation,
                        scoring,
                        &params,
                        bad_reads,
                    )?;
                    assert_eq!(
                        expected.to_bits(),
                        actual.to_bits(),
                        "trial {trial} {orientation:?} mask {bad_reads}: {expected} versus {actual}"
                    );
                    compared += 1;
                }
            }
        }
    }
    assert_eq!(compared, 160 * 2 * 4 * 3);
    Ok(())
}

#[test]
fn bounded_mapping_hypotheses_are_skipped_without_changing_the_likelihood() -> Result<()> {
    use super::super::methylation::Scoring;
    let params = params(false);
    let mut state = 0xdead_beef_u32;
    let sequence: Vec<u8> = (0..200)
        .map(|_| {
            state = state.wrapping_mul(1664525).wrapping_add(1013904223);
            b"ACGT".get((state >> 30) as usize).copied().expect("DNA")
        })
        .collect();
    let reference = reference(&sequence);
    // Sites at both ends keep the realignment window over each mate's own bases.
    let sites = [
        Site::new(1030, reference.base_at(1030).ok_or_eyre("site")?, Base::KNOWN)?,
        Site::new(1150, reference.base_at(1150).ok_or_eyre("site")?, Base::KNOWN)?,
    ];
    let far = reference.base_at(1150).ok_or_eyre("far site")?;
    let positions =
        |range: std::ops::Range<u32>| range.map(SequencePosition::reference).collect::<Vec<_>>();
    // A clean overlapping pair can only lose ~1e-6 of its weight to a
    // mismapped mate. A mate full of mismatches makes its own hypothesis the
    // likely one, and drags the all-mapped score down to where the bound,
    // which knows nothing of the template, can exclude neither.
    for (mismatches, expected_skips_per_plan) in [(0, 2), (40, 0)] {
        let first = sequence.get(20..120).ok_or_eyre("first mate")?.to_vec();
        let mut second = sequence.get(60..160).ok_or_eyre("second mate")?.to_vec();
        for base in second.iter_mut().step_by(2).take(mismatches) {
            *base = match *base {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
        let fragment = Fragment {
            reads: vec![
                read(&first, &positions(1020..1120), Strand::OT),
                read(&second, &positions(1060..1160), Strand::OT),
            ],
        };
        let prepared = Prepared::new(&reference, &sites, None, &[&fragment], &params)?.localize();
        let mut lazy = Scorer::new(&params);
        let mut exhaustive = Scorer::new(&params).exhaustive();
        for base in Base::KNOWN {
            let haplotype = Haplotype { alleles: vec![(1030, base), (1150, far)], indel: None };
            let template = prepared.template(&reference, &haplotype, &[1030, 1150], &params)?;
            let plan = prepared
                .plan(0, &fragment, &template, &params, &Scoring::Fixed(0.3))?
                .ok_or_eyre("realigned fragment")?;
            let [skipping] = lazy.run(std::slice::from_ref(&plan), true)?[..] else {
                color_eyre::eyre::bail!("one value per plan")
            };
            let [complete] = exhaustive.run(std::slice::from_ref(&plan), true)?[..] else {
                color_eyre::eyre::bail!("one value per plan")
            };
            assert_eq!(skipping.to_bits(), complete.to_bits(), "{mismatches} mismatches, {base:?}");
        }
        assert_eq!(lazy.skipped, expected_skips_per_plan * Base::KNOWN.len());
        assert_eq!(exhaustive.skipped, 0);
    }
    Ok(())
}
