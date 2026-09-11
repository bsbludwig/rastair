//! Tests for genotyping with ML scores
use crate::{call::tests::utils::*, call::variant_calling::GenotypeTag, pileups, vcf_assert};
use seqair_types::{Base::*, Probability};
use std::num::NonZeroU8;

#[test]
fn test_a_to_t_high_ml_score() -> Result<()> {
    // A→T transition with high ML score should be genotyped
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    let records = reprocess(records)?;

    // 3 T reads, 1 A read -> heterozygous (0/1)
    let expected_vcf = vcf_assert![
        (A T) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_a_to_g_high_ml_score() -> Result<()> {
    // A→G transition with high ML score should be genotyped
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ G ] OT,
        [ G ] OT,
        [ G ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    // 3 G reads, 1 A read -> heterozygous (0/1)
    let expected_vcf = vcf_assert![
        (A G) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_a_to_c_high_ml_score() -> Result<()> {
    // A→C transition with high ML score should be genotyped
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ C ] OT,
        [ C ] OT,
        [ C ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C);
    let records = reprocess(records)?;

    // 3 C reads, 1 A read -> heterozygous (0/1)
    let expected_vcf = vcf_assert![
        (A C) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_t_to_a_high_ml_score() -> Result<()> {
    // T→A transition with high ML score should be genotyped
    let (segment, pileups) = pileups!(
        [ T ] Ref,
        [ A ] OT,
        [ A ] OT,
        [ A ] OB,
        [ T ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], A);
    let records = reprocess(records)?;

    // 3 A reads, 1 T read -> heterozygous (0/1)
    let expected_vcf = vcf_assert![
        (T A) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_t_to_g_high_ml_score() -> Result<()> {
    // T→G transition with high ML score should be genotyped
    let (segment, pileups) = pileups!(
        [ T ] Ref,
        [ G ] OT,
        [ G ] OT,
        [ G ] OB,
        [ T ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    // 3 G reads, 1 T read -> heterozygous (0/1)
    let expected_vcf = vcf_assert![
        (T G) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_t_to_c_high_ml_score() -> Result<()> {
    // T→C transition with high ML score should be genotyped
    let (segment, pileups) = pileups!(
        [ T ] Ref,
        [ C ] OT,
        [ C ] OT,
        [ C ] OB,
        [ T ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C);
    let records = reprocess(records)?;

    // 3 C reads, 1 T read -> heterozygous (0/1)
    let expected_vcf = vcf_assert![
        (T C) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_c_to_a_high_ml_score() -> Result<()> {
    // C→A transition (non-methylation-confounded) with high ML score
    let (segment, pileups) = pileups!(
        [ C ] Ref,
        [ A ] OT,
        [ A ] OT,
        [ A ] OB,
        [ C ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], A);
    let records = reprocess(records)?;

    // 3 A reads, 1 C read -> heterozygous (0/1)
    let expected_vcf = vcf_assert![
        (C A) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_c_to_g_high_ml_score() -> Result<()> {
    // C→G transition (non-methylation-confounded) with high ML score
    let (segment, pileups) = pileups!(
        [ C ] Ref,
        [ G ] OT,
        [ G ] OT,
        [ G ] OB,
        [ C ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    // 3 G reads, 1 C read -> heterozygous (0/1)
    let expected_vcf = vcf_assert![
        (C G) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_c_to_t_outside_cpg_counts_both_strands() -> Result<()> {
    // C→T outside CpG should count both strands for genotyping.
    let (segment, pileups) = pileups!(
        [ C ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ C ] OB,
        [ C ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    let records = reprocess(records)?;

    // 2 T reads, 2 C reads -> heterozygous (0/1)
    let expected_vcf = vcf_assert![
        (C T) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_g_to_t_high_ml_score() -> Result<()> {
    // G→T transition (non-methylation-confounded) with high ML score
    let (segment, pileups) = pileups!(
        [ G ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    let records = reprocess(records)?;

    // 3 T reads, 1 G read -> heterozygous (0/1)
    let expected_vcf = vcf_assert![
        (G T) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_g_to_c_high_ml_score() -> Result<()> {
    // G→C transition (non-methylation-confounded) with high ML score
    let (segment, pileups) = pileups!(
        [ G ] Ref,
        [ C ] OT,
        [ C ] OT,
        [ C ] OB,
        [ G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C);
    let records = reprocess(records)?;

    // 3 C reads, 1 G read -> heterozygous (0/1)
    let expected_vcf = vcf_assert![
        (G C) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_multi_allelic_site_picks_highest_ml() -> Result<()> {
    // Multi-allelic site: A with both T and G alts
    // Should genotype based on highest ML score
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ G ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T); // ML = 1.0

    // Set G to have lower ML but still passing
    let g_filters = records[0].alt_filters_mut(G).unwrap();
    g_filters.ml = Some(Probability::new_panicky(0.85));
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (A T,G) PASS,  // Both alts passing, combined in one row
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_all_alts_below_ml_threshold() -> Result<()> {
    // All alts have ML scores below threshold - should call HomRef (0/0)
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ G ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_fail(&mut records[0], T); // ML = 0.0
    set_fail(&mut records[0], G); // ML = 0.0
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (A T) FAIL GT="0/0",
        (A G) FAIL GT="0/0",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_het_vs_hom_alt_genotyping_for_a_base() -> Result<()> {
    // A→T with reads suggesting heterozygous genotype (balanced ref/alt)
    // Note: With exactly balanced reads (2 ref, 2 alt), the binomial model
    // may return HomRef (0/0) depending on the error model parameters
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ A ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    let records = reprocess(records)?;

    // 2 T reads, 2 A reads -> with perfectly balanced reads
    let expected_vcf = vcf_assert![
        (A T) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_het_genotyping_with_unbalanced_reads() -> Result<()> {
    // A→T with unbalanced reads (more alt than ref) should produce het (0/1)
    // This demonstrates that het genotyping works when reads clearly favor alt
    // but not overwhelmingly enough for hom-alt
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ A ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    let records = reprocess(records)?;

    // 3 T reads, 2 A reads -> heterozygous (0/1)
    let expected_vcf = vcf_assert![
        (A T) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_hom_alt_genotyping_for_t_base() -> Result<()> {
    // T→G with mostly alt reads suggesting homozygous alt genotype (1/1)
    let (segment, pileups) = pileups!(
        [ T ] Ref,
        [ G ] OT,
        [ G ] OT,
        [ G ] OB,
        [ G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    // 4 G reads, 0 T reads -> homozygous alt (1/1)
    let expected_vcf = vcf_assert![
        (T G) PASS GT="1/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_compound_het_with_balanced_alts() -> Result<()> {
    // A reference with two alt alleles (T and G), both passing ML threshold
    // with roughly balanced read support should produce compound het (1/2)
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ G ] OT,
        [ G ] OT,
        [ G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    // 3 T reads, 3 G reads, 0 A reads -> compound heterozygous (1/2)
    let expected_vcf = vcf_assert![
        (A T,G) PASS GT="1/2",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_compound_het_with_slightly_unbalanced_alts() -> Result<()> {
    // Two alts passing ML threshold with slightly unbalanced reads
    // Should still call compound het if both have >20% of alt reads
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ T ] OB,
        [ G ] OT,
        [ G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    // 4 T reads (67%), 2 G reads (33%), 0 A reads
    // Both alts have >20% -> should call compound het (1/2)
    let expected_vcf = vcf_assert![
        (A T,G) PASS GT="1/2",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_two_alts_passing_but_one_dominant() -> Result<()> {
    // Two alts pass ML threshold but one has <20% of alt reads
    // Should fall back to single alt genotyping with the dominant alt
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ T ] OB,
        [ T ] OB,
        [ G ] OT,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    // 5 T reads (83%), 1 G read (17%), 1 A read
    // G has <20% of alt reads -> should call het with T only (0/1)
    let expected_vcf = vcf_assert![
        (A T,G) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_compound_het_vs_single_het_with_ref() -> Result<()> {
    // Two alts pass ML, but significant ref reads suggest 0/1 instead of 1/2
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ G ] OT,
        [ G ] OB,
        [ A ] OT,
        [ A ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    // 3 T reads, 2 G reads, 3 A reads
    // Both alts have reasonable support but ref reads suggest het not compound het
    // The likelihood model should determine 0/1 or 0/2 is more likely than 1/2
    // With 3 T, 2 G, 3 A reads, 0/1 (ref + T) is most likely
    let expected_vcf = vcf_assert![
        (A T,G) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_compound_het_for_c_base_with_non_methylation_alts() -> Result<()> {
    // C reference with A and G alts (not methylation-confounded)
    // Should use all reads, not strand-specific
    let (segment, pileups) = pileups!(
        [ C ] Ref,
        [ A ] OT,
        [ A ] OT,
        [ A ] OB,
        [ G ] OT,
        [ G ] OB,
        [ G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], A);
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    // 3 A reads, 3 G reads, 0 C reads -> compound heterozygous (1/2)
    let expected_vcf = vcf_assert![
        (C A,G) PASS GT="1/2",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_three_alts_passing_uses_top_two() -> Result<()> {
    // Three alts pass ML threshold, should consider top 2 by ML score
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ G ] OT,
        [ G ] OB,
        [ C ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;

    // Set T with highest ML, G with medium ML, C with lowest ML
    set_pass(&mut records[0], T);
    records[0].alt_filters_mut(T).unwrap().ml = Some(Probability::new_panicky(0.95));
    set_pass(&mut records[0], G);
    records[0].alt_filters_mut(G).unwrap().ml = Some(Probability::new_panicky(0.90));
    set_pass(&mut records[0], C);
    records[0].alt_filters_mut(C).unwrap().ml = Some(Probability::new_panicky(0.85));
    let records = reprocess(records)?;

    // Should genotype based on top 2 (T and G)
    // 2 T reads, 2 G reads, 1 C read -> compound het with T and G
    let expected_vcf = vcf_assert![
        (A T,G,C) PASS GT="1/2",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// Test that de-novo CpG genotyping ignores REF reads from converted strand.
#[test]
fn test_denovo_cpg_genotyping_ignores_converted_strand_ref_reads() -> Result<()> {
    // Position 1: A with next base G (A→C creates de-novo CpG)
    // Position 2: G (partner of de-novo CpG)
    let (segment, pileups) = pileups!(
        [ A G ] Ref,
        // OT strand: A (ref) and T (methylated C)
        [ A G ] OT,
        [ A G ] OT,
        [ T G ] OT,  // methylated
        [ T G ] OT,  // methylated
        // OB strand: C (alt) only - no ref reads
        [ C A ] OB,
        [ C A ] OB,
        [ C A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C); // C is the real variant creating de-novo CpG
    set_fail(&mut records[0], T); // T is methylation evidence, not a real variant
    set_fail(&mut records[1], A); // A on second position is methylation evidence
    let records = reprocess(records)?;

    // Genotyping for the C alt should only use OB strand
    let expected_vcf = vcf_assert![
        (A C) PASS GT="0/1",
        (A T) FAIL,
        (G .) PASS,  // partner position of de-novo CpG
        (G A) FAIL,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// The writer must report the genotype the pipeline stored, not one it derives
/// again for itself. `--error-model` only reaches the stored estimate, so a
/// writer that re-estimates with the default model contradicts the BED output
/// and the methylation calls, which both read `pos_metrics.extended.genotype`.
#[test]
fn vcf_genotype_is_the_stored_estimate() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    let mut records = reprocess(records)?;

    // Stand-in for what another error model concludes from the same counts.
    let stored = records[0].pos_metrics.extended.genotype.as_mut().expect("genotype");
    stored.genotype = GenotypeTag::hom_alt(NonZeroU8::new(1).expect("1 > 0"));

    let expected_vcf = vcf_assert![
        (A T) PASS GT="1/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}
