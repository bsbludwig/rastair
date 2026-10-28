use crate::{
    call::{
        tests::pileup::{SimpleRead, SimpleReads},
        tests::utils::*,
    },
    pileups, vcf_assert,
};
use seqair_types::Base::*;
use std::sync::Arc;

#[test]
fn test_cpg_context() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ C G ] OT,
        [ C A ] OB,
        [ C G ] OB,
    );

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=2./3.,
        (C T) FAIL,
        (G .) PASS M5mC=1./2.,
        (G A) FAIL,
    ];

    let records = test_call(segment, pileups, RecordFilters::all())?;
    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_mixed_methylation_and_real_variants() -> Result<()> {
    // Test C with both methylation evidence (C->T low ML) and real variant (C->A high ML)
    // Assumption: When we have both types, we output ref->. plus both alt rows
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ A G ] OB,
        [ A G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_fail(&mut records[0], T); // Low ML - methylation evidence
    set_pass(&mut records[0], A); // High ML - real variant
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C A) PASS M5mC=1. GT="0/1",  // T reads (methylation evidence) count as C allele → genotype C/A, beta 1.0
        (C T) FAIL,  // methylation evidence T fails as separate record
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_non_methylation_transitions() -> Result<()> {
    // Test transitions that are NOT methylation-related (C->A, C->G, G->C, G->T)
    // Assumption: These should never produce ref->. rows based on ML score alone
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ A C ] OT,
        [ A C ] OT,
        [ A C ] OB,
        [ A C ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], A);
    set_pass(&mut records[1], C);
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C A) PASS,
        (G C) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_multiple_methylation_transitions_same_position() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OB,
        [ T G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_fail(&mut records[0], T); // Force low ML
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=1.,  // genotype C/C, all reads show T (fully methylated)
        (C T) FAIL,  // the methylation transition
        (G .) PASS M5mC=None,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_reverse_strand_methylation() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ C G ] OT,
        [ C G ] OT,
        [ C A ] OB,
        [ C A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_fail(&mut records[1], A); // Low ML - methylation evidence on reverse strand
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=0.,
        (G .) PASS M5mC=1.,  // methylation evidence
        (G A) FAIL,  // low confidence methylation transition
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_empty_alts_only_produces_one_row() -> Result<()> {
    // Test that positions with no alts produce exactly one ref->. row
    // Assumption: No alts means reference matches, output single ref->. PASS row
    // Using proper CpG context (C followed by G) to ensure record passes filter
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ C G ] OT,
        [ C G ] OT,
        [ C G ] OB,
        [ C G ] OB,
    );

    let records = test_call(segment, pileups, RecordFilters::all())?;

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=0.,  // Single row for matching reference
        (G .) PASS M5mC=0.,  // Single row for matching reference
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_all_methylation_transitions_failing() -> Result<()> {
    // Test CpG site where both positions have failing methylation transitions
    // Assumption: Both C->T and G->A with low ML produce ref->. rows
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ C A ] OB,
        [ C A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_fail(&mut records[0], T);
    set_fail(&mut records[1], A);
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=1.,  // genotype C/C, all reads show T (fully methylated)
        (C T) FAIL,  // methylation transition
        (G .) PASS M5mC=1.,  // genotype G/G, all reads show A (fully methylated)
        (G A) FAIL,  // methylation transition
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_all_methylation_transitions_passing_as_variants() -> Result<()> {
    // Test CpG site where both positions have failing methylation transitions
    // Assumption: Both C->T and G->A with low ML produce ref->. rows
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T A ] OT,
        [ T A ] OT,
        [ T A ] OB,
        [ T A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    set_pass(&mut records[1], A);
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C T) PASS,
        (G A) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_non_covered_cpg_filtering() -> Result<()> {
    let (segment, pileups) = create_pileups(vec![C, G], vec![]);

    let records = test_call(segment.clone(), pileups.clone(), RecordFilters::cpgs())?;
    let vcf_records = metrics_to_vcf(&records, RecordFilters::cpgs())?;
    assert!(vcf_records.is_empty(), "Expected no VCF records for non-covered CpG site");

    let records = test_call(segment, pileups, RecordFilters::all())?;
    let expected_vcf = vcf_assert![
        (C .) PASS,
        (G .) PASS,
    ];
    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_methylation_transition_is_called_by_its_ml_score() -> Result<()> {
    // Test C->T with and without an ML score: the score decides
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OB,
        [ T G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T); // High ML score
    let mut records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C T) PASS,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    // Now, without ML
    records[0].find_alt_mut(T).unwrap().verdict = crate::metrics::AltVerdict::Unscored;
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS,
        (C T) FAIL,  // FAIL because ML is None
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_single_end_like_methylation_calling() -> Result<()> {
    // Single-end input has no pair tags; all reads are treated as first-in-pair.
    // OT/OB comes from alignment direction only.
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ C G ] OT,
        [ C G ] OT,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_fail(&mut records[0], T);
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=1./2.,
        (C T) FAIL,
        (G .) PASS M5mC=None,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_c_with_t_and_a_both_low_ml() -> Result<()> {
    // Test C with both T and A alts, both failing ML threshold
    // Assumption: C->T is methylation transition, C->A is not
    // Expected: ref->. row for methylation evidence, then both alts fail separately
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ A G ] OB,
        [ A G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_fail(&mut records[0], T); // Low ML - methylation evidence
    set_fail(&mut records[0], A); // Low ML - but not methylation transition
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=1.,  // methylation evidence from C->T
        (C T) FAIL,
        (C A) FAIL,
        (G .) PASS M5mC=None,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_g_with_a_and_c_alts() -> Result<()> {
    // Test G->A methylation evidence + G->C in a GC context (reverse CpG)
    // Context: G followed by C (GC = CpG on reverse strand), where G has alts A and C
    // Assumption: Multiple biological signals at one position should all be reported
    let (segment, pileups) = pileups!(
        [ G C ] Ref,
        [ A C ] OT,
        [ A C ] OT,
        [ C C ] OB,
        [ C C ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_fail(&mut records[0], A); // G->A low ML - methylation evidence
    set_pass(&mut records[0], C); // G->C high ML - non-methylation variant

    // Reprocess to recalculate genotypes with the modified ML scores
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (G C) PASS,  // All alts combined: A fails, C passes (any passing = PASS)
        (G A) FAIL,  // All alts combined: A fails, C passes (any passing = PASS)
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_cpg_islands_multiple_positions() -> Result<()> {
    // Test CGCGCG reference with various methylation patterns
    // Assumption: Multiple consecutive CpG sites should each be handled independently
    // Pattern: C1 methylated, G1 not, C2 not, G2 methylated, C3 methylated, G3 not
    let (segment, pileups) = pileups!(
        [ C G C G C G ] Ref,
        [ T G C A T G ] OT,
        [ T G C A T G ] OT,
        [ T G C A T G ] OB,
        [ T G C A T G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    // C1->T: low ML (methylation)
    set_fail(&mut records[0], T);
    // G1: matches ref (no alt)
    // C2: matches ref (no alt)
    // G2->A: low ML (methylation)
    set_fail(&mut records[3], A);
    // C3->T: low ML (methylation)
    set_fail(&mut records[4], T);
    // G3: matches ref (no alt)
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=1.,  // C1: genotype C/C, all reads show T (fully methylated)
        (C T) FAIL,  // C1->T methylation transition
        (G .) PASS M5mC=None,  // G1 matches ref
        (C .) PASS M5mC=None,  // C2 matches ref
        (G .) PASS M5mC=1.,  // G2: genotype G/G, all reads show A (fully methylated)
        (G A) FAIL,  // G2->A methylation transition
        (C .) PASS M5mC=1.,  // C3: genotype C/C, all reads show T (fully methylated)
        (C T) FAIL,  // C3->T methylation transition
        (G .) PASS M5mC=None,  // G3 matches ref
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_strand_bias_in_methylation() -> Result<()> {
    // Test C->T only on OT strand, not OB
    // Assumption: All variants from one strand should still be reported correctly
    // In TAPS, OT (original top) strand shows C->T for unmethylated C
    // OB (original bottom) reads matching C suggest the position is methylated on forward strand
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ C A ] OB,
        [ C A ] OB,
        [ C A ] OB,
    );
    let records = test_call(segment, pileups, RecordFilters::all())?;
    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=1.,  // methylation evidence
        (C T) FAIL,  // strand-biased methylation transition
        (G .) PASS M5mC=1.,  // G matches on both strands
        (G A) FAIL,  // G matches on both strands
    ];
    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    // And now the other way around!
    // Test C->T only on OB strand -- i.e., not methylation evidence
    // (and same for G->A)
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OB,
        [ T G ] OB,
        [ T G ] OB,
        [ C A ] OT,
        [ C A ] OT,
        [ C A ] OT,
    );
    let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
    set_pass(&mut records[0], T);
    set_pass(&mut records[1], A);
    let records = reprocess(records)?;
    let expected_vcf = vcf_assert![
        (C T) PASS M5mC=None,
        (G A) PASS M5mC=None,
    ];
    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_cpg_both_positions_written() -> Result<()> {
    // - C has variant C->A (passes)
    // - G has no alts (matches reference)
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ A G ] OT,
        [ A G ] OT,
        [ A G ] OB,
        [ A G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
    set_pass(&mut records[0], A);
    let records = reprocess(records)?;

    // Always write both CpG positions
    let expected_vcf = vcf_assert![
        (C A) PASS,
        (G .) PASS
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::cpgs())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_other_pos_in_cpg_passes_flag() -> Result<()> {
    // Test the other_pos_in_cpg_passes flag behavior
    // When one position in a CpG passes, the other position's ref call should also pass
    // even if it would normally fail filters
    //
    // Assumption: C position fails filters, but G passes
    // Expected: C position passes (ref->. row) due to other_pos_in_cpg_passes flag
    // However, individual alts (C->T) still fail based on their own filters
    // The position-level flag makes the position pass, but not necessarily individual alts
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T A ] OT,
        [ T A ] OT,
        [ C G ] OB,
        [ C G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
    set_fail(&mut records[0], T);
    set_pass(&mut records[1], A);
    let records = reprocess(records)?;

    // After `propagate_cpg_pass_flags` runs, C position gets
    // other_pos_in_cpg_passes=true because G passes. This makes the position
    // pass (so ref->. is PASS), but C->T alt still fails on its own merits

    let expected_vcf = vcf_assert![
        (C .) PASS,  // Position passes due to other_pos_in_cpg_passes
        (C T) FAIL,  // Alt still fails based on its own filters (low ML, etc.)
        (G A) PASS,  // G->A passes normally
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// A het C/A at the C of a CpG: the `A` reads on either strand carry no CpG,
/// so β at both sides is over the CpG chromosome's reads alone, the `T`s on
/// OT and the `A`s on OB after a `C`: all of them methylated.
#[test]
fn a_non_conversion_het_leaves_beta_to_the_cpg_reads() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ A G ] OT,
        [ A G ] OT,
        [ A G ] OB,
        [ A G ] OB,
        [ C A ] OB,
        [ C A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
    set_fail(&mut records[0], T); // methylation evidence
    set_pass(&mut records[0], A); // real variant
    set_fail(&mut records[1], A); // methylation evidence
    let records = reprocess(records)?;

    let vcf_records = metrics_to_vcf(&records, RecordFilters::cpgs())?;
    let expected_vcf = vcf_assert![
        (C A) PASS M5mC=1.0 GT="0/1",
        (G .) PASS M5mC=1.0,
    ];
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn uncovered_reference_cpg_reported_with_all_and_cpgs_only() -> Result<()> {
    // A reference CpG whose C column has zero surviving reads (e.g. all reads
    // filtered out) must still be reported when `--all` and `--cpgs-only` are
    // combined, matching the contract documented on `RecordFilters`. Only
    // `--cpgs-only` on its own keeps requiring coverage.
    let (segment, mut pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ C G ] OT,
        [ C G ] OB,
    );
    pileups[0].reads = SimpleReads(Arc::from(Vec::<SimpleRead>::new()));

    let records = test_call(segment, pileups, RecordFilters::all())?;
    assert_eq!(records.len(), 2, "both reference CpG positions should survive pre-filtering");
    let c_record = &records[0];
    let g_record = &records[1];
    assert!(
        c_record.reference_cpg().is_some() && !c_record.covered(),
        "C column should be an uncovered CpG"
    );
    assert!(
        g_record.reference_cpg().is_some() && g_record.covered(),
        "G column should be a covered CpG"
    );

    let all_cpgs = RecordFilters { vcf_all: true, cpgs_only: true };
    let cpgs_only = RecordFilters { vcf_all: false, cpgs_only: true };

    assert!(all_cpgs.matches(c_record), "--all --cpgs-only must report the uncovered CpG");
    assert!(!cpgs_only.matches(c_record), "--cpgs-only alone must keep dropping uncovered CpGs");
    assert!(all_cpgs.matches(g_record), "--all --cpgs-only must report the covered CpG");
    assert!(cpgs_only.matches(g_record), "--cpgs-only must report the covered CpG");

    // The uncovered CpG should also make it into the VCF as a reference-only record.
    let vcf_records = metrics_to_vcf(&records, all_cpgs)?;
    let expected_vcf = vcf_assert![
        (C .),
        (G .),
    ];
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// The model passes both `T` and `A` at a CpG `C`, but the `T` reads are all
/// on the converted strand and the `C` reads on the other: the engine's
/// genotype is `C/A`, holds no copy of `T`, so `T` is not called and its
/// reads count as methylation, as if it had failed.
#[test]
fn a_passed_conversion_the_genotype_holds_no_copy_of_is_methylation() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ C G ] OB,
        [ C G ] OB,
        [ C G ] OB,
        [ A G ] OT,
        [ A G ] OB,
        [ A G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    set_pass(&mut records[0], A);
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C A) PASS M5mC=1. GT="0/1",
        (C T) FAIL,
        (G .) PASS,
    ];
    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    let filters: Vec<_> = vcf_records.iter().map(|r| r.filters.join(";")).collect();
    assert!(filters.iter().any(|f| f.contains("no_genotype_copy")), "{filters:?}");
    expected_vcf.matches(vcf_records)?;

    Ok(())
}
