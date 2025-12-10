use crate::{call::tests::utils::*, pileups, vcf::lowDp, vcf_assert};
use rastair_types::Base::*;

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
    let vcf_records = metrics_to_vcf(&records)?;
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

    // Reprocess to recalculate genotypes with the modified ML scores
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .,A) PASS M5mC=1.,  // methylation evidence, passing A
        (C T) FAIL,  // Other: T fails as real variant
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records)?;
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
    set_fail(&mut records[0], A); // Low ML, but C->A is not methylation transition
    set_fail(&mut records[1], C); // Low ML, but G->C is not methylation transition

    let expected_vcf = vcf_assert![
        (C .) FAIL,  // No ref->. row because C->A is not methylation evidence
        (C A) PASS,  // FIXME: Skip because C->A is not methylation evidence
        (G .) FAIL,  // FIXME: Skip because G->C is not methylation evidence
        (G C) PASS,  // No ref->. row because G->C is not methylation evidence
    ];

    let vcf_records = metrics_to_vcf(&records)?;
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

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=1.,  // methylation evidence (only one ref->. row)
        (C T) FAIL,  // the methylation transition
        (G .) PASS M5mC=0.,
    ];

    let vcf_records = metrics_to_vcf(&records)?;
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

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=0.,
        (G .) PASS M5mC=1.,  // methylation evidence
        (G A) FAIL,  // low confidence methylation transition
    ];

    let vcf_records = metrics_to_vcf(&records)?;
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

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_all_methylation_transitions_failing() -> Result<()> {
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
    set_fail(&mut records[0], T);
    set_fail(&mut records[1], A);

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=1.,  // C is methylated
        (C T) FAIL,  // methylation transition
        (G .) PASS M5mC=1.,  // G is methylated
        (G A) FAIL,  // methylation transition
    ];

    let vcf_records = metrics_to_vcf(&records)?;
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

    let expected_vcf = vcf_assert![
        (C .) FAIL M5mC=None,
        (C T) PASS,
        (G .) FAIL M5mC=None,
        (G A) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

// #[test]
// fn test_methylation_transition_with_other_filters_failing() -> Result<()> {
//     // Test C->T with high ML but other filters fail (e.g., low depth, low quality)
//     // Current implementation trusts ML score if it is set, so should output PASS
//     let (segment, pileups) = pileups!(
//         [ C G ] Ref,
//         [ T G ] OT,
//         [ T G ] OT,
//         [ T G ] OB,
//         [ T G ] OB,
//     );

//     let mut records = test_call(segment, pileups, RecordFilters::all())?;
//     set_pass(&mut records[0], T); // High ML score
//     records[0].alt_filters_mut(T).unwrap().filters.add(lowDp, || true);

//     let expected_vcf = vcf_assert![
//         (C .) PASS,
//         (C T) PASS,  // PASS despite other filters failing
//         (G .) PASS,
//     ];

//     let vcf_records = metrics_to_vcf(&records)?;
//     expected_vcf.matches(vcf_records)?;

//     // Now, without ML
//     records[0].alt_filters_mut(T).unwrap().ml = None;

//     let expected_vcf = vcf_assert![
//         (C .) PASS,
//         (C T) FAIL,  // FAIL because ML is None and other filters fail
//         (G .) PASS,
//     ];

//     let vcf_records = metrics_to_vcf(&records)?;
//     expected_vcf.matches(vcf_records)?;

//     Ok(())
// }

// #[test]
// fn test_c_with_t_and_a_both_low_ml() -> Result<()> {
//     // Test C with both T and A alts, both failing ML threshold
//     // Assumption: C->T is methylation transition, C->A is not
//     // Expected: ref->. row for methylation evidence, then both alts fail separately
//     let (segment, pileups) = pileups!(
//         [ C G ] Ref,
//         [ T G ] OT,
//         [ T G ] OT,
//         [ A G ] OB,
//         [ A G ] OB,
//     );

//     let mut records = test_call(segment, pileups, RecordFilters::all())?;
//     set_fail(&mut records[0], T); // Low ML - methylation evidence
//     set_fail(&mut records[0], A); // Low ML - but not methylation transition

//     // Reprocess to recalculate genotypes with the modified ML scores
//     let records = reprocess(records)?;

//     let expected_vcf = vcf_assert![
//         (C .) PASS M5mC=1.,  // methylation evidence from C->T
//         (C T) FAIL,
//         (C A) FAIL,
//         (G .) PASS M5mC=0.,
//     ];

//     let vcf_records = metrics_to_vcf(&records)?;
//     expected_vcf.matches(vcf_records)?;

//     Ok(())
// }

// #[test]
// fn test_g_with_a_and_c_alts() -> Result<()> {
//     // Test G->A methylation evidence + G->C in a GC context (reverse CpG)
//     // Context: G followed by C (GC = CpG on reverse strand), where G has alts A and C
//     // Assumption: Multiple biological signals at one position should all be reported
//     let (segment, pileups) = pileups!(
//         [ G C ] Ref,
//         [ A C ] OT,
//         [ A C ] OT,
//         [ C C ] OB,
//         [ C C ] OB,
//     );

//     let mut records = test_call(segment, pileups, RecordFilters::all())?;
//     set_fail(&mut records[0], A); // G->A low ML - methylation evidence
//     set_pass(&mut records[0], C); // G->C high ML - non-methylation variant

//     // Reprocess to recalculate genotypes with the modified ML scores
//     let records = reprocess(records)?;

//     let expected_vcf = vcf_assert![
//         (G C) PASS,  // All alts combined: A fails, C passes (any passing = PASS)
//         (G A) FAIL,  // All alts combined: A fails, C passes (any passing = PASS)
//     ];

//     let vcf_records = metrics_to_vcf(&records)?;
//     expected_vcf.matches(vcf_records)?;

//     Ok(())
// }

// #[test]
// fn test_cpg_islands_multiple_positions() -> Result<()> {
//     // Test CGCGCG reference with various methylation patterns
//     // Assumption: Multiple consecutive CpG sites should each be handled independently
//     // Pattern: C1 methylated, G1 not, C2 not, G2 methylated, C3 methylated, G3 not
//     let (segment, pileups) = pileups!(
//         [ C G C G C G ] Ref,
//         [ T G C A T G ] OT,
//         [ T G C A T G ] OT,
//         [ T G C A T G ] OB,
//         [ T G C A T G ] OB,
//     );

//     let mut records = test_call(segment, pileups, RecordFilters::all())?;
//     // C1->T: low ML (methylation)
//     set_fail(&mut records[0], T);
//     // G1: matches ref (no alt)
//     // C2: matches ref (no alt)
//     // G2->A: low ML (methylation)
//     set_fail(&mut records[3], A);
//     // C3->T: low ML (methylation)
//     set_fail(&mut records[4], T);
//     // G3: matches ref (no alt)

//     let expected_vcf = vcf_assert![
//         (C .) PASS M5mC=1.,  // C1 methylated
//         (C T) FAIL,  // C1->T methylation transition
//         (G .) PASS M5mC=0.,  // G1 matches ref
//         (C .) PASS M5mC=0.,  // C2 matches ref
//         (G .) PASS M5mC=1.,  // G2 methylated
//         (G A) FAIL,  // G2->A methylation transition
//         (C .) PASS M5mC=1.,  // C3 methylated
//         (C T) FAIL,  // C3->T methylation transition
//         (G .) PASS M5mC=0.,  // G3 matches ref
//     ];

//     let vcf_records = metrics_to_vcf(&records)?;
//     expected_vcf.matches(vcf_records)?;

//     Ok(())
// }

// #[test]
// fn test_strand_bias_in_methylation() -> Result<()> {
//     // Test C->T only on OT strand, not OB
//     // Assumption: All variants from one strand should still be reported correctly
//     // In TAPS, OT (original top) strand shows C->T for unmethylated C
//     // OB (original bottom) reads matching C suggest the position is methylated on forward strand
//     let (segment, pileups) = pileups!(
//         [ C G ] Ref,
//         [ T G ] OT,
//         [ T G ] OT,
//         [ T G ] OT,
//         [ C A ] OB,
//         [ C A ] OB,
//         [ C A ] OB,
//     );
//     let records = test_call(segment, pileups, RecordFilters::all())?;
//     let expected_vcf = vcf_assert![
//         (C .) PASS M5mC=1.,  // methylation evidence
//         (C T) FAIL,  // strand-biased methylation transition
//         (G .) PASS M5mC=1.,  // G matches on both strands
//         (G A) FAIL,  // G matches on both strands
//     ];
//     let vcf_records = metrics_to_vcf(&records)?;
//     expected_vcf.matches(vcf_records)?;

//     // And now the other way around!
//     // Test C->T only on OB strand -- i.e., not methylation evidence
//     // (and same for G->A)
//     let (segment, pileups) = pileups!(
//         [ C G ] Ref,
//         [ T G ] OB,
//         [ T G ] OB,
//         [ T G ] OB,
//         [ C A ] OT,
//         [ C A ] OT,
//         [ C A ] OT,
//     );
//     let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
//     set_pass(&mut records[0], T);
//     set_pass(&mut records[1], A);
//     let expected_vcf = vcf_assert![
//         (C .,T) PASS M5mC=0.,
//         (G .,A) PASS M5mC=0.,
//     ];
//     let vcf_records = metrics_to_vcf(&records)?;
//     expected_vcf.matches(vcf_records)?;

//     Ok(())
// }

// #[test]
// fn test_cpg_both_positions_written() -> Result<()> {
//     // - C has variant C->A (passes)
//     // - G has no alts (matches reference)
//     let (segment, pileups) = pileups!(
//         [ C G ] Ref,
//         [ A G ] OT,
//         [ A G ] OT,
//         [ A G ] OB,
//         [ A G ] OB,
//     );

//     let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
//     set_pass(&mut records[0], A);

//     // Always write both CpG positions
//     let expected_vcf = vcf_assert![
//         (C .,A) PASS,  // C with potential methylation evidence
//         (G .) PASS   // Should write G to maintain CpG context
//     ];

//     let vcf_records = metrics_to_vcf(&records)?;
//     expected_vcf.matches(vcf_records)?;

//     Ok(())
// }

// #[test]
// fn test_other_pos_in_cpg_passes_flag() -> Result<()> {
//     // Test the other_pos_in_cpg_passes flag behavior
//     // When one position in a CpG passes, the other position's ref call should also pass
//     // even if it would normally fail filters
//     //
//     // Assumption: C position fails filters, but G passes
//     // Expected: C position passes (ref->. row) due to other_pos_in_cpg_passes flag
//     // However, individual alts (C->T) still fail based on their own filters
//     // The position-level flag makes the position pass, but not necessarily individual alts
//     let (segment, pileups) = pileups!(
//         [ C G ] Ref,
//         [ T A ] OT,
//         [ T A ] OT,
//         [ C G ] OB,
//         [ C G ] OB,
//     );

//     let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;

//     // Make C->T fail
//     set_fail(&mut records[0], T);

//     // Make G->A pass
//     set_pass(&mut records[1], A);

//     // After `propagate_cpg_pass_flags` runs, C position gets
//     // other_pos_in_cpg_passes=true because G passes. This makes the position
//     // pass (so ref->. is PASS), but C->T alt still fails on its own merits

//     let expected_vcf = vcf_assert![
//         (C .) PASS,  // Position passes due to other_pos_in_cpg_passes
//         (C T) FAIL,  // Alt still fails based on its own filters (low ML, etc.)
//         (G .,A) PASS,  // G->A passes normally
//     ];

//     let vcf_records = metrics_to_vcf(&records)?;
//     expected_vcf.matches(vcf_records)?;

//     Ok(())
// }
