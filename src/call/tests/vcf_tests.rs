//! Tests for VCF output from rastair call
use super::utils::*;
use crate::{pileups, vcf};
use rastair_types::Base::*;

#[test]
fn test_simple_variant() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ A C G T ] Ref,
        [ A A G T ] OB,
        [ A A G T ] OB,
        [ A C G T ] OT,
        [ A C G T ] OT,
    );

    let expected_vcf = vcf![
        (C A),
        (G .),
    ];

    let records = test_call(segment, pileups, RecordFilters::variants())?;
    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_all_matching_ref() -> Result<()> {
    // All reads match reference
    let (segment, pileups) = pileups!(
        [ A C G T ] Ref,
        [ A C G T ] OT,
        [ A C G T ] OT,
        [ A C G T ] OB,
    );

    // Test with variants filter: There are no variants
    let records = test_call(segment.clone(), pileups.clone(), RecordFilters::variants())?;
    // FIXME: Why are we printing CpGs anyway?
    let expected_vcf = vcf![
        (C .),
        (G .),
    ];
    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    // Now test with CpG filter: Should match the middle two positions
    let records = test_call(segment, pileups, RecordFilters::cpgs())?;
    let expected_vcf = vcf![
        (C .),
        (G .),
    ];
    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
#[ignore = "TODO: combining multiple alts into one row"]
fn test_multiple_alt_alleles() -> Result<()> {
    // Multiple different alt alleles at the same position
    let (segment, pileups) = pileups!(
        [ C ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ A ] OB,
        [ G ] OB,
    );
    let expected_vcf = vcf![
        (C T,A,G),
    ];

    let records = test_call(segment, pileups, RecordFilters::variants())?;
    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_cpg_context() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ C G ] OB,
        [ C G ] OB,
    );

    let expected_vcf = vcf![
        (C .) PASS,
        (C T) FAIL,
        (G .) PASS,
    ];

    let records = test_call(segment, pileups, RecordFilters::cpgs())?;
    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_filter_status_matching() -> Result<()> {
    // Test that we can match PASS/FAIL status on records
    // C->T and G->A are methylation transitions, so when they fail ML threshold,
    // we output both ref->. (PASS) and ref->alt (FAIL) rows
    let (segment, pileups) = pileups!(
        [ C G T ] Ref,
        [ T G T ] OT,
        [ T G T ] OT,
        [ T A T ] OB,
        [ T A T ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;

    // Manipulate ML scores to make some records fail
    // When C->T and G->A fail, they become methylation evidence
    set_fail(&mut records[0], T);
    set_fail(&mut records[1], A);

    let expected_vcf = vcf![
        (C .) PASS,  // methylation evidence: C is methylated
        (C T) FAIL,  // low confidence variant
        (G .) PASS,  // methylation evidence: G is methylated
        (G A) FAIL,  // low confidence variant
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    // Now again but passing - no methylation evidence, just variants
    set_pass(&mut records[0], T);
    set_pass(&mut records[1], A);

    let expected_vcf = vcf![
        (C T) PASS,
        (G A) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_denovo_cpg() -> Result<()> {
    // First position: G with alts C and T
    // Second position: G with alt A (methylation transition G->A)
    let (segment, pileups) = pileups!(
        [ G G ] Ref,
        [ C G ] OT,
        [ C G ] OT,
        [ T A ] OB,
        [ T A ] OB,
    );

    let expected_vcf = vcf![
        (G C) PASS,  // de-novo CpG candidate
        (G T) FAIL,  // other alt at first position
        (G .) PASS,  // second position: no alts match ref, or methylation evidence
        (G A) FAIL,  // second position: G->A methylation transition with low ML
    ];

    let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
    set_pass(&mut records[0], C); // set C alt to pass, creating de-novo CpG for sure

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_c_to_t_high_ml_score() -> Result<()> {
    // Test C->T transition with HIGH ML score (true variant, not methylation)
    // Assumption: C->T with high ML score is a real variant, not methylation evidence
    let (segment, pileups) = pileups!(
        [ C ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ T ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::variants())?;
    set_pass(&mut records[0], T); // High ML score - real variant

    let expected_vcf = vcf![
        (C T) PASS,  // Only the variant row, no ref->. row
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_g_to_a_high_ml_score() -> Result<()> {
    // Test G->A transition with HIGH ML score (true variant, not methylation)
    // Assumption: G->A with high ML score is a real variant, not methylation evidence
    let (segment, pileups) = pileups!(
        [ G ] Ref,
        [ A ] OT,
        [ A ] OT,
        [ A ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::variants())?;
    set_pass(&mut records[0], A); // High ML score - real variant

    let expected_vcf = vcf![
        (G A) PASS,  // Only the variant row, no ref->. row
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_mixed_methylation_and_real_variants() -> Result<()> {
    // Test C with both methylation evidence (C->T low ML) and real variant (C->A high ML)
    // Assumption: When we have both types, we output ref->. plus both alt rows
    let (segment, pileups) = pileups!(
        [ C ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ A ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::variants())?;
    set_fail(&mut records[0], T); // Low ML - methylation evidence
    set_pass(&mut records[0], A); // High ML - real variant

    let expected_vcf = vcf![
        (C .) PASS,  // methylation evidence
        (C T) FAIL,  // low confidence methylation transition
        (C A) PASS,  // real variant
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

    let mut records = test_call(segment, pileups, RecordFilters::variants())?;
    set_fail(&mut records[0], A); // Low ML, but C->A is not methylation transition
    set_fail(&mut records[1], C); // Low ML, but G->C is not methylation transition

    let expected_vcf = vcf![
        (C A) FAIL,  // No ref->. row because C->A is not methylation evidence
        (G C) FAIL,  // No ref->. row because G->C is not methylation evidence
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_multiple_methylation_transitions_same_position() -> Result<()> {
    // Edge case: Both C->T with low ML (shouldn't happen in practice, but test it)
    // Assumption: Multiple failing methylation transitions still produce one ref->. row
    let (segment, pileups) = pileups!(
        [ C ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ T ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::variants())?;
    set_fail(&mut records[0], T); // Low ML - methylation evidence

    let expected_vcf = vcf![
        (C .) PASS,  // methylation evidence (only one ref->. row)
        (C T) FAIL,  // the methylation transition
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_a_and_t_bases_never_methylation_evidence() -> Result<()> {
    // Test that A and T reference bases never produce methylation evidence rows
    // Assumption: Only C and G can be methylated in TAPS, so A/T never get ref->. from methylation
    let (segment, pileups) = pileups!(
        [ A T ] Ref,
        [ C G ] OT,
        [ C G ] OT,
        [ C G ] OB,
        [ C G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::variants())?;
    set_fail(&mut records[0], C); // Low ML, but A->C is not methylation
    set_fail(&mut records[1], G); // Low ML, but T->G is not methylation

    let expected_vcf = vcf![
        (A C) FAIL,  // No ref->. because A cannot show methylation evidence
        (T G) FAIL,  // No ref->. because T cannot show methylation evidence
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_reverse_strand_methylation() -> Result<()> {
    // Test G->A methylation evidence on reverse strand
    // Assumption: G->A with low ML on reverse strand is methylation, just like C->T on forward
    let (segment, pileups) = pileups!(
        [ G ] Ref,
        [ A ] OT,
        [ A ] OT,
        [ G ] OB,
        [ G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::variants())?;
    set_fail(&mut records[0], A); // Low ML - methylation evidence on reverse strand

    let expected_vcf = vcf![
        (G .) PASS,  // methylation evidence
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

    let records = test_call(segment, pileups, RecordFilters::cpgs())?;

    let expected_vcf = vcf![
        (C .),  // Single row for matching reference
        (G .),  // Single row for matching reference
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

    let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
    set_fail(&mut records[0], T);
    set_fail(&mut records[1], A);

    let expected_vcf = vcf![
        (C .) PASS,  // C is methylated
        (C T) FAIL,  // methylation transition
        (G .) PASS,  // G is methylated
        (G A) FAIL,  // methylation transition
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
#[ignore = "TODO: what happens when ml is None"]
fn test_c_to_t_with_no_ml_score() -> Result<()> {
    // Test C->T transition when ML model wasn't run
    // Question: When ml is None but we have C->T, does it fall back to filter checks?
    // Should we still treat it as methylation evidence, or only when ML explicitly fails?
    let (segment, pileups) = pileups!(
        [ C ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ T ] OB,
    );

    let records = test_call(segment, pileups, RecordFilters::variants())?;

    // What should we expect here?
    let expected_vcf = vcf![
        // (C .) PASS,  // Should we output this?
        // (C T) ???,   // What filter status?
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
#[ignore = "TODO: C->T passes ML but fails depth/quality filters"]
fn test_methylation_transition_with_other_filters_failing() -> Result<()> {
    // Test C->T with high ML but other filters fail (e.g., low depth, low quality)
    // Question: Should we treat this as methylation evidence because it's C->T,
    // or as a regular failing variant because ML score is high?
    // Does methylation evidence logic only apply when ML specifically fails?
    let (segment, pileups) = pileups!(
        [ C ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ T ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::variants())?;
    set_pass(&mut records[0], T); // High ML score

    // Now manually add other filters (low depth, etc.) to the alt
    // How do we do this in tests? Need helper function?

    // What should we expect here?
    let expected_vcf = vcf![
        // (C .) PASS,  // Or no ref->. because ML passed?
        // (C T) FAIL,  // Fails due to other filters
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_c_with_t_and_a_both_low_ml() -> Result<()> {
    // Test C with both T and A alts, both failing ML threshold
    // Assumption: C->T is methylation transition, C->A is not
    // Expected: ref->. row for methylation evidence, then both alts fail separately
    let (segment, pileups) = pileups!(
        [ C ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ A ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::variants())?;
    set_fail(&mut records[0], T); // Low ML - methylation evidence
    set_fail(&mut records[0], A); // Low ML - but not methylation transition

    let expected_vcf = vcf![
        (C .) PASS,  // methylation evidence from C->T
        (C T) FAIL,  // methylation transition with low ML
        (C A) FAIL,  // non-methylation transition with low ML, no ref->. for this
    ];

    let vcf_records = metrics_to_vcf(&records)?;
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

    let expected_vcf = vcf![
        (G .) PASS,  // methylation evidence from G->A
        (G A) FAIL,  // methylation transition with low ML
        (G C) PASS,  // non-methylation variant
    ];

    let vcf_records = metrics_to_vcf(&records)?;
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

    let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
    // C1->T: low ML (methylation)
    set_fail(&mut records[0], T);
    // G1: matches ref (no alt)
    // C2: matches ref (no alt)
    // G2->A: low ML (methylation)
    set_fail(&mut records[3], A);
    // C3->T: low ML (methylation)
    set_fail(&mut records[4], T);
    // G3: matches ref (no alt)

    let expected_vcf = vcf![
        (C .) PASS,  // C1 methylated
        (C T) FAIL,  // C1->T methylation transition
        (G .) PASS,  // G1 matches ref
        (C .) PASS,  // C2 matches ref
        (G .) PASS,  // G2 methylated
        (G A) FAIL,  // G2->A methylation transition
        (C .) PASS,  // C3 methylated
        (C T) FAIL,  // C3->T methylation transition
        (G .) PASS,  // G3 matches ref
    ];

    let vcf_records = metrics_to_vcf(&records)?;
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
        [ C G ] OB,
        [ C G ] OB,
        [ C G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
    set_fail(&mut records[0], T); // C->T with low ML from OT strand only

    let expected_vcf = vcf![
        (C .) PASS,  // methylation evidence
        (C T) FAIL,  // strand-biased methylation transition
        (G .) PASS,  // G matches on both strands
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
#[ignore = "TODO: write both CpG positions"]
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

    // What should we expect?
    let expected_vcf = vcf![
        (C A) PASS,  // C has variant
        (G .) PASS   // Should write G to maintain CpG context?
    ];

    let vcf_records = metrics_to_vcf(&records)?;
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

    // Make C->T fail
    set_fail(&mut records[0], T);

    // Make G->A pass
    set_pass(&mut records[1], A);

    // After cpg_sites.rs propagate_cpg_pass_flags runs,
    // C position gets other_pos_in_cpg_passes=true because G passes
    // This makes the position pass (so ref->. is PASS), but C->T alt still fails on its own merits

    let expected_vcf = vcf![
        (C .) PASS,  // Position passes due to other_pos_in_cpg_passes
        (C T) FAIL,  // Alt still fails based on its own filters (low ML, etc.)
        (G A) PASS,  // G->A passes normally
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
#[ignore = "TODO"]
fn test_denovo_cpg_both_positions_with_methylation() -> Result<()> {
    // Test complex scenario combining multiple TODOs:
    // - A->C creates de-novo CpG (now have CG)
    // - The newly created C shows methylation evidence (C->T on some reads)
    // - Should write both positions of the de-novo CpG
    //
    // Assumption: When de-novo CpG is created and shows methylation,
    // both positions should be written to maintain biological context
    //
    // Example: Reference is AG
    // - Some reads show CG (A->C, creating de-novo CpG)
    // - Some reads show TG (suggesting unmethylated C)
    // Expected output: Both C and G positions with methylation info
    let (segment, pileups) = pileups!(
        [ A G ] Ref,
        [ C G ] OT,  // Creates de-novo CpG
        [ C G ] OT,
        [ T G ] OB,  // Methylation evidence (unmethylated C->T)
        [ T G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
    // A->C creates de-novo CpG and should pass
    set_pass(&mut records[0], C);
    // A->T is methylation transition with low ML
    set_fail(&mut records[0], T);

    // What should we expect?
    let expected_vcf = vcf![
        // (A .) ???,  // Should we show unmethylated A?
        // (A C) PASS,  // De-novo CpG creation
        // (A T) FAIL,  // Methylation evidence
        // (G .) ???,  // Should we write the G from the de-novo CpG?
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}
