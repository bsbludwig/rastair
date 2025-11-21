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
#[ignore = "TODO: Future feature combining multiple alts into one row"]
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
