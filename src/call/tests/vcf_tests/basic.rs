use crate::{call::tests::utils::*, pileups, vcf_assert};
use rastair_types::Base::*;

#[test]
fn test_simple_variant() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ A T ] Ref,
        [ A C ] OB,
        [ A C ] OB,
        [ A C ] OT,
        [ A C ] OT,
    );

    let expected_vcf = vcf_assert![
        (A .) PASS,
        (T C) PASS,
    ];

    let records = test_call(segment, pileups, RecordFilters::all())?;
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
    // We are printing CpGs anyway
    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=0.0,
        (G .) PASS M5mC=0.0,
    ];
    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    // Now test with CpG filter: Should match the middle two positions
    // so same as above
    let records = test_call(segment.clone(), pileups.clone(), RecordFilters::cpgs())?;
    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=0.0,
        (G .) PASS M5mC=0.0,
    ];
    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    // Now with all filter: prints all positions
    let records = test_call(segment, pileups, RecordFilters::all())?;
    let expected_vcf = vcf_assert![
        (A .) PASS M5mC=None,
        (C .) PASS M5mC=0.0,
        (G .) PASS M5mC=0.0,
        (T .) PASS M5mC=None,
    ];
    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_multiple_alt_alleles() -> Result<()> {
    // Multiple different alt alleles at the same position
    let (segment, pileups) = pileups!(
        [ C ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ A ] OB,
        [ G ] OB,
    );
    let expected_vcf = vcf_assert![
        (C T,G) PASS,
        (C A) FAIL,
    ];

    let mut records = test_call(segment, pileups, RecordFilters::variants())?;
    set_pass(&mut records[0], T);
    set_pass(&mut records[0], G);
    set_fail(&mut records[0], A);
    let records = reprocess(records)?;

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_c_to_t_outside_cpg() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C C ] Ref,
        [ T C ] OT,
        [ T C ] OT,
        [ T C ] OB,
        [ T C ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::variants())?;
    set_pass(&mut records[0], T);
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C T) PASS M5mC=None,  // Only the variant row, no ref->. row
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_g_to_a_outside_cpg() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ G G ] Ref,
        [ A G ] OT,
        [ A G ] OT,
        [ A G ] OB,
        [ A G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::variants())?;
    set_pass(&mut records[0], A);
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (G A) PASS,  // Only the variant row, no ref->. row
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_c_to_t_outside_cpg_with_no_ml() -> Result<()> {
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

    let records = test_call(segment, pileups, RecordFilters::all())?;

    let expected_vcf = vcf_assert![
        (C T) PASS M5mC=None
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}
