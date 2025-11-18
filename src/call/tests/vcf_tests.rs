//! Tests for VCF output from rastair call
use super::utils::*;
use crate::{pileups, vcf};
use rastair_types::{Base::*, Probability};

#[test]
fn test_simple_variant() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ A C G T ] Ref,
        [ A T G T ] OB,
        [ A T G T ] OB,
        [ A C G T ] OT,
        [ A C G T ] OT,
    );

    let expected_vcf = vcf![
        (C T),
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
        (C T),
        (G .),
    ];

    let records = test_call(segment, pileups, RecordFilters::cpgs())?;
    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_filter_status_matching() -> Result<()> {
    // Test that we can match PASS/FAIL status on records
    let (segment, pileups) = pileups!(
        [ C G T ] Ref,
        [ T G T ] OT,
        [ T G T ] OT,
        [ T A T ] OB,
        [ T A T ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;

    // Manipulate ML scores to make some records fail
    records[0].alt_filters_mut(T).unwrap().ml = Some(Probability::ZERO);
    records[1].alt_filters_mut(A).unwrap().ml = Some(Probability::ZERO);

    let expected_vcf = vcf![
        (C T) FAIL,
        (G A) FAIL,
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    // Now again but passing
    records[0].alt_filters_mut(T).unwrap().ml = Some(Probability::ONE);
    records[1].alt_filters_mut(A).unwrap().ml = Some(Probability::ONE);

    let expected_vcf = vcf![
        (C T) PASS,
        (G A) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_with_vcf_matcher() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ A C ] Ref,
        [ A T ] OB,
        [ G T ] OB,
        [ A C ] OT,
    );

    let expected_vcf = vcf![
        (A G),
        (C T),
    ];

    let records = test_call(segment, pileups, RecordFilters::all())?;
    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}
