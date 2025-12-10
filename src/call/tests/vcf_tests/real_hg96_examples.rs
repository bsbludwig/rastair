use crate::{call::tests::utils::*, pileups, vcf_assert};
use rastair_types::Base::*;

#[test]
fn test_denovo_cpg_that_is_variant_hg96_chr20_75254() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C C ] Ref,
        [ C C ] OT,
        [ C C ] OB,
        [ A C ] OT,
        [ A C ] OB,
        [ A C ] OB,
        [ A G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::variants())?;
    set_fail(&mut records[0], A);

    // Reprocess to recalculate genotypes with the modified ML scores
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C A) FAIL,  // First position with A alt
        (C G) FAIL,  // Second position with G alt
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_denovo_cpg_hg96_chr20_76962() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ T G ] Ref,
        [ C G ] OT,
        [ C G ] OT,
        [ C G ] OT,
        [ C G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
    set_pass(&mut records[0], C);

    let expected_vcf = vcf_assert![
        (T C) PASS,  // De-novo CpG creation
        (G .) PASS,  // Write the G from the de-novo CpG
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_cpg_variant_hg96_chr20_65899() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ C A ] OT,
        [ C A ] OB,
        [ C A ] OB,
        [ C A ] OB,
    );

    let records = test_call(segment, pileups, RecordFilters::all())?;

    let expected_vcf = vcf_assert![
        (C .) PASS,
        (G A) PASS, // Actual variant
    ];

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}
