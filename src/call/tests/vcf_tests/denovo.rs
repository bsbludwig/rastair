use crate::{call::tests::utils::*, pileups, vcf_assert};
use rastair_types::Base::*;

#[test]
fn test_denovo_cpg() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ G G ] Ref,
        [ C G ] OT,
        [ C G ] OT,
        [ C G ] OB,
    );

    let expected_vcf = vcf_assert![
        (G C) PASS M5mC=0.,  // de-novo CpG candidate
        (G .) PASS M5mC=0.,  // second position
    ];

    let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
    set_pass(&mut records[0], C); // set C alt to pass, creating de-novo CpG for sure
    let records = reprocess(records)?; // to propagate de-novo CpG flags

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_denovo_cpg_methylated() -> Result<()> {
    // First position: G with alts C and T
    // Second position: G with alt A (methylation transition G->A)
    let (segment, pileups) = pileups!(
        [ G G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ C A ] OB,
        [ C A ] OB,
        [ G A ] OB,
    );

    let expected_vcf = vcf_assert![
        (G C) PASS M5mC=1.,  // All alts at first position combined: T fails, C passes (any passing = PASS)
        (G .) PASS M5mC=1.,  // other de-novo CpG position
    ];

    let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
    set_pass(&mut records[0], C); // set C alt to pass, creating de-novo CpG for sure
    set_fail(&mut records[0], T); // set T alt to fail, it's a methylation evidence
    set_fail(&mut records[1], A); // set C alt to pass, creating de-novo CpG for sure
    let records = reprocess(records)?; // to recalculate genotypes and propagate de-novo CpG flags

    let vcf_records = metrics_to_vcf(&records)?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

// /// Edge case: CCC reference with heterozygous C>G at positions 1 and 2
// ///
// /// Position 1: C>G creates de-novo CpG at (1,2)
// /// Position 2: C>G creates de-novo CpG at (2,3)
// /// Position 2 has dual role:
// ///   - It's the "matching C" for position 1's de-novo CpG (`ThisIsTheMatchingC`)
// ///   - It also has C>G creating its own de-novo CpG (`ThisBecomesG`)
// #[test]
// fn test_adjacent_denovo_cpgs_dual_role_middle_position() -> Result<()> {
//     let (segment, pileups) = pileups!(
//         [ C C G ] Ref,
//         [ G G G ] OT,
//         [ G G G ] OT,
//         [ G G G ] OB,
//         [ C T G ] OT,
//         [ C C A ] OB,
//     );

//     let expected_vcf = vcf_assert![
//         (C .) PASS M5mC=0., // Position 1: C with de-novo CpG G
//         (C G) PASS M5mC=None,
//         (C .) PASS M5mC=1., // Weird but maybe correct: one T that fails ML makes C methylated
//         (C G) PASS M5mC=1., // FIXME: This is the methylation evidence of the G but it is wrong, should be 0.0
//         (C T) FAIL M5mC=None,
//         (G .) PASS M5mC=0.5,
//         (G A) FAIL M5mC=None,
//     ];

//     let records = test_call(segment, pileups, RecordFilters::cpgs())?;
//     let vcf_records = metrics_to_vcf(&records)?;
//     expected_vcf.matches(vcf_records)?;

//     Ok(())
// }

// #[test]
// fn test_adjacent_denovo_cpgs_dual_role_middle_position2() -> Result<()> {
//     let (segment, pileups) = pileups!(
//         [ C C C ] Ref,
//         [ G G G ] OT,
//         [ G G G ] OT,
//         [ G C G ] OB,
//         [ C T G ] OT,
//         [ C C A ] OB,
//     );

//     let expected_vcf = vcf_assert![
//         (C .) FAIL, // FIXME: Should pass as C in de-novo CpG
//         (C G) PASS M5mC=None,
//         (C .) PASS M5mC=1., // Weird but maybe correct: one T that fails ML makes C methylated
//         (C G) PASS M5mC=1., // FIXME: This is the methylation evidence of the G but it is wrong, should be 0.0
//         (C T) PASS, // FIXME: This should be FAIL, only 1 T
//         (C G) PASS M5mC=0.5,
//         (C A) FAIL M5mC=None,
//     ];

//     let records = test_call(segment, pileups, RecordFilters::cpgs())?;
//     let vcf_records = metrics_to_vcf(&records)?;
//     expected_vcf.matches(vcf_records)?;

//     Ok(())
// }

// #[test]
// fn test_a_and_t_cant_form_denovo() -> Result<()> {
//     // Test that A and T reference bases don't produce methylation evidence rows
//     // Currently, de-novo CpG creation is only considered when one of the C or G is in the reference
//     let (segment, pileups) = pileups!(
//         [ A T ] Ref,
//         [ C G ] OT,
//         [ C G ] OT,
//         [ C G ] OB,
//         [ C G ] OB,
//     );

//     let records = test_call(segment, pileups, RecordFilters::all())?;

//     let expected_vcf = vcf_assert![
//         (A C) PASS M5mC=None,
//         (T G) PASS M5mC=None,
//     ];

//     let vcf_records = metrics_to_vcf(&records)?;
//     expected_vcf.matches(vcf_records)?;

//     Ok(())
// }
