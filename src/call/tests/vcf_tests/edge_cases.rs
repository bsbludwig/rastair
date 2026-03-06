//! Edge case tests for methylation calling after refactoring
//!
//! These tests cover subtle corner cases that were either untested or worked
//! by accident before the methylation calculation refactoring.

use crate::{call::tests::utils::*, pileups, vcf_assert};
use seqair_types::Base::*;

/// Edge case 1: Het with non-confounding alt on original CpG
///
/// ref=C, genotype=Het(C/A). T reads are purely methylation, no excess
/// correction needed because the het alt (A) is not the confounding base (T).
///
/// This was the old TODO case in the pre-refactor code.
#[test]
fn het_with_non_confounding_alt_original_cpg() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ C G ] OT,  // C reads (unmodified)
        [ C G ] OT,
        [ C G ] OT,
        [ A G ] OT,  // A reads (real variant)
        [ A G ] OT,
        [ A G ] OT,
        [ T G ] OT,  // T reads (methylation evidence)
        [ T G ] OT,
        [ T G ] OT,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], A); // A is a real variant
    set_fail(&mut records[0], T); // T is methylation evidence
    let records = reprocess(records)?;

    // Genotype should be Het(C/A)
    // Beta should be 3 T / (3 C + 3 T) = 0.5
    // No excess correction because A is not the confounding base T
    let expected_vcf = vcf_assert![
        (C A) PASS M5mC=0.5 GT="0/1",
        (C T) FAIL,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// Edge case 1b: Same but with G-side of CpG
///
/// ref=G, genotype=Het(G/C). A reads are purely methylation, no excess
/// correction needed because the het alt (C) is not the confounding base (A).
#[test]
fn het_with_non_confounding_alt_original_cpg_g_side() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ C G ] OB,  // G reads (unmodified)
        [ C G ] OB,
        [ C G ] OB,
        [ C C ] OB,  // C reads (real variant)
        [ C C ] OB,
        [ C C ] OB,
        [ C A ] OB,  // A reads (methylation evidence)
        [ C A ] OB,
        [ C A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], C); // C is a real variant
    set_fail(&mut records[1], A); // A is methylation evidence
    let records = reprocess(records)?;

    // Genotype should be Het(G/C)
    // Beta should be 3 A / (3 G + 3 A) = 0.5
    // No excess correction because C is not the confounding base A
    let expected_vcf = vcf_assert![
        (C .) PASS,
        (G C) PASS M5mC=0.5 GT="0/1",
        (G A) FAIL,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// Edge case 2: HomAlt on original CpG
///
/// ref=C, genotype=HomAlt(T). No C allele → beta=0.0.
/// This is not tested in the unit tests, only implicitly in integration tests.
#[test]
fn homalt_on_original_cpg() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OB,
        [ T G ] OB,
        [ T G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T); // T is a real variant (passes ML)
    let records = reprocess(records)?;

    // Genotype should be HomAlt(T/T)
    // Beta should be 0.0 because there's no C allele to methylate
    let expected_vcf = vcf_assert![
        (C T) PASS M5mC=0.0 GT="1/1",
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// Edge case 2b: Same but with G-side of CpG
///
/// ref=G, genotype=HomAlt(A). No G allele → beta=0.0.
#[test]
fn homalt_on_original_cpg_g_side() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ C A ] OT,
        [ C A ] OT,
        [ C A ] OT,
        [ C A ] OB,
        [ C A ] OB,
        [ C A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], A); // A is a real variant (passes ML)
    let records = reprocess(records)?;

    // Genotype should be HomAlt(A/A)
    // Beta should be 0.0 because there's no G allele to methylate
    let expected_vcf = vcf_assert![
        (C .) PASS,
        (G A) PASS M5mC=0.0 GT="1/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// Edge case 3: HomAlt on denovo CpG
///
/// ref=T, alt=G denovo, genotype=HomAlt(G). CpG on both chromosomes → normal beta.
/// This differs from edge case 2 because it's a denovo CpG, not an original one.
/// The T→G creates a CpG with the previous C, so this is a G-side denovo CpG.
#[test]
fn homalt_on_denovo_cpg() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C T ] Ref,
        [ C G ] OB,  // All G reads (denovo CpG on G-side)
        [ C G ] OB,
        [ C G ] OB,
        [ C A ] OB,  // A reads (methylation on G-side)
        [ C A ] OB,
        [ C A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], G); // G is a real variant creating denovo CpG
    set_fail(&mut records[1], A); // A is methylation evidence
    let records = reprocess(records)?;

    // Genotype should be HomAlt(G/G) at position 1
    // Beta should be 3 A / (3 A + 3 G) = 0.5 (NOT 0.0 like original CpG)
    // Because the CpG exists on both chromosomes via the G alt
    let expected_vcf = vcf_assert![
        (C .) PASS,
        (T G) PASS M5mC=0.5 GT="1/1",
        (T A) FAIL,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// Edge case 3b: Same but with C-side denovo CpG
///
/// ref=A, alt=C denovo, genotype=HomAlt(C). CpG on both chromosomes → normal beta.
/// The A→C creates a CpG with the next G, so this is a C-side denovo CpG.
#[test]
fn homalt_on_denovo_cpg_c_side() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ A G ] Ref,
        [ C G ] OT,  // All C reads (denovo CpG on C-side)
        [ C G ] OT,
        [ C G ] OT,
        [ T G ] OT,  // T reads (methylation on C-side)
        [ T G ] OT,
        [ T G ] OT,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C); // C is a real variant creating denovo CpG
    set_fail(&mut records[0], T); // T is methylation evidence
    let records = reprocess(records)?;

    // Genotype should be HomAlt(C/C) at position 0
    // Beta should be 3 T / (3 T + 3 C) = 0.5 (NOT 0.0 like original CpG)
    // Because the CpG exists on both chromosomes via the C alt
    let expected_vcf = vcf_assert![
        (A C) PASS M5mC=0.5 GT="1/1",
        (A T) FAIL,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// Edge case 4: Denovo with unmod reads but no confounding reads
///
/// ref=A, alt=G denovo, only G reads on OB, no A reads.
/// Should be beta=0.0, not NoEvidence (the CpG exists but is unmethylated).
#[test]
fn denovo_with_unmod_reads_no_confounding() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C A ] Ref,
        [ C G ] OB,  // G reads only (unmodified)
        [ C G ] OB,
        [ C G ] OB,
        [ C G ] OB,
        [ C G ] OB,
        // No A reads (confounding base for G-side)
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], G); // G is a real variant creating denovo CpG
    let records = reprocess(records)?;

    // Beta should be 0.0 (not NoEvidence) because:
    // - The CpG exists (we have G reads)
    // - There's no methylation evidence (no A reads)
    // - For denovo CpGs, having only unmodified reads is meaningful
    let expected_vcf = vcf_assert![
        (C .) PASS,
        (A G) PASS M5mC=0.0,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// Edge case 4b: Same but with C-side denovo CpG
///
/// ref=T, alt=C denovo, only C reads on OT, no T reads.
/// Should be beta=0.0, not NoEvidence.
#[test]
fn denovo_with_unmod_reads_no_confounding_c_side() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ T G ] Ref,
        [ C G ] OT,  // C reads only (unmodified)
        [ C G ] OT,
        [ C G ] OT,
        [ C G ] OT,
        [ C G ] OT,
        // No T reads (confounding base for C-side)
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C); // C is a real variant creating denovo CpG
    let records = reprocess(records)?;

    // Beta should be 0.0 (not NoEvidence) because:
    // - The CpG exists (we have C reads)
    // - There's no methylation evidence (no T reads)
    // - For denovo CpGs, having only unmodified reads is meaningful
    let expected_vcf = vcf_assert![
        (T C) PASS M5mC=0.0,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// Edge case 5: Het denovo where ref is NOT the confounding base
///
/// ref=A, alt=C denovo, genotype=Het(A/C). Confounding=T, neither allele is T
/// → no excess correction.
#[test]
fn het_denovo_ref_not_confounding() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ A G ] Ref,
        [ A G ] OT,  // A reads (ref)
        [ A G ] OT,
        [ A G ] OT,
        [ C G ] OT,  // C reads (denovo alt)
        [ C G ] OT,
        [ C G ] OT,
        [ T G ] OT,  // T reads (methylation evidence)
        [ T G ] OT,
        [ T G ] OT,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C); // C is a real variant creating denovo CpG
    set_fail(&mut records[0], T); // T is methylation evidence
    let records = reprocess(records)?;

    // Genotype should be Het(A/C)
    // Beta should be 3 T / (3 C + 3 T) = 0.5
    // No excess correction because neither A nor C is the confounding base T
    let expected_vcf = vcf_assert![
        (A C) PASS M5mC=0.5 GT="0/1",
        (A T) FAIL,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// Edge case 5b: Same but with G-side denovo CpG
///
/// ref=T, alt=G denovo, genotype=Het(T/G). Confounding=A, neither allele is A
/// → no excess correction.
#[test]
fn het_denovo_ref_not_confounding_g_side() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C T ] Ref,
        [ C T ] OB,  // T reads (ref)
        [ C T ] OB,
        [ C T ] OB,
        [ C G ] OB,  // G reads (denovo alt)
        [ C G ] OB,
        [ C G ] OB,
        [ C A ] OB,  // A reads (methylation evidence)
        [ C A ] OB,
        [ C A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], G); // G is a real variant creating denovo CpG
    set_fail(&mut records[1], A); // A is methylation evidence
    let records = reprocess(records)?;

    // Genotype should be Het(T/G)
    // Beta should be 3 A / (3 G + 3 A) = 0.5
    // No excess correction because neither T nor G is the confounding base A
    let expected_vcf = vcf_assert![
        (C .) PASS,
        (T G) PASS M5mC=0.5 GT="0/1",
        (T A) FAIL,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}
