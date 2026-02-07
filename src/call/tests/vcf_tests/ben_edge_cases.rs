//! Edge cases identified by Ben: CGG and CAG contexts with variants at middle position
//!
//! These tests document complex scenarios where positions have dual roles in
//! both original and denovo CpG contexts, or where het variants create different
//! denovo CpGs on each chromosome.

use crate::{call::tests::utils::*, pileups, vcf_assert};
use rastair_types::Base::*;

// =============================================================================
// CGG Context: Middle G with G→C variant
// =============================================================================
// In CGG context (C at pos 0, G at pos 1, G at pos 2):
// - Original CpG at positions 0-1 (C-side at 0, G-side at 1)
// - When position 1 has G→C variant, creates denovo CpG at positions 1-2
// - Position 1 becomes DUAL ROLE: G-side of original CpG AND C-side of denovo CpG

/// CGG with G→C real variant at middle position, no methylation
///
/// Tests the basic dual-role scenario: position 1 is both G-side of original
/// CpG (0-1) and C-side of denovo CpG (1-2).
#[test]
fn cgg_middle_g_to_c_no_methylation() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G G ] Ref,
        [ C C G ] OT,  // Position 1: C reads (real variant)
        [ C C G ] OT,
        [ C C G ] OT,
        [ C C G ] OT,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], C); // C is real variant at position 1
    let records = reprocess(records)?;

    // Position 0: C-side of original CpG, no methylation
    // Position 1: Dual role - G-side of original + C-side of denovo
    // Position 2: G-side of denovo CpG
    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=vec![0.0, 0.0],  // Position 0: no methylation on original CpG
        (G C) PASS M5mC=vec![0.0, 0.0],  // Position 1: dual role, no methylation either side
        (G .) PASS M5mC=0.0,              // Position 2: G-side of denovo, no methylation
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// CGG with G→C variant, methylation on original CpG (G-side at position 1)
///
/// Original CpG at 0-1: methylation shows as A reads on OB strand at position 1
#[test]
fn cgg_middle_g_to_c_methylation_on_original_cpg() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G G ] Ref,
        [ C C G ] OT,  // Position 1: C reads (real variant)
        [ C C G ] OT,
        [ C C G ] OT,
        [ C A G ] OB,  // Position 1: A reads (methylation on original G-side)
        [ C A G ] OB,
        [ C A G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], C); // C is real variant
    set_fail(&mut records[1], A); // A is methylation evidence
    let records = reprocess(records)?;

    // Position 1 shows methylation on ORIGINAL CpG (G-side)
    let expected_vcf = vcf_assert![
        (C .) PASS,
        (G C) PASS M5mC=vec![1.0, 0.0],  // Original: 1.0, Denovo: 0.0
        (G A) FAIL,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// CGG with G→C variant, methylation on denovo CpG (C-side at position 1)
///
/// Denovo CpG at 1-2: methylation shows as T reads on OT strand at position 1
#[test]
fn cgg_middle_g_to_c_methylation_on_denovo_cpg() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G G ] Ref,
        [ C T G ] OT,  // Position 1: T reads (methylation on denovo C-side)
        [ C T G ] OT,
        [ C T G ] OT,
        [ C C G ] OT,  // Position 1: C reads (real variant)
        [ C C G ] OT,
        [ C C G ] OT,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], C); // C is real variant
    set_fail(&mut records[1], T); // T is methylation evidence
    let records = reprocess(records)?;

    // Position 1 shows methylation on DENOVO CpG (C-side)
    let expected_vcf = vcf_assert![
        (C .) PASS,
        (G C) PASS M5mC=vec![0.0, 0.5],  // Original: 0.0, Denovo: 0.5
        (G T) FAIL,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// CGG with G→C variant, methylation on BOTH original and denovo CpGs
///
/// Shows both A reads (original G-side) and T reads (denovo C-side)
#[test]
fn cgg_middle_g_to_c_methylation_on_both() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G G ] Ref,
        [ C T G ] OT,  // Position 1: T reads (denovo C-side methylation)
        [ C T G ] OT,
        [ C C G ] OT,  // Position 1: C reads (real variant)
        [ C C G ] OT,
        [ C A G ] OB,  // Position 1: A reads (original G-side methylation)
        [ C A G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], C); // C is real variant
    set_fail(&mut records[1], T); // T is methylation evidence
    set_fail(&mut records[1], A); // A is methylation evidence
    let records = reprocess(records)?;

    // Position 1 shows methylation on BOTH CpG contexts
    let expected_vcf = vcf_assert![
        (C .) PASS,
        (G C) PASS M5mC=vec![1.0, 0.5],  // Original: 1.0, Denovo: 0.5
        (G T) FAIL,
        (G A) FAIL,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// CGG with Het(G/C) at middle position
///
/// Tests genotype Het(G/C) with methylation
#[test]
fn cgg_middle_het_g_c_with_methylation() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G G ] Ref,
        [ C G G ] OT,  // Position 1: G reads (ref allele)
        [ C G G ] OT,
        [ C G G ] OT,
        [ C C G ] OT,  // Position 1: C reads (alt allele)
        [ C C G ] OT,
        [ C C G ] OT,
        [ C T G ] OT,  // Position 1: T reads (methylation on C allele)
        [ C T G ] OT,
        [ C A G ] OB,  // Position 1: A reads (methylation on G allele)
        [ C A G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], C); // C is real variant
    set_fail(&mut records[1], T); // T is methylation evidence
    set_fail(&mut records[1], A); // A is methylation evidence
    let records = reprocess(records)?;

    // Position 1: Het(G/C) with methylation on both alleles
    let expected_vcf = vcf_assert![
        (C .) PASS,
        (G C) PASS GT="0/1",  // Should have both M5mC values
        (G T) FAIL,
        (G A) FAIL,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// CGG with HomAlt(C/C) at middle position
///
/// Tests HomAlt genotype where all chromosomes have the variant
#[test]
fn cgg_middle_homalt_c_with_methylation() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G G ] Ref,
        [ C C G ] OT,  // Position 1: All C reads (HomAlt)
        [ C C G ] OT,
        [ C C G ] OT,
        [ C C G ] OT,
        [ C C G ] OT,
        [ C C G ] OT,
        [ C T G ] OT,  // Position 1: T reads (methylation on C allele)
        [ C T G ] OT,
        [ C T G ] OT,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], C); // C is real variant
    set_fail(&mut records[1], T); // T is methylation evidence
    let records = reprocess(records)?;

    // Position 1: HomAlt(C/C)
    // Original CpG: beta=0.0 (no G allele to methylate)
    // Denovo CpG: normal beta calculation
    let expected_vcf = vcf_assert![
        (C .) PASS,
        (G C) PASS M5mC=vec![0.0, 0.333333] GT="1/1",  // Original: 0.0, Denovo: 3/(6+3)
        (G T) FAIL,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

// =============================================================================
// CAG Context: Middle A with het variants creating different denovo CpGs
// =============================================================================
// In CAG context (C at pos 0, A at pos 1, G at pos 2):
// - No original CpG (CA and AG are not CpGs)
// - If position 1 has A→C variant: creates denovo CpG at 1-2 (C-side at 1, G-side at 2)
// - If position 1 has A→G variant: creates denovo CpG at 0-1 (C-side at 0, G-side at 1)
// - With Het(C/G): BOTH denovo CpGs exist, on different chromosomes!

/// CAG with A→C variant at middle position
///
/// Creates denovo CpG at positions 1-2 (C at 1, G at 2)
#[test]
fn cag_middle_a_to_c_creates_denovo_at_1_2() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C A G ] Ref,
        [ C C G ] OT,  // Position 1: C reads (creates denovo CpG at 1-2)
        [ C C G ] OT,
        [ C C G ] OT,
        [ C T G ] OT,  // Position 1: T reads (methylation on C-side)
        [ C T G ] OT,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], C); // C is real variant
    set_fail(&mut records[1], T); // T is methylation evidence
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS,
        (A C) PASS M5mC=0.4,  // Denovo CpG at 1-2, methylation on C-side
        (A T) FAIL,
        (G .) PASS M5mC=0.0,  // G-side of denovo CpG
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// CAG with A→G variant at middle position
///
/// Creates denovo CpG at positions 0-1 (C at 0, G at 1)
#[test]
fn cag_middle_a_to_g_creates_denovo_at_0_1() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C A G ] Ref,
        [ C G G ] OB,  // Position 1: G reads (creates denovo CpG at 0-1)
        [ C G G ] OB,
        [ C G G ] OB,
        [ C A G ] OB,  // Position 1: A reads (methylation on G-side)
        [ C A G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], G); // G is real variant
    set_fail(&mut records[1], A); // A is methylation evidence
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=0.0,  // C-side of denovo CpG
        (A G) PASS M5mC=0.4,  // Denovo CpG at 0-1, methylation on G-side
        (A A) FAIL,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// CAG with Het(A/C) at middle position
///
/// One chromosome: A (no denovo)
/// Other chromosome: C (denovo CpG at 1-2)
#[test]
fn cag_middle_het_a_c() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C A G ] Ref,
        [ C A G ] OT,  // Position 1: A reads (ref allele)
        [ C A G ] OT,
        [ C A G ] OT,
        [ C C G ] OT,  // Position 1: C reads (alt allele, denovo at 1-2)
        [ C C G ] OT,
        [ C C G ] OT,
        [ C T G ] OT,  // Position 1: T reads (methylation on C)
        [ C T G ] OT,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], C); // C is real variant
    set_fail(&mut records[1], T); // T is methylation evidence
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS,
        (A C) PASS GT="0/1",  // Het with denovo CpG on C allele
        (A T) FAIL,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// CAG with Het(A/G) at middle position
///
/// One chromosome: A (no denovo)
/// Other chromosome: G (denovo CpG at 0-1)
///
/// NOTE: Methylation on G-side would show as A reads, but A is also the ref base,
/// so this creates a confounding situation. This test shows the genotype call
/// without attempting to distinguish methylation from ref allele.
#[test]
fn cag_middle_het_a_g() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C A G ] Ref,
        [ C A G ] OB,  // Position 1: A reads (ref allele - confounded with methylation!)
        [ C A G ] OB,
        [ C A G ] OB,
        [ C A G ] OB,
        [ C A G ] OB,
        [ C G G ] OB,  // Position 1: G reads (alt allele, denovo at 0-1)
        [ C G G ] OB,
        [ C G G ] OB,
        [ C G G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], G); // G is real variant
    let records = reprocess(records)?;

    // A reads are confounded: could be ref allele OR methylation evidence
    // The genotype estimator should see Het(A/G)
    let expected_vcf = vcf_assert![
        (C .) PASS,
        (A G) PASS GT="0/1",  // Het with denovo CpG on G allele
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// CAG with Het(C/G) at middle position - THE EXTREMELY WEIRD CASE
///
/// One chromosome: C (denovo CpG at 1-2)
/// Other chromosome: G (denovo CpG at 0-1)
/// BOTH denovo CpGs exist simultaneously, but at DIFFERENT positions!
#[test]
fn cag_middle_het_c_g_dual_denovo() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C A G ] Ref,
        [ C C G ] OT,  // Position 1: C reads (denovo CpG at 1-2)
        [ C C G ] OT,
        [ C C G ] OT,
        [ C G G ] OB,  // Position 1: G reads (denovo CpG at 0-1)
        [ C G G ] OB,
        [ C G G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], C); // C is real variant
    set_pass(&mut records[1], G); // G is real variant
    let records = reprocess(records)?;

    // This is the weird case: TWO different denovo CpGs
    // Position 0: C-side of denovo CpG (0-1) from G allele
    // Position 1: Dual denovo role:
    //   - C-side of denovo CpG (1-2) from C allele
    //   - G-side of denovo CpG (0-1) from G allele
    // Position 2: G-side of denovo CpG (1-2) from C allele
    let expected_vcf = vcf_assert![
        (C .) PASS,  // Should have denovo methylation info from position 0-1
        (A C) PASS GT="1/2",  // Compound het: both alleles create different denovo CpGs
        (A G) PASS GT="1/2",
        (G .) PASS,  // Should have denovo methylation info from position 1-2
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// CAG with Het(C/G), methylation on C allele (denovo at 1-2)
#[test]
fn cag_middle_het_c_g_methylation_on_c_allele() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C A G ] Ref,
        [ C C G ] OT,  // Position 1: C reads (denovo CpG at 1-2)
        [ C C G ] OT,
        [ C C G ] OT,
        [ C T G ] OT,  // Position 1: T reads (methylation on C allele)
        [ C T G ] OT,
        [ C G G ] OB,  // Position 1: G reads (denovo CpG at 0-1)
        [ C G G ] OB,
        [ C G G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], C); // C is real variant
    set_pass(&mut records[1], G); // G is real variant
    set_fail(&mut records[1], T); // T is methylation evidence on C allele
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS,
        (A C) PASS GT="1/2",  // Methylation on C allele (denovo at 1-2)
        (A G) PASS GT="1/2",
        (A T) FAIL,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// CAG with Het(C/G), methylation on G allele (denovo at 0-1)
#[test]
fn cag_middle_het_c_g_methylation_on_g_allele() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C A G ] Ref,
        [ C C G ] OT,  // Position 1: C reads (denovo CpG at 1-2)
        [ C C G ] OT,
        [ C C G ] OT,
        [ C G G ] OB,  // Position 1: G reads (denovo CpG at 0-1)
        [ C G G ] OB,
        [ C G G ] OB,
        [ C A G ] OB,  // Position 1: A reads (methylation on G allele)
        [ C A G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], C); // C is real variant
    set_pass(&mut records[1], G); // G is real variant
    set_fail(&mut records[1], A); // A is methylation evidence on G allele
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS,
        (A C) PASS GT="1/2",
        (A G) PASS GT="1/2",  // Methylation on G allele (denovo at 0-1)
        (A A) FAIL,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// CAG with Het(C/G), methylation on BOTH alleles
///
/// Most complex case: methylation on both denovo CpGs at different positions
/// NOTE: A reads (methylation on G allele) are confounded with the ref base A,
/// so we can't cleanly separate them. This test focuses on T reads (methylation
/// on C allele) which are unambiguous.
#[test]
fn cag_middle_het_c_g_methylation_on_both_alleles() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C A G ] Ref,
        [ C C G ] OT,  // Position 1: C reads (denovo CpG at 1-2)
        [ C C G ] OT,
        [ C T G ] OT,  // Position 1: T reads (methylation on C allele)
        [ C T G ] OT,
        [ C G G ] OB,  // Position 1: G reads (denovo CpG at 0-1)
        [ C G G ] OB,
        [ C A G ] OB,  // Position 1: A reads (confounded: ref OR methylation on G allele)
        [ C A G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], C); // C is real variant
    set_pass(&mut records[1], G); // G is real variant
    set_fail(&mut records[1], T); // T is methylation evidence on C allele
    // NOTE: Not setting A as fail because it's the ref base
    let records = reprocess(records)?;

    // Both denovo CpGs exist, but only C allele methylation is unambiguous
    // Position 0: C-side of denovo 0-1 (from G allele)
    // Position 1: dual denovo role
    // Position 2: G-side of denovo 1-2 (from C allele)
    let expected_vcf = vcf_assert![
        (C .) PASS,
        (A C) PASS GT="1/2",  // Methylation on C allele visible via T reads
        (A G) PASS GT="1/2",  // G allele methylation confounded with ref
        (A T) FAIL,
        (G .) PASS,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}
