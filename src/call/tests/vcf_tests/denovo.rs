use crate::{call::tests::utils::*, pileups, vcf_assert};
use seqair_types::Base::*;

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

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
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
        (G C) PASS M5mC=1.,  // Real variant C, methylation info preserved in M5mC
        (G T) FAIL,
        (G .) PASS M5mC=1.,  // other de-novo CpG position
        (G A) FAIL,
    ];

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C); // set C alt to pass, creating de-novo CpG for sure
    set_fail(&mut records[0], T); // set T alt to fail, it's a methylation evidence
    set_fail(&mut records[1], A); // methylation evidence
    let records = reprocess(records)?; // to recalculate genotypes and propagate de-novo CpG flags

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// Edge case: CCC reference with heterozygous C>G at positions 1 and 2
///
/// Position 1: C>G creates de-novo CpG at (1,2)
/// Position 2: C>G creates de-novo CpG at (2,3)
/// Position 2 has dual role:
///   - It's the "matching C" for position 1's de-novo CpG (`ThisIsTheMatchingC`)
///   - It also has C>G creating its own de-novo CpG (`ThisBecomesG`)
#[test]
fn test_adjacent_denovo_cpgs_dual_role_middle_position() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C C G ] Ref,
        [ G T G ] OT,
        [ G G G ] OT,
        [ G G G ] OB,
        [ C C G ] OT,
        [ C C A ] OB,
    );

    let expected_vcf = vcf_assert![
        (C G) PASS M5mC=0., // Position 1: C with de-novo CpG G
        (C G) PASS M5mC=0.5, // Position 2: Real variant G - original CpG only (no OB evidence for denovo)
        (C T) FAIL,
        (G .) PASS M5mC=1.0,
        (G A) FAIL,
    ];

    let records = test_call(segment, pileups, RecordFilters::cpgs())?;
    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_adjacent_denovo_cpgs_dual_role_middle_position2() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C C C ] Ref,
        [ G G G ] OT,
        [ G G G ] OT,
        [ C T G ] OT,
        [ G C A ] OB,
        [ C C A ] OB,
        [ C G A ] OB,
    );

    let expected_vcf = vcf_assert![
        (C G) PASS M5mC=0.0 GT="0/1", // Position 1: C with de-novo CpG G, and no methylation evidence
        (C G) PASS M5mC=vec![1.0, 0.0] GT="0/1", // Position 2: Real variant G, genotype G/G, no C allele to methylate
        (C T) FAIL,
        (C G) PASS M5mC=1.0 GT="1/1", // Position 3: Real variant G, methylation info for that G>A
        (C A) FAIL,
    ];

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], G);
    set_pass(&mut records[1], G);
    set_fail(&mut records[1], T); // methylation evidence
    set_pass(&mut records[2], G); // real variant
    set_fail(&mut records[2], A); // methylation evidence
    let records = reprocess(records)?; // to recalculate genotypes and propagate de-novo CpG flags

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_a_and_t_cant_form_denovo() -> Result<()> {
    // Test that A and T reference bases don't produce methylation evidence rows
    // Currently, de-novo CpG creation is only considered when one of the C or G is in the reference
    let (segment, pileups) = pileups!(
        [ A T ] Ref,
        [ C G ] OT,
        [ C G ] OT,
        [ C G ] OB,
        [ C G ] OB,
    );

    let records = test_call(segment, pileups, RecordFilters::all())?;

    let expected_vcf = vcf_assert![
        (A C) PASS M5mC=None,
        (T G) PASS M5mC=None,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// C[G>C]G scenario: position 1 (G-ref) is both the G-side of the original CpG (C-G at
/// positions 0-1) and the new de-novo CpG C-side (G->C SNP creates C-G at positions 1-2).
#[test]
fn test_ref_cpg_g_to_c_snp_no_original_methylation() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G G ] Ref,
        [ C T G ] OT,  // T at pos1 = methylated new C after G->C SNP (TAPS: mC → T on OT)
        [ C T G ] OT,
        [ C T G ] OT,
        [ C C G ] OT,  // C at pos1 = unmethylated new C (the SNP allele, unmodified)
        [ C G A ] OB,  // G at pos1 = original G ref (no methylation), A at pos2 = methylated new G-side
        [ C G A ] OB,
        [ C G A ] OB,
        [ C G G ] OB,  // G at pos2 = unmethylated new G-side
    );

    let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
    set_pass(&mut records[1], C); // G->C at pos1 is a real variant creating a de-novo CpG
    set_fail(&mut records[1], T); // T at pos1 is methylation evidence for the new C-side
    set_fail(&mut records[2], A); // A at pos2 is methylation evidence for the new G-side
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=0.,    // pos0: original CpG C-side, no T reads → NoEvidence
        (G C) PASS M5mC=vec![0.0, 0.75],  // pos1: original G-side beta=0.0, de-novo C-side beta=0.75
        (G .) PASS M5mC=0.,  // pos2: no OB reads with before=C → NoEvidence
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::cpgs())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_cpg_that_is_also_denovo() -> Result<()> {
    // Test that A and T reference bases don't produce methylation evidence rows
    // Currently, de-novo CpG creation is only considered when one of the C or G is in the reference
    let (segment, pileups) = pileups!(
        [ C C G ] Ref,
        [ T G G ] OT,
        [ T G G ] OT,
        [ C T G ] OT,
        [ C T G ] OT,
        [ C C A ] OB,
        [ C C A ] OB,
        [ C A G ] OB,
        [ C A G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::cpgs())?;
    set_fail(&mut records[0], T); // make this methylation evidence
    set_pass(&mut records[1], G); // make this G in a denovo
    set_fail(&mut records[1], A); // make this methylation evidence
    set_fail(&mut records[2], A); // make this methylation evidence
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=1.0,
        (C G) PASS M5mC=vec![1.0, 1.0], // Real variant G - both original and de-novo CpG beta values
        (G .) PASS M5mC=1.0,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::cpgs())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_rejected_denovo_candidate_reports_no_methylation() -> Result<()> {
    // pos0 is a plain CpH C; pos1 has a low-VAF T>G that would create a CpG.
    let (segment, pileups) = pileups!(
        [ C T ] Ref,
        [ C T ] OT,
        [ C T ] OT,
        [ T G ] OT,
        [ C T ] OB,
        [ C T ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_fail(&mut records[1], G);
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C T) FAIL M5mC=None,
        (T G) FAIL M5mC=None,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_accepted_denovo_candidate_reports_methylation() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C T ] Ref,
        [ C T ] OT,
        [ C T ] OT,
        [ T G ] OT,
        [ C T ] OB,
        [ C T ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[1], G);
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=1.0,  // matching C of the new CpG: the single TG read is methylated
        (C T) FAIL,
        (T G) PASS M5mC=0.0,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}
