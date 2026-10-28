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
///   - Its reference C is the partner of position 1's de-novo CpG
///   - It also has C>G creating its own de-novo CpG
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
        (C G) PASS M5mC=None, // Position 1: C with de-novo CpG G
        (C G) PASS M5mC=vec![Some(0.5), None], // Position 2: Real variant G - original CpG only (no OB evidence for denovo)
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
        (C G) PASS M5mC=None GT="0/1", // Position 1: C with de-novo CpG G, and no methylation evidence
        (C G) PASS M5mC=vec![0.0, 0.0] GT="1/1", // Position 2: Real variant G, genotype G/G, no C allele to methylate
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
        (C .) PASS M5mC=None,    // pos0: original CpG C-side, no T reads → NoEvidence
        (G C) PASS M5mC=vec![0.0, 0.75],  // pos1: original G-side beta=0.0, de-novo C-side beta=0.75
        (G .) PASS M5mC=None,  // pos2: no OB reads with before=C → NoEvidence
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
        (T G) PASS M5mC=None,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// `TA>CG` on both chromosomes: two variants that form a CpG together, as at chr12:20,868,616
/// in NA12878. The reference context alone sees no CpG, so unless the pair is recognised the
/// methylated `C` (read `T` on OT) and `G` (read `A` on OB) count as reference reads and make
/// both sites `0/1`.
#[test]
fn two_variants_forming_a_cpg_are_a_denovo_cpg() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ T A ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ C G ] OT,
        [ C A ] OB,
        [ C A ] OB,
        [ C A ] OB,
        [ C G ] OB,
    );

    let expected_vcf = vcf_assert![
        (T C) PASS GT="1/1" M5mC=0.75,
        (A G) PASS GT="1/1" M5mC=0.75,
    ];

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C);
    set_pass(&mut records[1], G);
    let records = reprocess(records)?;

    let vcf_records = metrics_to_vcf(&records, RecordFilters::variants())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// The `C` on one chromosome and the `G` on the other (`CA` / `TG`) make no
/// CpG: the pair-HMM's reads put them in trans.
#[test]
fn two_variants_in_trans_make_no_cpg() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ T A ] Ref,
        [ C A ] OT,
        [ C A ] OT,
        [ C A ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ C A ] OB,
        [ C A ] OB,
        [ C A ] OB,
        [ T G ] OB,
        [ T G ] OB,
        [ T G ] OB,
    );

    let expected_vcf = vcf_assert![
        (T C) PASS GT="0|1" M5mC=None,
        (A G) PASS GT="1|0" M5mC=None,
    ];

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C);
    set_pass(&mut records[1], G);
    let records = reprocess(records)?;
    assert!(records.iter().all(|r| !r.in_denovo_cpg()));

    let vcf_records = metrics_to_vcf(&records, RecordFilters::variants())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// A `T>C` before a reference `G` whose `G>A` sits on the same chromosome
/// (`CA` / `TG`) makes no CpG: the `C`'s chromosome lost its `G`.
#[test]
fn a_cpg_partner_replaced_on_the_same_chromosome_makes_no_cpg() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ T G ] Ref,
        [ C A ] OT,
        [ C A ] OT,
        [ C A ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ C A ] OB,
        [ C A ] OB,
        [ C A ] OB,
        [ T G ] OB,
        [ T G ] OB,
        [ T G ] OB,
    );

    let expected_vcf = vcf_assert![
        (T C) PASS GT="0|1" M5mC=None,
        (G A) PASS GT="0|1" M5mC=None,
    ];

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C);
    set_pass(&mut records[1], A);
    let records = reprocess(records)?;
    assert!(records.iter().all(|r| !r.in_denovo_cpg()));

    let vcf_records = metrics_to_vcf(&records, RecordFilters::variants())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// `TA>CG` on one chromosome, half methylated, beside the reference `TA`.
/// The reference chromosome's `T` reads have `A` beside them, so β never
/// counts them and no het correction applies: β is the CpG chromosome's.
#[test]
fn a_cis_het_pair_cpg_is_not_corrected_for_the_reference_t() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ T A ] Ref,
        [ T A ] OT,
        [ T A ] OT,
        [ T A ] OT,
        [ T A ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ C G ] OT,
        [ C G ] OT,
        [ T A ] OB,
        [ T A ] OB,
        [ T A ] OB,
        [ T A ] OB,
        [ C A ] OB,
        [ C A ] OB,
        [ C G ] OB,
        [ C G ] OB,
    );

    let expected_vcf = vcf_assert![
        (T C) PASS GT="0|1" M5mC=0.5,
        (A G) PASS GT="0|1" M5mC=0.5,
    ];

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C);
    set_pass(&mut records[1], G);
    let records = reprocess(records)?;

    let vcf_records = metrics_to_vcf(&records, RecordFilters::variants())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}
