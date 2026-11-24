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
/// The reads put the two `G` alts on one chromosome (`GG`) and the two
/// reference `C`s on the other (`CC`), so no chromosome reads `C` then `G`
/// there: position 2 is in the reference CpG its `C` keeps with the `G`
/// after it, and in no de-novo one.
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
        (C G) PASS M5mC=None,
        (C G) PASS M5mC=0.5,
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
        (C G) PASS M5mC=None GT="0/1", // Position 1: its C is the C of a de-novo CpG with position 2's G
        (C G) PASS M5mC=0.0 GT="1/1", // Position 2: G on both chromosomes, after a C on one
        (C T) FAIL,
        (C G) PASS M5mC=None GT="1/1", // Position 3: G on both, after position 2's G on both: no CpG
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

/// `TA>CG` on one chromosome, half methylated, beside the reference `TA`.
/// The reference chromosome's `T` reads have `A` beside them, so β never
/// counts them and takes no half of the reads for them: β is the CpG
/// chromosome's.
#[test]
fn a_cpg_of_two_hets_on_one_chromosome_keeps_all_its_reads() -> Result<()> {
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

/// A het `T>C` before a `G` whose own het `G>A` sits on the same chromosome:
/// the chromosomes read `CA` and `TG`, neither a CpG, and the `TG` reads are
/// no methylation.
#[test]
fn a_c_whose_partner_is_replaced_on_its_chromosome_makes_no_cpg() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ T G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ C A ] OT,
        [ C A ] OT,
        [ C A ] OT,
        [ T G ] OB,
        [ T G ] OB,
        [ T G ] OB,
        [ C A ] OB,
        [ C A ] OB,
        [ C A ] OB,
    );

    let expected_vcf = vcf_assert![
        (T C) PASS GT="0|1" M5mC=None,
        (G A) PASS GT="0|1" M5mC=None,
    ];

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C);
    set_pass(&mut records[1], A);
    let records = reprocess(records)?;
    assert!(!records[0].in_cpg() && !records[1].in_cpg());

    let vcf_records = metrics_to_vcf(&records, RecordFilters::variants())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// A reference CpG with a het `C>T` on one chromosome and a het `G>A` on the
/// other: they read `TG` and `CA`, the CpG is on neither, and both positions
/// keep their record without a β.
#[test]
fn a_reference_cpg_split_across_the_chromosomes_has_no_beta() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ C A ] OT,
        [ C A ] OT,
        [ C A ] OT,
        [ T G ] OB,
        [ T G ] OB,
        [ T G ] OB,
        [ C A ] OB,
        [ C A ] OB,
        [ C A ] OB,
    );

    let expected_vcf = vcf_assert![
        (C T) PASS GT="0|1" M5mC=None,
        (G A) PASS GT="1|0" M5mC=None,
    ];

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    set_pass(&mut records[1], A);
    let records = reprocess(records)?;
    assert!(!records[0].in_cpg() && !records[1].in_cpg());

    let vcf_records = metrics_to_vcf(&records, RecordFilters::variants())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// A het `C>T` at a reference CpG beside a homozygous `G`: the `T`
/// chromosome's reads have the `G` beside them and look methylated on OT, so
/// β there takes half the reads as its. On OB each read shows which
/// chromosome it is from, and the `G`'s β counts the CpG chromosome's alone.
#[test]
fn a_het_c_to_t_confounds_the_c_and_leaves_the_g_alone() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ C G ] OT,
        [ C G ] OT,
        [ T G ] OB,
        [ T G ] OB,
        [ T G ] OB,
        [ T G ] OB,
        [ C A ] OB,
        [ C A ] OB,
        [ C G ] OB,
        [ C G ] OB,
    );

    let expected_vcf = vcf_assert![
        (C T) PASS GT="0/1" M5mC=0.5,
        (G .) PASS M5mC=0.5,
        (G A) FAIL,
    ];

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    let records = reprocess(records)?;

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// A het `T>C` and a het `A>G` on different chromosomes: they read `CA` and
/// `TG`, so the two alts make no CpG.
#[test]
fn two_alts_on_different_chromosomes_make_no_cpg() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ T A ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ C A ] OT,
        [ C A ] OT,
        [ C A ] OT,
        [ T G ] OB,
        [ T G ] OB,
        [ T G ] OB,
        [ C A ] OB,
        [ C A ] OB,
        [ C A ] OB,
    );

    let expected_vcf = vcf_assert![
        (T C) PASS GT="0|1" M5mC=None,
        (A G) PASS GT="1|0" M5mC=None,
    ];

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C);
    set_pass(&mut records[1], G);
    let records = reprocess(records)?;
    assert!(!records[0].in_cpg() && !records[1].in_cpg());

    let vcf_records = metrics_to_vcf(&records, RecordFilters::variants())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// A homozygous `T>C` before a het `A>G`: the chromosomes read `CG`, half
/// methylated, and `CA`. The `CA` reads have no `G` beside the `C`, so the
/// `C`'s β is the CpG chromosome's; on OB they read as a methylated `G`, so
/// the `G`'s β takes half the reads as theirs.
#[test]
fn a_hom_c_before_a_het_g_confounds_the_g_alone() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ T A ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ C G ] OT,
        [ C G ] OT,
        [ C A ] OT,
        [ C A ] OT,
        [ C A ] OT,
        [ C A ] OT,
        [ C A ] OB,
        [ C A ] OB,
        [ C G ] OB,
        [ C G ] OB,
        [ C A ] OB,
        [ C A ] OB,
        [ C A ] OB,
        [ C A ] OB,
    );

    let expected_vcf = vcf_assert![
        (T C) PASS GT="1/1" M5mC=0.5,
        (A G) PASS GT="0/1" M5mC=0.5,
    ];

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C);
    set_pass(&mut records[1], G);
    let records = reprocess(records)?;

    let vcf_records = metrics_to_vcf(&records, RecordFilters::variants())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}
