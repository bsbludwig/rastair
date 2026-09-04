use crate::{call::tests::utils::*, pileups, vcf_assert};
use seqair_types::Base::*;

/// A non-CpG variant emits a record with no methylation context. The M5mC
/// FORMAT field for such a record must round-trip through BCF without making
/// htslib-based float readers panic: previously it was encoded as a zero-length
/// value (`n == 0`), which makes `rust_htslib`'s `Format::float()` panic on
/// `chunks(0)`. The field must be either absent or carry at least one value.
#[test]
fn non_methylation_record_m5mc_is_readable_in_bcf() -> Result<()> {
    use rust_htslib::bcf::{self, Read as _};
    use std::io::Write as _;

    // pos1 is a non-CpG A->G variant; it has no methylation context.
    let (segment, pileups) = pileups!(
        [ A A ] Ref,
        [ A G ] OT,
        [ A G ] OT,
        [ A G ] OB,
        [ A G ] OB,
    );

    let records = test_call(segment, pileups, RecordFilters::all())?;
    let bcf_bytes = metrics_to_bcf(&records, RecordFilters::all())?;

    let mut tmp = tempfile::NamedTempFile::new()?;
    tmp.write_all(&bcf_bytes)?;
    tmp.flush()?;

    let mut reader = bcf::Reader::from_path(tmp.path())?;
    let mut saw_record = false;
    for rec in reader.records() {
        let rec = rec?;
        saw_record = true;
        // Must not panic. A present M5mC field must have >= 1 value per sample.
        if let Ok(values) = rec.format(b"M5mC").float() {
            for sample in values.iter() {
                assert!(
                    !sample.is_empty(),
                    "M5mC must not be a zero-length FORMAT value at pos {}",
                    rec.pos() + 1
                );
            }
        }
    }
    assert!(saw_record, "expected at least one emitted record");

    Ok(())
}

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
        // (A .) PASS, // FIXME: include a?
        (T C) PASS,
    ];

    let records = test_call(segment, pileups, RecordFilters::all())?;
    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
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
    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    // Now test with CpG filter: Should match the middle two positions
    // so same as above
    let records = test_call(segment.clone(), pileups.clone(), RecordFilters::cpgs())?;
    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=0.0,
        (G .) PASS M5mC=0.0,
    ];
    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    // Now with `--all`: only the CpG positions are emitted. The reference-only
    // non-CpG positions (A and T) must stay suppressed even under `--all` —
    // emitting them would surface a bare M5mC value without the CPG/CPGnovo
    // tags set.
    let records = test_call(segment, pileups, RecordFilters::all())?;
    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=0.0,
        (G .) PASS M5mC=0.0,
    ];
    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
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

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
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

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
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

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
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

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_single_end_like_variant_calling() -> Result<()> {
    // Single-end input has no pair flags; orientation comes from read direction.
    // At this test layer, all reads are effectively "first in pair".
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ C ] OT,
        [ C ] OT,
        [ C ] OB,
        [ C ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C);
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (A C) PASS GT="1/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}
