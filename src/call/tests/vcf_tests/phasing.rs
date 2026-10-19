//! What `--phase` puts in the VCF: `GT` with a `|` separator and a `PS` tag.
//!
//! These tests set the phase call by hand, which also lets them state cases
//! the solver cannot produce, such as a call on a homozygous genotype. The
//! solver itself is tested end to end in `tests/phasing_cli.rs`.

use crate::{
    call::{
        phasing::{ColumnPhase, fixtures::call},
        tests::utils::*,
    },
    metrics::PileupMetrics,
    pileups, vcf_assert,
};
use seqair_types::Base::*;

fn with_call(record: &mut PileupMetrics, set: u32, first: u8) {
    record.phase = Some(Box::new(ColumnPhase::Phased(call(set, first))));
}

/// A `0/1` het, unphased.
fn het_record() -> Result<Vec<PileupMetrics>> {
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ G ] OT,
        [ G ] OT,
        [ G ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], G);
    reprocess(records)
}

/// A phased `GT` lists haplotype 1's allele first: `0|1` where it carries the
/// reference, `1|0` where it carries the alt.
#[test]
fn a_phased_het_lists_haplotype_one_first() -> Result<()> {
    for (first, gt) in [(0, "0|1"), (1, "1|0")] {
        let mut records = het_record()?;
        with_call(&mut records[0], 1001, first);

        let expected_vcf = vcf_assert![
            (A G) PASS GT=gt PS="1001",
        ];
        expected_vcf.matches(metrics_to_vcf(&records, RecordFilters::all())?)?;
    }

    Ok(())
}

/// Without a phase call the record is unphased and carries no `PS` at all,
/// rather than an empty one.
#[test]
fn an_unphased_het_keeps_its_slash_and_no_phase_set() -> Result<()> {
    let expected_vcf = vcf_assert![
        (A G) PASS GT="0/1" PS=".",
    ];
    expected_vcf.matches(metrics_to_vcf(&het_record()?, RecordFilters::all())?)?;

    Ok(())
}

/// `0|0` and `n|n` say nothing a phase set could add, so a stray call on one is
/// ignored rather than written out.
#[test]
fn a_homozygous_genotype_is_never_phased() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ G ] OT,
        [ G ] OT,
        [ G ] OB,
        [ G ] OB,
    );
    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], G);
    let mut records = reprocess(records)?;
    with_call(&mut records[0], 1001, 1);

    let expected_vcf = vcf_assert![
        (A G) PASS GT="1/1" PS=".",
    ];
    expected_vcf.matches(metrics_to_vcf(&records, RecordFilters::all())?)?;

    Ok(())
}

/// A compound het is why the phase call names an allele rather than saying
/// "haplotype 1 is the reference": neither of its alleles is.
#[test]
fn a_compound_het_phases_its_two_alts() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ C ] OT,
        [ C ] OT,
        [ G ] OB,
        [ G ] OB,
    );
    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C);
    set_pass(&mut records[0], G);
    let mut records = reprocess(records)?;
    assert!(
        records[0].pos_metrics.extended.genotype.expect("genotype").genotype.is_heterozygous(),
        "the fixture is meant to produce a compound het"
    );
    with_call(&mut records[0], 1001, 2);

    let expected_vcf = vcf_assert![
        (A C, G) PASS GT="2|1" PS="1001",
    ];
    expected_vcf.matches(metrics_to_vcf(&records, RecordFilters::all())?)?;

    Ok(())
}

/// The solver names alleles by their index in `pileup.alts`, which also holds
/// alts that were not called; the VCF numbers only the called ones. Here the
/// rejected `T`, the most frequent alt, takes slot 1, so the called `G` is
/// slot 2 internally and `1` in the VCF, and haplotype 1 carries it.
#[test]
fn a_phase_call_is_renumbered_like_the_genotype() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ A ] OT,
        [ A ] OB,
        [ A ] OT,
        [ G ] OT,
        [ G ] OB,
        [ G ] OB,
        [ T ] OT,
        [ T ] OB,
        [ T ] OT,
        [ T ] OB,
    );
    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_fail(&mut records[0], T);
    set_pass(&mut records[0], G);
    let mut records = reprocess(records)?;
    assert_eq!(
        records[0].alts.iter().map(|alt| alt.base).collect::<Vec<_>>(),
        [T, G],
        "the uncalled alt must come first, or no renumbering happens"
    );
    with_call(&mut records[0], 1001, 2);

    let expected_vcf = vcf_assert![
        (A G) PASS GT="1|0" PS="1001",
    ];
    expected_vcf.matches(metrics_to_vcf(&records, RecordFilters::variants())?)?;

    Ok(())
}
