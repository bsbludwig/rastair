//! What `--phase` puts in the VCF: `GT` with a `|` separator and a `PS` tag,
//! the haplotype split, and the methylation it resolves.
//!
//! These tests set the phase call by hand, which also lets them state cases
//! the solver cannot produce, such as a call on a homozygous genotype. The
//! solver itself is tested end to end in `tests/phasing_cli.rs`.

use crate::{
    call::{
        phasing::{
            ColumnPhase, CpgPresence, HaplotypeSplit, RecordPhases, fixtures::call,
            phase_confounded_betas,
        },
        tests::utils::*,
    },
    metrics::PileupMetrics,
    pileups,
    utils::ByHaplotype,
    vcf_assert,
};
use seqair_types::Base::*;

fn with_call(record: &mut PileupMetrics, set: u32, first: u8) {
    record.phase = Some(Box::new(ColumnPhase::Phased(RecordPhases {
        snv: Some(call(set, first)),
        ..RecordPhases::default()
    })));
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

/// A CpG whose `C` is a heterozygous `C>T`: four top-strand `T`s and two `C`s,
/// and on the bottom strand, which reads the alleles apart, two `C`s and two
/// `T`s. Unphased, the `T`s are assumed half the variant's: `M5mC` is
/// (4 - 3) / (2 + 1) at the `C`. The `G` counts only reads whose `C` is there
/// to pair with, so the variant's haplotype never dilutes it.
fn confounded_cpg() -> Result<Vec<PileupMetrics>> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ C G ] OT,
        [ C G ] OT,
        [ C A ] OB,
        [ C A ] OB,
        [ T G ] OB,
        [ T G ] OB,
    );
    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    reprocess(records)
}

/// A split in the block at 1001, with haplotype 2 lacking the CpG unless
/// `cpg` says otherwise. `counts` is `HPC`.
fn with_split(record: &mut PileupMetrics, counts: [u32; 4], cpg: [CpgPresence; 2]) {
    let [u1, c1, u2, c2] = counts;
    record.phase = Some(Box::new(ColumnPhase::Phased(RecordPhases {
        split: Some(HaplotypeSplit {
            set: call(1001, 0).set,
            unconverted: ByHaplotype::new(u1, u2),
            converted: ByHaplotype::new(c1, c2),
            cpg: ByHaplotype::new(cpg[0], cpg[1]),
            fragments: ByHaplotype::default(),
        }),
        ..RecordPhases::default()
    })));
}

const CARRIER_FIRST: [CpgPresence; 2] = [CpgPresence::Present, CpgPresence::Absent];

#[test]
fn unphased_a_confounded_cpg_assumes_half_its_converted_reads_are_the_variant() -> Result<()> {
    let mut records = confounded_cpg()?;
    phase_confounded_betas(&mut records);

    let expected_vcf = vcf_assert![
        (C T) PASS GT="0/1" M5mC=1./3. HPS=".",
        (G .) PASS M5mC=1. HPS=".",
        (G A) FAIL,
    ];
    expected_vcf.matches(metrics_to_vcf(&records, RecordFilters::all())?)?;

    Ok(())
}

/// Haplotype 1 carries the cytosine and shows two `C`s and two `T`s, so the
/// CpG is half methylated; haplotype 2's `T`s are the allele. The `G`'s split
/// would say 0.5 as well, but its beta was never confounded.
#[test]
fn phased_a_confounded_cpg_takes_its_methylation_from_the_carrier() -> Result<()> {
    let mut records = confounded_cpg()?;
    with_split(&mut records[0], [2, 2, 0, 2], CARRIER_FIRST);
    with_split(&mut records[1], [1, 1, 2, 0], CARRIER_FIRST);
    phase_confounded_betas(&mut records);

    let expected_vcf = vcf_assert![
        (C T) PASS GT="0/1" M5mC=1./2. HPS="1001" HPC="2,2,0,2" HPX="1,0",
        (G .) PASS M5mC=1. HPS="1001" HPX="1,0",
        (G A) FAIL HPS=".",  // the split is the column's, written once
    ];
    expected_vcf.matches(metrics_to_vcf(&records, RecordFilters::all())?)?;

    Ok(())
}

/// A carrier with no reads, or a split that cannot say which haplotype lost
/// the cytosine, has no methylation to report, and the unphased one stands.
#[test]
fn a_one_sided_or_unresolved_split_keeps_the_unphased_methylation() -> Result<()> {
    let unresolved = [CpgPresence::Unresolved, CpgPresence::Unresolved];
    for (counts, cpg) in [([0, 0, 0, 2], CARRIER_FIRST), ([2, 2, 0, 2], unresolved)] {
        let mut records = confounded_cpg()?;
        with_split(&mut records[0], counts, cpg);
        phase_confounded_betas(&mut records);

        let expected_vcf = vcf_assert![
            (C T) PASS GT="0/1" M5mC=1./3. HPS="1001",
            (G .) PASS M5mC=1.,
            (G A) FAIL,
        ];
        expected_vcf.matches(metrics_to_vcf(&records, RecordFilters::all())?)?;
    }

    Ok(())
}

/// Without a variant there is no assumption to replace, however differently
/// the two haplotypes are methylated.
#[test]
fn a_cpg_no_variant_confounds_keeps_its_methylation() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ C G ] Ref,
        [ T G ] OT,
        [ T G ] OT,
        [ T G ] OT,
        [ C G ] OT,
        [ C G ] OB,
        [ C G ] OB,
    );
    let mut records = reprocess(test_call(segment, pileups, RecordFilters::all())?)?;
    with_split(&mut records[0], [0, 3, 1, 0], [CpgPresence::Present, CpgPresence::Present]);
    phase_confounded_betas(&mut records);

    let expected_vcf = vcf_assert![
        (C .) PASS M5mC=3./4. HPS="1001" HPX="1,1",
        (C T) FAIL,
        (G .) PASS M5mC=0.,
    ];
    expected_vcf.matches(metrics_to_vcf(&records, RecordFilters::all())?)?;

    Ok(())
}
