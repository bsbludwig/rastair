//! Tests for genotyping with ML scores
use crate::utils::MappingQuality;
use crate::{
    call::tests::pileup::{SimpleRead, SimpleReads},
    metrics::PileupMetrics,
};
use crate::{call::tests::utils::*, call::variant_calling::GenotypeTag, pileups, vcf_assert};
use seqair_types::{Base, Base::*, Probability, Strand};
use std::num::NonZeroU8;

/// Three alt reads and one reference read of any substitution TAPS does not
/// confound: the engine calls the het.
#[test]
fn a_called_substitution_is_a_het_at_three_reads_in_four() -> Result<()> {
    for (reference, alt) in
        [(A, T), (A, G), (A, C), (T, A), (T, G), (T, C), (C, A), (C, G), (G, T), (G, C)]
    {
        let reads =
            [(alt, Strand::OT), (alt, Strand::OT), (alt, Strand::OB), (reference, Strand::OB)]
                .map(|(base, strand)| ReadLine { bases: vec![base], strand });
        let (segment, pileups) = create_pileups(vec![reference], reads.into());
        let mut records = test_call(segment, pileups, RecordFilters::all())?;
        set_pass(&mut records[0], alt);
        let records = reprocess(records)?;
        assert_eq!(genotype_at(&records, 0), Some(GenotypeTag::CT), "{reference:?}>{alt:?}");
    }
    Ok(())
}

#[test]
fn test_c_to_t_outside_cpg_counts_both_strands() -> Result<()> {
    // Outside a CpG a T over a C is no methylation, so both strands' Ts are the T allele.
    let (segment, pileups) = pileups!(
        [ C ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ C ] OB,
        [ C ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    let records = reprocess(records)?;

    // 2 T reads, 2 C reads -> heterozygous (0/1)
    let expected_vcf = vcf_assert![
        (C T) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_multi_allelic_site_picks_highest_ml() -> Result<()> {
    // Two called alts at one site are written in one record.
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ G ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T); // ML = 1.0

    // Set G to have lower ML but still passing
    let g_filters = records[0].find_alt_mut(G).unwrap();
    g_filters.verdict = crate::metrics::AltVerdict::Scored(Probability::new_panicky(0.85));
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (A T,G) PASS,  // Both alts passing, combined in one row
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_all_alts_below_ml_threshold() -> Result<()> {
    // All alts have ML scores below threshold - should call HomRef (0/0)
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ G ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_fail(&mut records[0], T); // ML = 0.0
    set_fail(&mut records[0], G); // ML = 0.0
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (A T) FAIL GT="0/0",
        (A G) FAIL GT="0/0",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_het_vs_hom_alt_genotyping_for_a_base() -> Result<()> {
    // Two reads of each allele are a het.
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ A ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (A T) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_het_genotyping_with_unbalanced_reads() -> Result<()> {
    // Three alt reads against two reference reads are a het.
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ A ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    let records = reprocess(records)?;

    // 3 T reads, 2 A reads -> heterozygous (0/1)
    let expected_vcf = vcf_assert![
        (A T) PASS GT="0/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_hom_alt_genotyping_for_t_base() -> Result<()> {
    // T→G with mostly alt reads suggesting homozygous alt genotype (1/1)
    let (segment, pileups) = pileups!(
        [ T ] Ref,
        [ G ] OT,
        [ G ] OT,
        [ G ] OB,
        [ G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    // 4 G reads, 0 T reads -> homozygous alt (1/1)
    let expected_vcf = vcf_assert![
        (T G) PASS GT="1/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_compound_het_with_balanced_alts() -> Result<()> {
    // A reference with two alt alleles (T and G), both passing ML threshold
    // with roughly balanced read support should produce compound het (1/2)
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ G ] OT,
        [ G ] OT,
        [ G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    // 3 T reads, 3 G reads, 0 A reads -> compound heterozygous (1/2)
    let expected_vcf = vcf_assert![
        (A T,G) PASS GT="1/2",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    // One `GQ` over the six genotypes of `A`, `T`, `G`, and their `PL` in
    // `GL` order, `1/2` (the fifth) the likeliest.
    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    let [record] = &vcf_records[..] else { panic!("one record: {vcf_records:?}") };
    let pl: Vec<i32> =
        record.format.get("PL").expect("a PL").split(',').map(|p| p.parse().expect("PL")).collect();
    assert_eq!(pl.len(), 6, "{pl:?}");
    assert_eq!(pl.get(4), Some(&0), "{pl:?}");
    let gq: i32 = record.format.get("GQ").expect("a GQ").parse()?;
    assert!(pl.iter().enumerate().filter(|&(i, _)| i != 4).all(|(_, &p)| p >= gq), "{pl:?} {gq}");

    Ok(())
}

#[test]
fn test_compound_het_with_slightly_unbalanced_alts() -> Result<()> {
    // Two called alts and no reference read: the compound het.
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ T ] OB,
        [ G ] OT,
        [ G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (A T,G) PASS GT="1/2",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_two_alts_passing_but_one_dominant() -> Result<()> {
    // Both alts are called, so the pair carries both: one reference read
    // among seven is the engine's artefact, not a chromosome.
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ T ] OB,
        [ T ] OB,
        [ G ] OT,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (A T,G) PASS GT="1/2",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_compound_het_vs_single_het_with_ref() -> Result<()> {
    // Two alts the model passes, but three reference reads: the pair carries
    // the reference and the alt with more reads, and the other is not called.
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ G ] OT,
        [ G ] OB,
        [ A ] OT,
        [ A ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    let expected_vcf = vcf_assert![
        (A T) PASS GT="0/1",
        (A G) FAIL GT="0/0",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_compound_het_for_c_base_with_non_methylation_alts() -> Result<()> {
    // A and G over a C are not confounded by TAPS: every read counts.
    let (segment, pileups) = pileups!(
        [ C ] Ref,
        [ A ] OT,
        [ A ] OT,
        [ A ] OB,
        [ G ] OT,
        [ G ] OB,
        [ G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], A);
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    // 3 A reads, 3 G reads, 0 C reads -> compound heterozygous (1/2)
    let expected_vcf = vcf_assert![
        (C A,G) PASS GT="1/2",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

#[test]
fn test_three_alts_passing_uses_top_two() -> Result<()> {
    // Three alts the model passes: the engine scores the two it scores
    // highest, and a diploid genotype holds no third.
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ G ] OT,
        [ G ] OB,
        [ C ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;

    // Set T with highest ML, G with medium ML, C with lowest ML
    set_pass(&mut records[0], T);
    records[0].find_alt_mut(T).unwrap().verdict =
        crate::metrics::AltVerdict::Scored(Probability::new_panicky(0.95));
    set_pass(&mut records[0], G);
    records[0].find_alt_mut(G).unwrap().verdict =
        crate::metrics::AltVerdict::Scored(Probability::new_panicky(0.90));
    set_pass(&mut records[0], C);
    records[0].find_alt_mut(C).unwrap().verdict =
        crate::metrics::AltVerdict::Scored(Probability::new_panicky(0.85));
    let records = reprocess(records)?;

    // 2 T reads, 2 G reads, 1 C read -> compound het with T and G; C is
    // written apart, not called.
    let expected_vcf = vcf_assert![
        (A T,G) PASS GT="1/2",
        (A C) FAIL,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// Test that de-novo CpG genotyping ignores REF reads from converted strand.
#[test]
fn test_denovo_cpg_genotyping_ignores_converted_strand_ref_reads() -> Result<()> {
    // Position 1: A with next base G (A→C creates de-novo CpG)
    // Position 2: G (partner of de-novo CpG)
    let (segment, pileups) = pileups!(
        [ A G ] Ref,
        // OT strand: A (ref) and T (methylated C)
        [ A G ] OT,
        [ A G ] OT,
        [ T G ] OT,  // methylated
        [ T G ] OT,  // methylated
        // OB strand: C (alt) only - no ref reads
        [ C A ] OB,
        [ C A ] OB,
        [ C A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], C); // C is the real variant creating de-novo CpG
    set_fail(&mut records[0], T); // T is methylation evidence, not a real variant
    set_fail(&mut records[1], A); // A on second position is methylation evidence
    let records = reprocess(records)?;

    // Genotyping for the C alt should only use OB strand
    let expected_vcf = vcf_assert![
        (A C) PASS GT="0/1",
        (A T) FAIL,
        (G .) PASS,  // partner position of de-novo CpG
        (G A) FAIL,
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// The writer must report the genotype the pipeline stored, not one it derives
/// again for itself. `--error-model` only reaches the stored estimate, so a
/// writer that re-estimates with the default model contradicts the BED output
/// and the methylation calls, which both read `pos_metrics.genotype`.
#[test]
fn vcf_genotype_is_the_stored_estimate() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ T ] OB,
        [ A ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_pass(&mut records[0], T);
    let mut records = reprocess(records)?;

    // Stand-in for what another error model concludes from the same counts.
    let stored = records[0].pos_metrics.genotype.as_mut().expect("genotype");
    stored.genotype = GenotypeTag::HomAlt(NonZeroU8::new(1).expect("1 > 0"));

    let expected_vcf = vcf_assert![
        (A T) PASS GT="1/1",
    ];

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    expected_vcf.matches(vcf_records)?;

    Ok(())
}

/// A two-base column from `(copies, bases, strand)` groups, too deep to spell out in `pileups!`.
fn column(
    reference: [Base; 2],
    groups: &[(usize, [Base; 2], Strand)],
) -> Result<Vec<PileupMetrics>> {
    column_with(reference, groups, |_| ())
}

/// [`column`] with every read passed through `edit` first.
fn column_with(
    reference: [Base; 2],
    groups: &[(usize, [Base; 2], Strand)],
    edit: impl Fn(&mut SimpleRead),
) -> Result<Vec<PileupMetrics>> {
    let reads = groups
        .iter()
        .flat_map(|&(n, bases, strand)| {
            std::iter::repeat_with(move || ReadLine { bases: bases.to_vec(), strand }).take(n)
        })
        .collect();
    let (segment, mut pileups) = create_pileups(reference.to_vec(), reads);
    for pileup in &mut pileups {
        let mut reads = pileup.reads.to_vec();
        reads.iter_mut().for_each(&edit);
        pileup.reads = SimpleReads(reads.into());
    }
    test_call(segment, pileups, RecordFilters::all())
}

fn genotype_at(records: &[PileupMetrics], idx: usize) -> Option<GenotypeTag> {
    records.get(idx)?.pos_metrics.genotype.as_ref().map(|gt| gt.genotype)
}

#[test]
fn unmethylated_reads_on_the_converting_strand_are_reference_evidence() -> Result<()> {
    // A het C>T at a reference CpG: OB, the informative strand, shows one C in ten, which alone
    // reads as hom-alt; the five OT C reads are unmethylated CpG Cs that only the C allele can
    // produce.
    let mut records = column(
        [C, G],
        &[
            (5, [C, G], Strand::OT),
            (5, [T, G], Strand::OT),
            (1, [C, G], Strand::OB),
            (9, [T, G], Strand::OB),
        ],
    )?;
    set_pass(&mut records[0], T);
    let records = reprocess(records)?;
    assert_eq!(genotype_at(&records, 0), Some(GenotypeTag::CT));
    Ok(())
}

#[test]
fn one_unconverted_read_at_q40_is_a_c_allele() -> Result<()> {
    // A C>T at a reference CpG with one unconverted OT C among eleven at Q40: only a C allele
    // produces it, and at that quality it outweighs six OB Ts, so the engine calls the het.
    let mut records = column(
        [C, G],
        &[(1, [C, G], Strand::OT), (10, [T, G], Strand::OT), (6, [T, G], Strand::OB)],
    )?;
    set_pass(&mut records[0], T);
    let records = reprocess(records)?;
    assert_eq!(genotype_at(&records, 0), Some(GenotypeTag::CT));
    Ok(())
}

#[test]
fn unmethylated_reads_of_a_de_novo_c_are_alt_evidence() -> Result<()> {
    // A het T>C forming a CpG: OB shows one C in ten, which alone reads as hom-ref; the five OT C
    // reads are unmethylated new Cs, which only the C allele can produce.
    let mut records = column(
        [T, G],
        &[
            (5, [C, G], Strand::OT),
            (5, [T, G], Strand::OT),
            (1, [C, G], Strand::OB),
            (9, [T, G], Strand::OB),
        ],
    )?;
    set_pass(&mut records[0], C);
    let records = reprocess(records)?;
    assert_eq!(genotype_at(&records, 0), Some(GenotypeTag::CT));
    Ok(())
}

#[test]
fn low_mapq_reference_reads_do_not_make_a_hom_alt_het() -> Result<()> {
    // Two MAPQ-10 A reads among twelve G: as errors at the NovaSeq rate they are a second
    // allele, but each is misplaced one time in ten.
    let groups = [
        (1, [A, T], Strand::OT),
        (1, [A, T], Strand::OB),
        (6, [G, T], Strand::OT),
        (6, [G, T], Strand::OB),
    ];
    let mut records = column_with([A, T], &groups, |read| {
        if read.base == A {
            read.mapq = MappingQuality::from_byte(10);
        }
    })?;
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;
    assert_eq!(genotype_at(&records, 0), Some(GenotypeTag::TT));

    // The same reads at MAPQ 60 are a het.
    let mut records = column([A, T], &groups)?;
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;
    assert_eq!(genotype_at(&records, 0), Some(GenotypeTag::CT));
    Ok(())
}

/// A CpG with no read on the strand that tells a SNP from methylation has no
/// genotype to weigh: `./.`, with no `GQ`, `PL` or `QUAL`, instead of a
/// confident `0/0`.
#[test]
fn a_cpg_without_informative_reads_is_no_call() -> Result<()> {
    use crate::call::tests::fixtures::column;

    let mut site = column(1500, C, &[], None, &[]);
    site.context.after_1 = Some(G);
    site.pos_metrics.genotype = site.hom_ref_genotype(Default::default());
    let records = metrics_to_vcf(&[site], RecordFilters::all())?;
    let [record] = &records[..] else { panic!("one record: {records:?}") };
    assert_eq!(record.format.get("GT").map(String::as_str), Some("./."));
    assert_eq!(record.format.get("GQ"), None);
    assert_eq!(record.format.get("PL"), None);
    assert_eq!(record.qual, ".");
    Ok(())
}

/// A called alt behind a rejected one in the column: its record's `PL` is
/// over `A`, `G` alone, and the filtered `T` record beside it, written `0/0`,
/// carries no `GQ` or `PL`, the called genotype's not being its own.
#[test]
fn a_filtered_alt_beside_a_called_one_has_no_genotype_quality() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ A ] OT,
        [ A ] OB,
        [ T ] OT,
        [ G ] OT,
        [ G ] OT,
        [ G ] OB,
    );
    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    set_fail(&mut records[0], T);
    set_pass(&mut records[0], G);
    let records = reprocess(records)?;

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    let called = vcf_records.iter().find(|r| r.alt == ["G"]).expect("the called record");
    let pl = called.format.get("PL").expect("a PL");
    assert_eq!(pl.split(',').count(), 3, "{pl}");
    assert!(called.format.contains_key("GQ"));
    let filtered = vcf_records.iter().find(|r| r.alt == ["T"]).expect("the filtered record");
    assert_eq!(filtered.format.get("GT").map(String::as_str), Some("0/0"));
    assert_eq!(filtered.format.get("GQ"), None);
    assert_eq!(filtered.format.get("PL"), None);
    Ok(())
}

/// An alt the model did not score has no ML value, not a `0` that reads as
/// "scored, rejected".
#[test]
fn test_unscored_alt_has_missing_ml() -> Result<()> {
    let (segment, pileups) = pileups!(
        [ A ] Ref,
        [ T ] OT,
        [ T ] OT,
        [ G ] OT,
        [ G ] OB,
    );

    let mut records = test_call(segment, pileups, RecordFilters::all())?;
    records[0].find_alt_mut(T).unwrap().verdict =
        crate::metrics::AltVerdict::Scored(Probability::new_panicky(0.5));
    // An alt whose features did not compute stays unscored.
    records[0].find_alt_mut(G).unwrap().verdict = crate::metrics::AltVerdict::Unscored;
    let records = reprocess(records)?;

    let vcf_records = metrics_to_vcf(&records, RecordFilters::all())?;
    let ml = |alt: &str| {
        let record = vcf_records.iter().find(|r| r.alt == [alt]);
        record.map(|r| r.format.get("ML").map(String::as_str))
    };
    assert_eq!(ml("T"), Some(Some("0.5")), "{vcf_records:?}");
    assert_eq!(ml("G"), Some(None), "{vcf_records:?}");

    Ok(())
}
