use super::*;

#[test]
fn fully_methylated() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [C G] OB,
    );
    let metrics = to_metrics(&ps[0], &seg, None);
    let result = call(&metrics).unwrap();

    assert_original(result, 1.0);
}

#[test]
fn unmethylated() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [C G] OT,
        [C G] OT,
        [C G] OT,
    );
    let metrics = to_metrics(&ps[0], &seg, None);
    let result = call(&metrics).unwrap();

    assert_original(result, 0.0);
}

#[test]
fn partially_methylated() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [C G] OT,
        [C G] OT,
    );
    let metrics = to_metrics(&ps[0], &seg, None);
    let result = call(&metrics).unwrap();

    assert_original(result, 0.6);
}

#[test]
fn het_snp() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [C G] OT,
    );

    let metrics = to_metrics(&ps[0], &seg, Some(ct_genotype()));
    let result = call(&metrics).unwrap();

    assert_original(result, 0.71428571);
}

#[test]
fn het_snp_fifty_fifty() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
    );
    let metrics = to_metrics(&ps[0], &seg, Some(ct_genotype()));
    let result = call(&metrics).unwrap();

    assert_original(result, 0.0);
}

#[test]
fn no_evidence() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [C G] OB,
    );
    let metrics = to_metrics(&ps[0], &seg, None);
    let result = call(&metrics).unwrap();

    assert_none(result);
}

#[test]
fn no_alt_t() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [C G] OT,
        [C G] OT,
    );
    let metrics = to_metrics(&ps[0], &seg, None);
    let result = call(&metrics).unwrap();

    assert_original(result, 0.0);
}

/// `1/2` of two alts: no chromosome keeps the reference `C`, so the
/// original CpG exists on neither, as under a homozygous alt.
#[test]
fn two_alts_leave_no_reference_cpg() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [A G] OT,
        [A G] OT,
        [A G] OT,
        [C G] OT,
    );
    let gt = EstimatedGenotype {
        genotype: GenotypeTag::alt_het(
            std::num::NonZeroU8::MIN,
            std::num::NonZeroU8::new(2).expect("two"),
        ),
        quality: Phred::from_phred(20),
        source: GenotypeSource::Binomial,
    };
    let metrics = to_metrics(&ps[0], &seg, Some(gt));
    let result = call(&metrics).unwrap();

    assert!(result.is_empty(), "no chromosome keeps the CpG: {result:?}");
}
