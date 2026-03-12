use super::*;

#[test]
fn standard_a_to_c() {
    let (seg, ps) = pileups!(
        [A G] Ref,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
    );
    let gt = EstimatedGenotype {
        genotype: GenotypeTag::HomRef,
        likelihood: Probability::new(0.8).unwrap(),
        confidence: Probability::new(0.99).unwrap(),
    };
    let metrics = to_metrics(&ps[0], &seg, Some(gt));
    let result = call(&metrics).unwrap();

    assert_denovo(result, 0.5);
}

#[test]
fn multi_allelic_warning() {
    let (seg, ps) = pileups!(
        [A G] Ref,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
        [T G] OB,
    );
    let metrics = to_metrics(&ps[0], &seg, None);
    let result = call(&metrics).unwrap();

    assert_denovo(result, 0.3);
}

#[test]
fn no_evidence() {
    let (seg, ps) = pileups!(
        [A G] Ref,
        [A G] OB,
    );
    let metrics = to_metrics(&ps[0], &seg, None);
    let result = call(&metrics).unwrap();

    assert_none(result);
}
