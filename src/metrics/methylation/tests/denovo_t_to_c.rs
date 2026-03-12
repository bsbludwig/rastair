use std::num::NonZeroU8;

use super::*;

#[test]
fn standard() {
    let (seg, ps) = pileups!(
        [T G] Ref,
        [T G] OT,
        [T G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
    );
    let metrics = to_metrics(&ps[0], &seg, None);
    let result = call(&metrics).unwrap();

    assert_denovo(result, 0.2);
}

#[test]
fn het_snp_adjustment() {
    let (seg, ps) = pileups!(
        [T G] Ref,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [C G] OT,
        [T G] OB,
        [T G] OB,
    );
    let gt = EstimatedGenotype {
        genotype: GenotypeTag::RefHet(NonZeroU8::new(1).unwrap()),
        likelihood: Probability::new(0.8).unwrap(),
        confidence: Probability::new(0.99).unwrap(),
    };

    let metrics = to_metrics(&ps[0], &seg, Some(gt));
    let result = call(&metrics).unwrap();

    assert_denovo(result, 0.6);
}

#[test]
fn no_evidence() {
    let (seg, ps) = pileups!(
        [T G] Ref,
        [T G] OB,
    );
    let metrics = to_metrics(&ps[0], &seg, None);
    let result = call(&metrics).unwrap();

    assert_none(result);
}
