use std::num::NonZeroU8;

use super::*;

#[test]
fn standard() {
    let (seg, ps) = pileups!(
        [C A] Ref,
        [C A] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
    );
    let metrics = to_metrics(&ps[1], &seg, None);
    let result = call(&metrics).unwrap();

    assert_denovo(result, 0.1);
}

#[test]
fn het_snp_adjustment() {
    let (seg, ps) = pileups!(
        [C A] Ref,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C G] OB,
        [C A] OT,
        [C A] OT,
    );
    let gt = EstimatedGenotype {
        genotype: GenotypeTag::RefHet(NonZeroU8::new(1).unwrap()),
        likelihood: Probability::new(0.8).unwrap(),
        confidence: Probability::new(0.99).unwrap(),
    };
    let metrics = to_metrics(&ps[1], &seg, Some(gt));
    let result = call(&metrics).unwrap();

    assert_denovo(result, 0.6);
}

#[test]
fn no_evidence() {
    let (seg, ps) = pileups!(
        [C A] Ref,
        [C A] OT,
    );
    let metrics = to_metrics(&ps[1], &seg, None);
    let result = call(&metrics).unwrap();

    assert_none(result);
}
