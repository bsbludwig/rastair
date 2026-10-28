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
        quality: Phred::from_phred(20),
        source: GenotypeSource::Binomial,
    };

    let metrics = to_metrics(&ps[0], &seg, Some(gt));
    let result = call(&metrics).unwrap();

    assert_denovo(result, 0.6);
}

/// A `1/2` of the de-novo `C` and an `A` leaves no chromosome with the
/// reference `T`, so every converted read is a methylated `C`: 6 of 10.
#[test]
fn compound_het_without_the_reference_is_not_confounded() {
    let (seg, ps) = pileups!(
        [T G] Ref,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [T G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
        [C G] OT,
        [A G] OB,
        [A G] OB,
        [A G] OB,
    );
    let gt = EstimatedGenotype {
        genotype: GenotypeTag::AltHet(NonZeroU8::new(1).unwrap(), NonZeroU8::new(2).unwrap()),
        quality: Phred::from_phred(20),
        source: GenotypeSource::Binomial,
    };

    let metrics = to_metrics(&ps[0], &seg, Some(gt));
    assert_eq!(metrics.alts.iter().map(|alt| alt.base).collect::<Vec<_>>(), [C, A]);
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
