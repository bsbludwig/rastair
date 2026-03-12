use super::*;

#[test]
fn fully_methylated() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C G] OT,
    );
    let metrics = to_metrics(&ps[1], &seg, None);
    let result = call(&metrics).unwrap();

    assert_original(result, 1.0);
}

#[test]
fn unmethylated() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [C G] OB,
        [C G] OB,
        [C G] OB,
    );
    let metrics = to_metrics(&ps[1], &seg, None);
    let result = call(&metrics).unwrap();

    assert_original(result, 0.0);
}

#[test]
fn partially_methylated() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
    );
    let metrics = to_metrics(&ps[1], &seg, None);
    let result = call(&metrics).unwrap();

    assert_original(result, 0.7);
}

#[test]
fn het_snp() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C G] OB,
    );
    let metrics = to_metrics(&ps[1], &seg, Some(ct_genotype()));
    let result = call(&metrics).unwrap();

    assert_original(result, 0.71428571);
}

#[test]
fn het_snp_fifty_fifty() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
    );
    let metrics = to_metrics(&ps[1], &seg, Some(ct_genotype()));
    let result = call(&metrics).unwrap();

    assert_original(result, 0.0);
}

#[test]
fn no_evidence() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [C G] OT,
    );
    let metrics = to_metrics(&ps[1], &seg, None);
    let result = call(&metrics).unwrap();

    assert_none(result);
}

#[test]
fn no_alt_a() {
    let (seg, ps) = pileups!(
        [C G] Ref,
        [C G] OB,
        [C G] OB,
    );
    let metrics = to_metrics(&ps[1], &seg, None);
    let result = call(&metrics).unwrap();

    assert_original(result, 0.0);
}

#[test]
fn filtered_denovo_alt_does_not_change_beta() {
    let (seg, ps) = pileups!(
        [C G G] Ref,
        [C T G] OT,
        [C T G] OT,
        [C C G] OB,
    );
    let mut metrics = to_metrics(&ps[1], &seg, None);

    let denovo_alt = metrics.alts.iter_mut().find(|a| a.base == C).expect("expected C alt");
    denovo_alt.call = AltCall::ReadError;

    let result = call(&metrics).unwrap();
    assert_none(result);
}
