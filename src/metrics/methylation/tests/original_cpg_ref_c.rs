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
