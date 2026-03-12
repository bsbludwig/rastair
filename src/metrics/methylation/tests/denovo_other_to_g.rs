use super::*;

#[test]
fn standard_t_to_g() {
    let (seg, ps) = pileups!(
        [C T] Ref,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C A] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
    );
    let metrics = to_metrics(&ps[1], &seg, None);
    let result = call(&metrics).unwrap();

    assert_denovo(result, 0.4);
}

#[test]
fn multi_allelic_warning() {
    let (seg, ps) = pileups!(
        [C T] Ref,
        [C A] OB,
        [C A] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C G] OB,
        [C A] OT,
    );
    let metrics = to_metrics(&ps[1], &seg, None);
    let result = call(&metrics).unwrap();

    assert_denovo(result, 0.2);
}

#[test]
fn no_evidence() {
    let (seg, ps) = pileups!(
        [C T] Ref,
        [C T] OT,
    );
    let metrics = to_metrics(&ps[1], &seg, None);
    let result = call(&metrics).unwrap();

    assert_none(result);
}
