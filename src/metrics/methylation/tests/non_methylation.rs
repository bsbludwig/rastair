use super::*;

#[test]
fn non_cpg_context() {
    let (seg, ps) = pileups!(
        [A T] Ref,
        [A T] OT,
    );
    let metrics = to_metrics(&ps[0], &seg, None);
    let result = call(&metrics).unwrap();

    assert_none(result);
}

#[test]
fn wrong_context() {
    let (seg, ps) = pileups!(
        [C A] Ref,
        [C A] OT,
    );
    let metrics = to_metrics(&ps[0], &seg, None);
    let result = call(&metrics).unwrap();

    assert_none(result);
}
