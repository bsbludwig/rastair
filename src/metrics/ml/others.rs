use super::utils::one_hot_encode_base;
use crate::{
    metrics::{MetricsForAlt, PileupMetrics},
    utils::IntoF64 as _,
};
use ndarray::{Array2, array};

pub fn others(
    current: &MetricsForAlt,
    _before: Option<&PileupMetrics>,
    _after: Option<&PileupMetrics>,
) -> Array2<f64> {
    let alt = current.alt;

    let PileupMetrics { pileup, pos_metrics: pos, ref_metrics: r, .. } = &current.metrics;

    let ref_base = pileup.reference_base;
    let depth = pos.read_depth as f64;

    let seq_ctx = &pileup.context;
    let (p1a, p1c, p1g, p1t) = one_hot_encode_base(seq_ctx.before_2);
    let (p2a, p2c, p2g, p2t) = one_hot_encode_base(seq_ctx.before_1);
    let (p4a, p4c, p4g, p4t) = one_hot_encode_base(seq_ctx.after_1);
    let (p5a, p5c, p5g, p5t) = one_hot_encode_base(seq_ctx.after_2);

    let (ref_a, ref_c, ref_g, ref_t) = one_hot_encode_base(ref_base);
    let (alt_a, alt_c, alt_g, alt_t) = one_hot_encode_base(alt.base);

    // Calculate strand bias ratios
    let sb_alt = (alt.strand_count.ot + 1).f() / (alt.strand_count.ob + 1).f();
    let sb_ref = (r.strand_count.ot + 1).f() / (r.strand_count.ob + 1).f();

    // Calculate alt_score
    let alt_score = (alt.depth.f() * alt.baseq + 1.).log2() - (r.depth.f() * r.baseq + 1.).log2();

    // Never change the order of these variables, as they are used in the model
    array![[
        ref_a,
        ref_c,
        ref_g,
        ref_t,
        alt_a,
        alt_c,
        alt_g,
        alt_t,
        pos.mapq.f(),
        pos.mapq0.f(),
        p1a,
        p1c,
        p1g,
        p1t,
        p2a,
        p2c,
        p2g,
        p2t,
        p4a,
        p4c,
        p4g,
        p4t,
        p5a,
        p5c,
        p5g,
        p5t,
        pos.region_entropy,
        r.depth.f() / depth,
        alt.depth.f() / depth,
        r.strand_count.ot.f() / depth,
        r.strand_count.ob.f() / depth,
        alt.strand_count.ot.f() / depth,
        alt.strand_count.ob.f() / depth,
        sb_alt,
        sb_ref,
        alt_score,
        r.baseq.f(),
        alt.baseq.f(),
        r.baseq_s.ot.f(),
        r.baseq_s.ob.f(),
        alt.baseq_s.ot.f(),
        alt.baseq_s.ob.f(),
        r.mapq_s.ot.f(),
        r.mapq_s.ob.f(),
        alt.mapq_s.ot.f(),
        alt.mapq_s.ob.f(),
        r.mapq.f(),
        alt.mapq.f(),
        r.position_in_read.f(),
        alt.position_in_read.f(),
        r.num_aligned_bases.f(),
        alt.num_aligned_bases.f(),
        r.num_indels.f(),
        alt.num_indels.f(),
    ]]
}
