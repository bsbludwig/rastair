use super::utils::one_hot_encode_base;
use crate::{
    metrics2::{MetricsForAlt, PileupMetrics},
    utils::conversion::IntoF64 as _,
};
use ndarray::{Array2, array};
use rastair_types::Base;
use tracing::trace;

/// Generate features for CpG mutation candidates
///
/// Call this with `MetricsForAlt` where the alt is a CpG methylation candidate
pub fn cpg(
    current: &MetricsForAlt,
    before: Option<&PileupMetrics>,
    after: Option<&PileupMetrics>,
) -> Array2<f64> {
    use Base::*;

    let alt = current.alt;
    let PileupMetrics { pileup, pos_metrics: pos, ref_metrics: r, .. } = &current.metrics;
    let ref_base = pileup.reference_base;

    assert!(pileup.is_cpg, "cpg called on non-CpG position");
    assert!(
        (ref_base == C && alt.base == T) || (ref_base == G && alt.base == A),
        "cpg called on non-methylation candidate"
    );

    let depth = pos.read_depth as f64;

    let seq_ctx = &pos.sequence_context;
    let (p1a, p1c, p1g, p1t) = one_hot_encode_base(seq_ctx.before_2);
    let (p2a, p2c, p2g, p2t) = one_hot_encode_base(seq_ctx.before_1);
    let (p4a, p4c, p4g, p4t) = one_hot_encode_base(seq_ctx.after_1);
    let (p5a, p5c, p5g, p5t) = one_hot_encode_base(seq_ctx.after_2);

    let (ref_a, ref_c, ref_g, ref_t) = one_hot_encode_base(ref_base);
    let (alt_a, alt_c, alt_g, alt_t) = one_hot_encode_base(alt.base);

    let alt_score = if ref_base == C {
        // For C: use "ob" (original bottom) strand data
        (alt.strand_count.ob.f() * alt.baseq_s.ob + 1.0).log2()
            - (r.strand_count.ob.f() * r.baseq_s.ob + 1.0).log2()
    } else {
        // For G: use "ot" (original top) strand data
        (alt.strand_count.ot.f() * alt.baseq_s.ot + 1.0).log2()
            - (r.strand_count.ot.f() * r.baseq_s.ot + 1.0).log2()
    };

    let AdjecentFeatures { beta_ratio, alt_ad_adj, alt_score_adj } =
        calculate_adjacent_features(current, before, after);

    // Never change the order of these variables, as they are used in the model
    array![[
        alt_ad_adj,
        alt_score_adj,
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
        beta_ratio
    ]]
}

struct AdjecentFeatures {
    beta_ratio: f64,
    alt_ad_adj: f64,
    alt_score_adj: f64,
}

fn calculate_adjacent_features(
    // current
    c: &MetricsForAlt,
    // before
    b: Option<&PileupMetrics>,
    // after
    a: Option<&PileupMetrics>,
) -> AdjecentFeatures {
    use Base::*;

    let ref_base = c.metrics.pileup.reference_base;
    if ref_base == C
        && let Some(after) = a
        && after.pileup.reference_base == G
    {
        let c_alt = c.alt;
        assert_eq!(c_alt.base, T);
        let c_r = &c.metrics.ref_metrics;
        let r = &after.ref_metrics;

        let beta_center =
            c_alt.strand_count.ot.f() / (c_alt.strand_count.ot.f() + c_r.strand_count.ot.f());

        if let Some(MetricsForAlt { alt, .. }) = after.alt_metrics(A) {
            let alt_ad = alt.depth.f();
            let depth = after.pos_metrics.read_depth.f();
            let alt_ad_adj = alt_ad / depth;

            // Calculate alt_score for G→A
            let alt_score_adj = (alt.strand_count.ot.f() * alt.baseq_s.ot + 1.0).log2()
                - (r.strand_count.ot.f() * c_r.baseq_s.ot + 1.0).log2();

            let beta_after =
                alt.strand_count.ob.f() / (alt.strand_count.ob + r.strand_count.ob).f();

            let beta_ratio = (beta_center + 1.0).log2() - (beta_after + 1.0).log2();

            AdjecentFeatures { beta_ratio, alt_ad_adj, alt_score_adj }
        } else {
            AdjecentFeatures {
                beta_ratio: (beta_center + 1.0).log2(),
                alt_ad_adj: 0.,
                alt_score_adj: 0.,
            }
        }
    } else if ref_base == G
        && let Some(before) = b
        && before.pileup.reference_base == C
    {
        let c_alt = c.alt;
        assert_eq!(c_alt.base, A);
        let c_r = &c.metrics.ref_metrics;
        let r = &before.ref_metrics;

        let beta_center =
            c_r.strand_count.ob.f() / (c_alt.strand_count.ob.f() + c_r.strand_count.ob.f());

        if let Some(MetricsForAlt { alt, .. }) = before.alt_metrics(T) {
            let alt_ad = alt.depth.f();
            let depth = before.pos_metrics.read_depth.f();
            let alt_ad_adj = alt_ad / depth;

            // Calculate alt_score for C→T
            let alt_score_adj = (alt.strand_count.ob.f() * alt.baseq_s.ob + 1.0).log2()
                - (r.strand_count.ob.f() * c_r.baseq_s.ob + 1.0).log2();

            let beta_before =
                alt.strand_count.ot.f() / (alt.strand_count.ot + r.strand_count.ot).f();

            let beta_ratio = (beta_center + 1.0).log2() - (beta_before + 1.0).log2();

            AdjecentFeatures { beta_ratio, alt_ad_adj, alt_score_adj }
        } else {
            AdjecentFeatures {
                beta_ratio: (beta_center + 1.0).log2(),
                alt_ad_adj: 0.,
                alt_score_adj: 0.,
            }
        }
    } else {
        trace!(%ref_base, "No adjacent evidence for methylation");
        AdjecentFeatures { beta_ratio: 0., alt_ad_adj: 0., alt_score_adj: 0. }
    }
}
