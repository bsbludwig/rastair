use super::utils::one_hot_encode_base;
use crate::{
    metrics2::{MetricsForAlt, PileupMetrics},
    utils::conversion::IntoF64 as _,
    vcf::DeNovoCpGCandidate,
};
use ndarray::{Array2, array};
use rastair_types::Base::{self, *};
use tracing::trace;

/// Generate features for denovo CpG mutation candidates
///
/// Call this with `MetricsForAlt` where the alt is a denovo CpG candidate
pub fn denovo_cpg(
    current: &MetricsForAlt,
    before: Option<&PileupMetrics>,
    after: Option<&PileupMetrics>,
) -> Array2<f64> {
    let alt = current.alt;
    assert!(
        matches!(alt.denovo, DeNovoCpGCandidate::Candidate { .. }),
        "denovo_cpg called on non-denovo candidate"
    );

    let PileupMetrics { pileup, pos_metrics: pos, ref_metrics: r, .. } = &current.metrics;

    let ref_base = pileup.reference_base;
    let depth = pos.read_depth as f64;

    let seq_ctx = &pos.sequence_context;
    let (p1a, p1c, p1g, p1t) = one_hot_encode_base(seq_ctx.before_2);
    let (p2a, p2c, p2g, p2t) = one_hot_encode_base(seq_ctx.before_1);
    let (p4a, p4c, p4g, p4t) = one_hot_encode_base(seq_ctx.after_1);
    let (p5a, p5c, p5g, p5t) = one_hot_encode_base(seq_ctx.after_2);

    let (ref_a, ref_c, ref_g, ref_t) = one_hot_encode_base(ref_base);
    let (alt_a, alt_c, alt_g, alt_t) = one_hot_encode_base(alt.base);

    let alt_score = if alt.base == C {
        (alt.strand_count.ob.f() * alt.baseq_s.ob + 1.).log2()
            - (r.strand_count.ob.f() * r.baseq_s.ob + 1.).log2()
    } else if alt.base == A {
        (alt.strand_count.ot.f() * alt.baseq_s.ot + 1.).log2()
            - (r.strand_count.ot.f() * r.baseq_s.ot + 1.).log2()
    } else {
        unreachable!("denovo CpG alt base must be A or C")
    };

    let AdjecentFeatures { beta_ratio, alt_ad_adj, alt_score_adj, sb_adj } =
        calculate_adjacent_features(current, before, after);

    // Never change the order of these variables, as they are used in the model
    array![[
        alt_ad_adj,
        alt_score_adj,
        sb_adj,
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
    sb_adj: f64,
}

fn calculate_adjacent_features(
    // current
    c: &MetricsForAlt,
    // before
    b: Option<&PileupMetrics>,
    // after
    a: Option<&PileupMetrics>,
) -> AdjecentFeatures {
    let c_alt = c.alt;

    match c_alt.denovo {
        DeNovoCpGCandidate::Candidate { alt_base: Base::C, .. } => {
            let beta_center = {
                let c_count = c.metrics.alt(C).map(|x| x.strand_count.ot).unwrap_or_default();
                let t_count = c.metrics.alt(T).map(|x| x.strand_count.ot).unwrap_or_default();
                if c_count + t_count == 0 { 0.0 } else { t_count.f() / (t_count + c_count).f() }
            };

            if let Some(after) = a
                && let Some(MetricsForAlt { alt, .. }) = after.alt_metrics(A)
            {
                assert_eq!(after.pileup.reference_base, G, "De-novo CpG not followed by G");
                let r = &after.ref_metrics;

                let beta_after = {
                    let g_count = r.strand_count.ob;
                    let a_count = alt.strand_count.ob;
                    if g_count + a_count == 0 { 0.0 } else { a_count.f() / (a_count + g_count).f() }
                };
                let beta_ratio = (beta_center + 1.0).log2() - (beta_after + 1.0).log2();

                let alt_ad_adj = alt.depth.f() / after.pos_metrics.read_depth.f();
                let alt_score_adj = (alt.strand_count.ot.f() * alt.baseq_s.ot + 1.).log2()
                    - (r.strand_count.ot.f() * r.baseq_s.ot + 1.).log2();
                let sb_adj = (alt.strand_count.ob + 1).f() / (alt.strand_count.ot + 1).f();

                AdjecentFeatures { beta_ratio, alt_ad_adj, alt_score_adj, sb_adj }
            } else {
                let beta_ratio = (beta_center + 1.0).log2();
                AdjecentFeatures { beta_ratio, alt_ad_adj: 0.0, alt_score_adj: 0.0, sb_adj: 0.0 }
            }
        }
        DeNovoCpGCandidate::Candidate { alt_base: Base::G, .. } => {
            let beta_center = {
                let g_count = c.metrics.alt(G).map(|x| x.strand_count.ot).unwrap_or_default();
                let a_count = c.metrics.alt(A).map(|x| x.strand_count.ot).unwrap_or_default();
                if g_count + a_count == 0 { 0.0 } else { a_count.f() / (a_count + g_count).f() }
            };

            if let Some(before) = b
                && let Some(alt) = before.alt(T)
            {
                assert_eq!(before.pileup.reference_base, C, "De-novo CpG not preceded by C");
                let r = &before.ref_metrics;

                let beta_before = {
                    let c_count = r.strand_count.ob;
                    let t_count = alt.strand_count.ob;
                    if c_count + t_count == 0 { 0.0 } else { t_count.f() / (t_count + c_count).f() }
                };
                let beta_ratio = (beta_center + 1.0).log2() - (beta_before + 1.0).log2();

                let alt_ad_adj = alt.depth.f() / before.pos_metrics.read_depth.f();
                let alt_score_adj = (alt.strand_count.ob.f() * alt.baseq_s.ob + 1.).log2()
                    - (r.strand_count.ob.f() * r.baseq_s.ob + 1.).log2();
                let sb_adj = (alt.strand_count.ot + 1).f() / (alt.strand_count.ob + 1).f();

                AdjecentFeatures { beta_ratio, alt_ad_adj, alt_score_adj, sb_adj }
            } else {
                let beta_ratio = (beta_center + 1.0).log2();
                AdjecentFeatures { beta_ratio, alt_ad_adj: 0.0, alt_score_adj: 0.0, sb_adj: 0.0 }
            }
        }
        _ => {
            // Not a denovo CpG candidate or unexpected alt base
            trace!("No denovo CpG context found for adjacent feature calculation");
            AdjecentFeatures { beta_ratio: 0.0, alt_ad_adj: 0.0, alt_score_adj: 0.0, sb_adj: 0.0 }
        }
    }
}
