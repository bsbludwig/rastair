use super::super::shared::{CommonFeatures, alt_score_methylation_aware};
use super::super::utils::safe_div;
use crate::{
    metrics::{MetricsForAlt, PileupMetrics},
    utils::IntoF64 as _,
};
use color_eyre::{
    Result,
    eyre::{Context as _, ensure, eyre},
};
use rastair_types::Base;
use tracing::trace;

pub const FEATURES: usize = 55;

/// Generate features for CpG mutation candidates
///
/// Call this with `MetricsForAlt` where the alt is a CpG methylation candidate
pub fn cpg(
    current: &MetricsForAlt,
    before: Option<&PileupMetrics>,
    after: Option<&PileupMetrics>,
) -> Result<[f64; FEATURES]> {
    use Base::*;

    let alt = current.alt;
    let PileupMetrics { pileup, pos_metrics: pos, ref_metrics: r, .. } = &current.metrics;
    let ref_base = pileup.reference_base;

    ensure!(*pos.cpg, "cpg called on non-CpG position");
    ensure!(
        (ref_base == C && alt.base == T) || (ref_base == G && alt.base == A),
        "cpg called on non-methylation candidate"
    );

    let common = CommonFeatures::extract(current);
    let alt_score = alt_score_methylation_aware(alt, r, ref_base);

    let AdjecentFeatures { beta_ratio, alt_ad_adj, alt_score_adj } =
        calculate_adjacent_features(current, before, after)
            .wrap_err("Failed to calculate adjacent features for CpG")?;

    // Never change the order of these variables, as they are used in the model
    let mut features = Vec::with_capacity(FEATURES);
    features.extend_from_slice(&[alt_ad_adj, alt_score_adj]);
    features.extend_from_slice(&common.base_encoding);
    features.extend_from_slice(&common.position_metrics);
    features.extend_from_slice(&common.context_encoding);
    features.push(common.region_entropy);
    features.extend_from_slice(&common.depth_ratios);
    features.push(alt_score);
    features.extend_from_slice(&common.base_quality_metrics);
    features.extend_from_slice(&common.mapping_quality_metrics);
    features.extend_from_slice(&common.read_metrics);
    features.push(beta_ratio);

    features.try_into().map_err(|_: Vec<f64>| eyre!("Expected {FEATURES} CpG features"))
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
) -> Result<AdjecentFeatures> {
    use Base::*;

    let ref_base = c.metrics.ref_base();
    let res = if ref_base == C
        && let Some(after) = a
        && after.ref_base() == G
    {
        let c_alt = c.alt;
        ensure!(c_alt.base == T);
        let c_r = &c.metrics.ref_metrics;
        let r = &after.ref_metrics;

        let beta_center = safe_div(
            c_alt.strand_count.ot.f(),
            c_alt.strand_count.ot.f() + c_r.strand_count.ot.f(),
        );

        if let Some(MetricsForAlt { alt, .. }) = after.alt_metrics(A) {
            let alt_ad = alt.depth.f();
            let depth = after.pos_metrics.depth.f();
            let alt_ad_adj = safe_div(alt_ad, depth);

            // Calculate alt_score for G→A
            let alt_score_adj = (alt.strand_count.ot.f() * alt.baseq_s.ot + 1.0).log2()
                - (r.strand_count.ot.f() * c_r.baseq_s.ot + 1.0).log2();

            let beta_after =
                safe_div(alt.strand_count.ob.f(), (alt.strand_count.ob + r.strand_count.ob).f());

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
        && before.ref_base() == C
    {
        let c_alt = c.alt;
        ensure!(c_alt.base == A);
        let c_r = &c.metrics.ref_metrics;
        let r = &before.ref_metrics;

        let beta_center =
            safe_div(c_r.strand_count.ob.f(), c_alt.strand_count.ob.f() + c_r.strand_count.ob.f());

        if let Some(MetricsForAlt { alt, .. }) = before.alt_metrics(T) {
            let alt_ad = alt.depth.f();
            let depth = before.pos_metrics.depth.f();
            let alt_ad_adj = safe_div(alt_ad, depth);

            // Calculate alt_score for C→T
            let alt_score_adj = (alt.strand_count.ob.f() * alt.baseq_s.ob + 1.0).log2()
                - (r.strand_count.ob.f() * c_r.baseq_s.ob + 1.0).log2();

            let beta_before =
                safe_div(alt.strand_count.ot.f(), (alt.strand_count.ot + r.strand_count.ot).f());

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
    };

    Ok(res)
}
