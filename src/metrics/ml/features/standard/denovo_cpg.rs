use super::super::shared::{CommonFeatures, alt_score_methylation_aware};
use super::super::utils::safe_div;
use crate::{
    metrics::{FormsDenovo, MetricsForAlt, PileupMetrics},
    utils::IntoF64 as _,
};
use color_eyre::{
    Result,
    eyre::{Context as _, bail, ensure, eyre},
};
use rastair_types::Base::*;
use tracing::trace;

pub const FEATURES: usize = 56;

/// Generate features for denovo CpG mutation candidates
///
/// Call this with `MetricsForAlt` where the alt is a denovo CpG candidate
pub fn denovo_cpg(
    current: &MetricsForAlt,
    before: Option<&PileupMetrics>,
    after: Option<&PileupMetrics>,
) -> Result<[f64; FEATURES]> {
    let alt = current.alt;
    ensure!(*alt.denovo, "denovo_cpg called on non-denovo candidate");

    let PileupMetrics { ref_metrics: r, .. } = &current.metrics;

    if alt.base != C && alt.base != G {
        bail!("denovo CpG alt base must be C or G")
    }

    let common = CommonFeatures::extract(current);
    let alt_score = alt_score_methylation_aware(alt, r, alt.base);

    let AdjecentFeatures { beta_ratio, alt_ad_adj, alt_score_adj, sb_adj } =
        calculate_adjacent_features(current, before, after)
            .wrap_err("Failed to calculate adjacent features for de-novo CpG")?;

    // Never change the order of these variables, as they are used in the model
    let mut features = Vec::with_capacity(FEATURES);
    features.extend_from_slice(&[alt_ad_adj, alt_score_adj, sb_adj]);
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

    features.try_into().map_err(|_: Vec<f64>| eyre!("Expected {FEATURES} denovo CpG features"))
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
) -> Result<AdjecentFeatures> {
    let c_alt = c.alt;

    Ok(match c_alt.denovo {
        FormsDenovo::ThisBecomesC => {
            let beta_center = {
                let c_count = c.metrics.alt(C).map(|x| x.strand_count.ot).unwrap_or_default();
                let t_count = c.metrics.alt(T).map(|x| x.strand_count.ot).unwrap_or_default();
                safe_div(t_count.f(), (t_count + c_count).f())
            };

            if let Some(after) = a
                && let Some(MetricsForAlt { alt, .. }) = after.alt_metrics(A)
            {
                ensure!(after.ref_base() == G, "De-novo CpG not followed by G");
                let r = &after.ref_metrics;

                let beta_after = {
                    let g_count = r.strand_count.ob;
                    let a_count = alt.strand_count.ob;
                    safe_div(a_count.f(), (a_count + g_count).f())
                };
                let beta_ratio = (beta_center + 1.0).log2() - (beta_after + 1.0).log2();

                let alt_ad_adj = safe_div(alt.depth.f(), after.pos_metrics.depth.f());
                let alt_score_adj = (alt.strand_count.ot.f() * alt.baseq_s.ot + 1.).log2()
                    - (r.strand_count.ot.f() * r.baseq_s.ot + 1.).log2();
                let sb_adj = (alt.strand_count.ob + 1).f() / (alt.strand_count.ot + 1).f();

                AdjecentFeatures { beta_ratio, alt_ad_adj, alt_score_adj, sb_adj }
            } else {
                let beta_ratio = (beta_center + 1.0).log2();
                AdjecentFeatures { beta_ratio, alt_ad_adj: 0.0, alt_score_adj: 0.0, sb_adj: 0.0 }
            }
        }
        FormsDenovo::ThisBecomesG => {
            let beta_center = {
                let g_count = c.metrics.alt(G).map(|x| x.strand_count.ot).unwrap_or_default();
                let a_count = c.metrics.alt(A).map(|x| x.strand_count.ot).unwrap_or_default();
                safe_div(a_count.f(), (a_count + g_count).f())
            };

            if let Some(before) = b
                && let Some(alt) = before.alt(T)
            {
                ensure!(before.ref_base() == C, "De-novo CpG not preceded by C");
                let r = &before.ref_metrics;

                let beta_before = {
                    let c_count = r.strand_count.ob;
                    let t_count = alt.strand_count.ob;
                    safe_div(t_count.f(), (t_count + c_count).f())
                };
                let beta_ratio = (beta_center + 1.0).log2() - (beta_before + 1.0).log2();

                let alt_ad_adj = safe_div(alt.depth.f(), before.pos_metrics.depth.f());
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
    })
}
