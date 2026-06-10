use super::super::shared::{
    CommonFeatures, CommonSectionA, CommonSectionB, alt_score_methylation_aware,
};
use super::super::utils::safe_div;
use crate::metrics::ml::features::define_features;
use crate::{
    metrics::{FormsDenovo, MetricsForAlt, PileupMetrics},
    utils::IntoF64 as _,
};
use color_eyre::{
    Result,
    eyre::{Context as _, bail, ensure},
};
use rastair_types::Base::*;
use tracing::trace;

define_features! {
    /// ML features for a de-novo CpG candidate (an X→C or X→G SNV that creates a CpG).
    pub struct DenovoCpgFeatures {
        /// Adjacent-position alt allele depth fraction.
        scalar alt_ad_adj;
        /// Adjacent-position methylation-aware alt score.
        scalar alt_score_adj;
        /// Adjacent-position strand bias.
        scalar sb_adj;
        /// First block of [`CommonFeatures`].
        flatten section_a: CommonSectionA;
        /// Methylation-aware alt score at this position.
        scalar alt_score;
        /// Second block of [`CommonFeatures`].
        flatten section_b: CommonSectionB;
        /// Log-ratio of this position's beta to the adjacent CpG partner's beta.
        scalar beta_ratio;
    }
}

impl DenovoCpgFeatures {
    /// Build features for a de-novo CpG mutation candidate.
    ///
    /// `current` must be a `MetricsForAlt` whose alt is a de-novo CpG candidate.
    pub fn extract(
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<DenovoCpgFeatures> {
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

        Ok(DenovoCpgFeatures {
            alt_ad_adj,
            alt_score_adj,
            sb_adj,
            section_a: CommonSectionA::from_common(&common),
            alt_score,
            section_b: CommonSectionB::from_common(&common),
            beta_ratio,
        })
    }
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
