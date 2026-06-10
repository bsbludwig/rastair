use super::super::shared::{
    CommonFeatures, CommonSectionA, CommonSectionB, alt_score_methylation_aware,
};
use super::super::utils::safe_div;
use crate::metrics::ml::features::define_features;
use crate::{
    metrics::{MetricsForAlt, PileupMetrics},
    utils::IntoF64 as _,
};
use color_eyre::{
    Result,
    eyre::{Context as _, ensure},
};
use rastair_types::Base;
use tracing::trace;

define_features! {
    /// ML features for a CpG methylation candidate (ref C→T on OT, or ref G→A on OB).
    pub struct CpgFeatures {
        /// Adjacent-position alt allele depth fraction.
        scalar alt_ad_adj;
        /// Adjacent-position methylation-aware alt score.
        scalar alt_score_adj;
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

impl CpgFeatures {
    /// Build features for a CpG mutation candidate.
    ///
    /// `current` must be a `MetricsForAlt` whose alt is a CpG methylation candidate.
    pub fn extract(
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<CpgFeatures> {
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

        Ok(CpgFeatures {
            alt_ad_adj,
            alt_score_adj,
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
