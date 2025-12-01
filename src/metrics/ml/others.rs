use super::shared::{CommonFeatures, alt_score_generic};
use crate::{
    metrics::{MetricsForAlt, PileupMetrics},
    utils::IntoF64 as _,
};
use color_eyre::{Result, eyre::Context as _};
use ndarray::Array2;

pub fn others(
    current: &MetricsForAlt,
    _before: Option<&PileupMetrics>,
    _after: Option<&PileupMetrics>,
) -> Result<Array2<f64>> {
    let alt = current.alt;
    let PileupMetrics { ref_metrics: r, .. } = &current.metrics;

    let common = CommonFeatures::extract(current);

    // Calculate strand bias ratios
    let sb_alt = (alt.strand_count.ot + 1).f() / (alt.strand_count.ob + 1).f();
    let sb_ref = (r.strand_count.ot + 1).f() / (r.strand_count.ob + 1).f();

    let alt_score = alt_score_generic(alt, r);

    // Never change the order of these variables, as they are used in the model
    let mut features = Vec::with_capacity(54);
    features.extend_from_slice(&common.base_encoding);
    features.extend_from_slice(&common.position_metrics);
    features.extend_from_slice(&common.context_encoding);
    features.push(common.region_entropy);
    features.extend_from_slice(&common.depth_ratios);
    features.extend_from_slice(&[sb_alt, sb_ref, alt_score]);
    features.extend_from_slice(&common.base_quality_metrics);
    features.extend_from_slice(&common.mapping_quality_metrics);
    features.extend_from_slice(&common.read_metrics);

    Array2::from_shape_vec((1, features.len()), features)
        .wrap_err("Failed to create feature array for non-CpG variants")
}
