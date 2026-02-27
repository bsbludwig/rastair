use super::super::shared::{CommonFeatures, alt_score_generic};
use crate::{
    metrics::{MetricsForAlt, PileupMetrics},
    utils::IntoF64 as _,
};
use color_eyre::Result;

pub const FEATURES: usize = 54;

pub fn others(
    current: &MetricsForAlt,
    _before: Option<&PileupMetrics>,
    _after: Option<&PileupMetrics>,
    buf: &mut [f64; FEATURES],
) -> Result<()> {
    let alt = current.alt;
    let PileupMetrics { ref_metrics: r, .. } = &current.metrics;

    let common = CommonFeatures::extract(current);
    let sb_alt = (alt.strand_count.ot + 1).f() / (alt.strand_count.ob + 1).f();
    let sb_ref = (r.strand_count.ot + 1).f() / (r.strand_count.ob + 1).f();
    let alt_score = alt_score_generic(alt, r);

    // Never change the order of these writes, as they are used in the model
    buf[0..8].copy_from_slice(&common.base_encoding);
    buf[8..10].copy_from_slice(&common.position_metrics);
    buf[10..26].copy_from_slice(&common.context_encoding);
    buf[26] = common.region_entropy;
    buf[27..33].copy_from_slice(&common.depth_ratios);
    buf[33..36].copy_from_slice(&[sb_alt, sb_ref, alt_score]);
    buf[36..42].copy_from_slice(&common.base_quality_metrics);
    buf[42..48].copy_from_slice(&common.mapping_quality_metrics);
    buf[48..54].copy_from_slice(&common.read_metrics);
    Ok(())
}
