use super::super::shared::{CommonFeatures, CommonSectionA, CommonSectionB, alt_score_generic};
use crate::metrics::ml::features::define_features;
use crate::{
    metrics::{MetricsForAlt, PileupMetrics},
    utils::IntoF32 as _,
};
use color_eyre::Result;

define_features! {
    /// ML features for "other" (non-methylation, non-denovo) SNV candidates.
    pub struct OthersFeatures {
        /// First block of [`CommonFeatures`].
        flatten section_a: CommonSectionA;
        /// Strand bias of the alt allele: (OT + 1) / (OB + 1).
        scalar sb_alt;
        /// Strand bias of the reference allele: (OT + 1) / (OB + 1).
        scalar sb_ref;
        /// Generic (strand-combined) alt score.
        scalar alt_score;
        /// Second block of [`CommonFeatures`].
        flatten section_b: CommonSectionB;
    }
}

impl OthersFeatures {
    pub fn extract(
        current: &MetricsForAlt,
        _before: Option<&PileupMetrics>,
        _after: Option<&PileupMetrics>,
    ) -> Result<OthersFeatures> {
        let alt = current.alt;
        let PileupMetrics { ref_metrics: r, .. } = &current.metrics;

        let common = CommonFeatures::extract(current);
        let sb_alt = (alt.strand_count.ot + 1).f() / (alt.strand_count.ob + 1).f();
        let sb_ref = (r.strand_count.ot + 1).f() / (r.strand_count.ob + 1).f();
        let alt_score = alt_score_generic(alt, r);

        Ok(OthersFeatures {
            section_a: CommonSectionA::from_common(&common),
            sb_alt,
            sb_ref,
            alt_score,
            section_b: CommonSectionB::from_common(&common),
        })
    }
}
