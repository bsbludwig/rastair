use super::define_features;
use super::utils::one_hot_encode_base;
use crate::{
    metrics::{AlleleMetrics, MetricsForAlt, PileupMetrics},
    utils::IntoF32 as _,
};
use seqair_types::Base;

/// Common features shared across all ML models.
///
/// The fields appear in two sections separated by per-model features:
/// - Section A (33): `base_encoding` + `position_metrics` + `context_encoding` + `region_entropy` + `depth_ratios`
/// - Section B (18): `base_quality_metrics` + `mapping_quality_metrics` + `read_metrics`
pub struct CommonFeatures {
    /// One-hot encoding of reference and alt bases (8 features)
    pub base_encoding: [f32; 8],
    /// Position-level mapping quality metrics (2 features)
    pub position_metrics: [f32; 2],
    /// One-hot encoding of sequence context (16 features)
    pub context_encoding: [f32; 16],
    /// Regional entropy (1 feature)
    pub region_entropy: f32,
    /// Depth ratios for ref and alt on both strands (6 features)
    pub depth_ratios: [f32; 6],
    /// Base quality metrics for ref and alt (6 features)
    pub base_quality_metrics: [f32; 6],
    /// Mapping quality metrics for ref and alt (6 features)
    pub mapping_quality_metrics: [f32; 6],
    /// Read position and alignment metrics (6 features)
    pub read_metrics: [f32; 6],
    // FIXME: consider adding `homopolymer_run` and `soft_clip_rate`
}

impl CommonFeatures {
    pub fn extract(current: &MetricsForAlt) -> Self {
        let alt = current.alt;
        let PileupMetrics { pileup, pos_metrics: pos, ref_metrics: r, .. } = &current.metrics;

        let ref_base = pileup.reference_base;
        let depth = pos.depth.f();
        // Prevent division by zero resulting in NaN values
        let depth = if depth > 0.0 { depth } else { 1.0 };

        let seq_ctx = &pileup.context;
        let (p1a, p1c, p1g, p1t) = one_hot_encode_base(seq_ctx.before_2);
        let (p2a, p2c, p2g, p2t) = one_hot_encode_base(seq_ctx.before_1);
        let (p4a, p4c, p4g, p4t) = one_hot_encode_base(seq_ctx.after_1);
        let (p5a, p5c, p5g, p5t) = one_hot_encode_base(seq_ctx.after_2);

        let (ref_a, ref_c, ref_g, ref_t) = one_hot_encode_base(ref_base);
        let (alt_a, alt_c, alt_g, alt_t) = one_hot_encode_base(alt.base);

        Self {
            base_encoding: [ref_a, ref_c, ref_g, ref_t, alt_a, alt_c, alt_g, alt_t],
            position_metrics: [pos.mapq.f(), pos.mapq0.f()],
            context_encoding: [
                p1a, p1c, p1g, p1t, p2a, p2c, p2g, p2t, p4a, p4c, p4g, p4t, p5a, p5c, p5g, p5t,
            ],
            region_entropy: pos.region_entropy.f(),
            depth_ratios: [
                r.depth.f() / depth,
                alt.depth.f() / depth,
                r.strand_count.ot.f() / depth,
                r.strand_count.ob.f() / depth,
                alt.strand_count.ot.f() / depth,
                alt.strand_count.ob.f() / depth,
            ],
            base_quality_metrics: [
                r.baseq.f(),
                alt.baseq.f(),
                r.baseq_s.ot.f(),
                r.baseq_s.ob.f(),
                alt.baseq_s.ot.f(),
                alt.baseq_s.ob.f(),
            ],
            mapping_quality_metrics: [
                r.mapq_s.ot.f(),
                r.mapq_s.ob.f(),
                alt.mapq_s.ot.f(),
                alt.mapq_s.ob.f(),
                r.mapq.f(),
                alt.mapq.f(),
            ],
            read_metrics: [
                r.position_in_read.f(),
                alt.position_in_read.f(),
                r.num_aligned_bases.f(),
                alt.num_aligned_bases.f(),
                r.num_indels.f(),
                alt.num_indels.f(),
            ],
        }
    }
}

define_features! {
    /// First contiguous block of [`CommonFeatures`] (33 features).
    ///
    /// The alt-based models interleave a model-specific scalar (`alt_score`)
    /// between this section and [`CommonSectionB`], which is why the common
    /// features are split into two `#[repr(C)]` pieces instead of one.
    pub struct CommonSectionA {
        array base_encoding: 8 = [
            "ref_A", "ref_C", "ref_G", "ref_T", "alt_A", "alt_C", "alt_G", "alt_T",
        ];
        array position_metrics: 2 = ["pos_mapq", "pos_mapq0"];
        array context_encoding: 16 = [
            "ctx_before_2_A", "ctx_before_2_C", "ctx_before_2_G", "ctx_before_2_T",
            "ctx_before_1_A", "ctx_before_1_C", "ctx_before_1_G", "ctx_before_1_T",
            "ctx_after_1_A", "ctx_after_1_C", "ctx_after_1_G", "ctx_after_1_T",
            "ctx_after_2_A", "ctx_after_2_C", "ctx_after_2_G", "ctx_after_2_T",
        ];
        scalar region_entropy;
        array depth_ratios: 6 = [
            "depth_ratio_ref", "depth_ratio_alt",
            "depth_ratio_ref_ot", "depth_ratio_ref_ob",
            "depth_ratio_alt_ot", "depth_ratio_alt_ob",
        ];
    }
}

define_features! {
    /// Second contiguous block of [`CommonFeatures`] (18 features).
    pub struct CommonSectionB {
        array base_quality_metrics: 6 = [
            "baseq_ref", "baseq_alt",
            "baseq_ref_ot", "baseq_ref_ob", "baseq_alt_ot", "baseq_alt_ob",
        ];
        array mapping_quality_metrics: 6 = [
            "mapq_ref_ot", "mapq_ref_ob", "mapq_alt_ot", "mapq_alt_ob",
            "mapq_ref", "mapq_alt",
        ];
        array read_metrics: 6 = [
            "pos_in_read_ref", "pos_in_read_alt",
            "num_aligned_ref", "num_aligned_alt",
            "num_indels_ref", "num_indels_alt",
        ];
    }
}

impl CommonSectionA {
    pub fn from_common(c: &CommonFeatures) -> Self {
        Self {
            base_encoding: c.base_encoding,
            position_metrics: c.position_metrics,
            context_encoding: c.context_encoding,
            region_entropy: c.region_entropy,
            depth_ratios: c.depth_ratios,
        }
    }
}

impl CommonSectionB {
    pub fn from_common(c: &CommonFeatures) -> Self {
        Self {
            base_quality_metrics: c.base_quality_metrics,
            mapping_quality_metrics: c.mapping_quality_metrics,
            read_metrics: c.read_metrics,
        }
    }
}

/// Calculate `alt_score` for methylation-aware contexts (CpG and denovo CpG)
///
/// This uses strand-specific calculations:
/// - For C bases: uses OB (original bottom) strand
/// - For G bases: uses OT (original top) strand
pub fn alt_score_methylation_aware(
    alt: &AlleleMetrics,
    r: &AlleleMetrics,
    check_base: Base,
) -> f32 {
    use Base::*;

    if check_base == C {
        (alt.strand_count.ob.f() * alt.baseq_s.ob.f() + 1.0).log2()
            - (r.strand_count.ob.f() * r.baseq_s.ob.f() + 1.0).log2()
    } else {
        // For G or any other base, use OT strand
        (alt.strand_count.ot.f() * alt.baseq_s.ot.f() + 1.0).log2()
            - (r.strand_count.ot.f() * r.baseq_s.ot.f() + 1.0).log2()
    }
}

/// Calculate `alt_score` for non-methylation contexts
///
/// This uses combined strand data (total depth * average base quality)
pub fn alt_score_generic(alt: &AlleleMetrics, r: &AlleleMetrics) -> f32 {
    (alt.depth.f() * alt.baseq.f() + 1.0).log2() - (r.depth.f() * r.baseq.f() + 1.0).log2()
}
