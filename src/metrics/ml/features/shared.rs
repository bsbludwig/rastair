use super::utils::one_hot_encode_base;
use crate::{
    metrics::{AlleleMetrics, MetricsForAlt, PileupMetrics},
    utils::IntoF64 as _,
};
use rastair_types::Base;

/// Common features shared across all ML models
///
/// This struct groups the feature sets that appear in the same order
/// across `cpg`, `denovo_cpg`, and `others` models.
pub struct CommonFeatures {
    /// One-hot encoding of reference and alt bases (8 features)
    pub base_encoding: [f64; 8],
    /// Position-level mapping quality metrics (2 features)
    pub position_metrics: [f64; 2],
    /// One-hot encoding of sequence context (16 features)
    pub context_encoding: [f64; 16],
    /// Regional entropy (1 feature)
    pub region_entropy: f64,
    /// Depth ratios for ref and alt on both strands (6 features)
    pub depth_ratios: [f64; 6],
    /// Base quality metrics for ref and alt (6 features)
    pub base_quality_metrics: [f64; 6],
    /// Mapping quality metrics for ref and alt (6 features)
    pub mapping_quality_metrics: [f64; 6],
    /// Read position and alignment metrics (6 features)
    pub read_metrics: [f64; 6],
}

impl CommonFeatures {
    /// Extract all common features from the current position
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
            region_entropy: pos.region_entropy,
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

/// Calculate `alt_score` for methylation-aware contexts (CpG and denovo CpG)
///
/// This uses strand-specific calculations:
/// - For C bases: uses OB (original bottom) strand
/// - For G bases: uses OT (original top) strand
pub fn alt_score_methylation_aware(
    alt: &AlleleMetrics,
    r: &AlleleMetrics,
    check_base: Base,
) -> f64 {
    use Base::*;

    if check_base == C {
        (alt.strand_count.ob.f() * alt.baseq_s.ob + 1.0).log2()
            - (r.strand_count.ob.f() * r.baseq_s.ob + 1.0).log2()
    } else {
        // For G or any other base, use OT strand
        (alt.strand_count.ot.f() * alt.baseq_s.ot + 1.0).log2()
            - (r.strand_count.ot.f() * r.baseq_s.ot + 1.0).log2()
    }
}

/// Calculate `alt_score` for non-methylation contexts
///
/// This uses combined strand data (total depth * average base quality)
pub fn alt_score_generic(alt: &AlleleMetrics, r: &AlleleMetrics) -> f64 {
    (alt.depth.f() * alt.baseq + 1.0).log2() - (r.depth.f() * r.baseq + 1.0).log2()
}
