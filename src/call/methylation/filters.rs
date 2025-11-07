use crate::{
    call::methylation::ThresholdParams,
    metrics::{Filters, MetricsForAlt},
    utils::IntoF64,
    vcf,
};
use color_eyre::Result;

pub fn add_filters(config: &ThresholdParams, current: &MetricsForAlt) -> Result<Filters> {
    let mut filters = Filters::default();

    if !current.is_evidence_for_methylation() {
        // Not a CpG site, skipping filters
        return Ok(filters);
    }

    let alt = &current.alt;
    let r = &current.metrics.pos_metrics;

    filters.add(vcf::lowDp, || alt.depth < config.m_min_depth);
    filters.add(vcf::m_vaf, || alt.allele_frequency.f() < config.m_vaf_min);
    filters.add(vcf::m_bq_ratio, || {
        let quality_ratio = (alt.depth.f() * alt.baseq.f() + 1.) / (r.depth.f() * r.baseq.f() + 1.);
        quality_ratio < config.m_bq_ratio_min
    });
    filters.add(vcf::m_pos, || {
        alt.position_in_read.f() < config.m_read_position_min
            || alt.position_in_read.f() > config.m_read_position_max
    });
    filters.add(vcf::m_highDp, || alt.depth > config.m_max_coverage);

    Ok(filters)
}
