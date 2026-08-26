use crate::{
    call::{denovo_cpg, methylation, variant_calling::VariantCallingParams},
    metrics::{Filters, PileupMetrics},
    utils::logging::ThisIsABug as _,
    vcf::RastairFilter,
};
use color_eyre::eyre::{ContextCompat as _, Result};
use tracing::instrument;

pub struct ThresholdFilterParams {
    pub variant_calling: VariantCallingParams,
    pub methylation: methylation::ThresholdParams,
    pub denovo_cpg: denovo_cpg::DenovoParams,
}

#[instrument(level = "debug", skip_all)]
pub fn apply_threshold_filters(
    current: &mut PileupMetrics,
    params: &ThresholdFilterParams,
) -> Result<()> {
    current.pos_filters.add(RastairFilter::LowDp, || {
        current.pos_metrics.depth < params.variant_calling.v_min_depth
    });

    for alt_base in current.alts() {
        let alt =
            current.alt_metrics(alt_base).wrap_err("Failed to get alt metrics").this_is_a_bug()?;

        let generic = {
            let mut filters = Filters::default();
            filters
                .add(RastairFilter::LowDp, || alt.alt.depth < params.variant_calling.v_min_depth);
            filters
        };
        let m_filters = methylation::add_filters(&params.methylation, &alt)?;
        let denovo_filters = denovo_cpg::add_filters(&params.denovo_cpg, &alt)?;

        let alt_filters = current
            .alt_filters_mut(alt_base)
            .wrap_err("Failed to get mutable alt metrics")
            .this_is_a_bug()?;

        alt_filters.filters.merge(generic);
        alt_filters.filters.merge(m_filters);
        alt_filters.filters.merge(denovo_filters);
    }
    Ok(())
}
