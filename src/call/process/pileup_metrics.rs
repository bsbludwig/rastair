use crate::{
    call::{
        methylation::params::MethylationCallingParams, pileup::Pileup,
        variant_calling::VariantCallingParams,
    },
    metrics::{self, PileupMetrics, PositionMetricsExt},
    sequence::Segment,
};
use color_eyre::eyre::{Result, WrapErr};
use tracing::instrument;

pub struct PileupMetricsParams {
    pub variant_calling: VariantCallingParams,
    pub methylation: MethylationCallingParams,
}

#[instrument(level = "info", skip_all)]
pub fn calculate_pileup_metrics(
    pileups: impl Iterator<Item = Pileup>,
    segment: &Segment,
    params: &PileupMetricsParams,
) -> impl Iterator<Item = Result<PileupMetrics>> {
    pileups.into_iter().map(PileupMetrics::new).map(move |metrics| {
        // Set "extended" metrics that depend on the pileup and some external params
        let mut current = metrics?;

        let genotype = current.pileup.estimate_genotype(params.variant_calling.error_model);
        let methylated = metrics::methylation::call(&params.methylation.thresholds, &current)?
            .unwrap_or_default();

        let region_entropy = segment
            .entropy_around::<100>(current.pileup.idx())
            .wrap_err("Failed to calculate region entropy")?;

        let ext = PositionMetricsExt {
            genotype,
            methylated,
            region_entropy,
            denovo_adj: metrics::DenovoAdjecent::No,
        };
        current.set_extended_metrics(ext);

        Ok(current)
    })
}
