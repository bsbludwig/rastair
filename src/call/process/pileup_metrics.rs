use crate::{
    call::{
        methylation::params::ThresholdParams, pileup::Pileup, variant_calling::VariantCallingParams,
    },
    metrics::{self, PileupMetrics, PositionMetricsExt},
    sequence::Segment,
};
use color_eyre::eyre::{Result, WrapErr};
use tracing::instrument;

pub struct PileupMetricsParams {
    pub variant_calling: VariantCallingParams,
    pub methylation: ThresholdParams,
}

#[instrument(level = "info", skip_all)]
pub fn calculate_pileup_metrics(
    pileups: impl Iterator<Item = Pileup>,
    segment: &Segment,
    params: &PileupMetricsParams,
) -> impl Iterator<Item = Result<PileupMetrics>> {
    pileups.into_iter().map(move |pileup| {
        let mut current =
            PileupMetrics::new(pileup).wrap_err("Failed to calculate pileup metrics")?;

        // Set "extended" metrics that depend on the segment and params. This is
        // done in a separate step since it also uses the pileup we just
        // constructed.
        current.set_extended_metrics({
            let genotype = current.pileup.estimate_genotype(params.variant_calling.error_model);
            let methylated =
                metrics::methylation::call(&params.methylation, &current)?.unwrap_or_default();

            let region_entropy = segment
                .entropy_around::<100>(current.pileup.idx())
                .wrap_err("Failed to calculate region entropy")?;

            PositionMetricsExt {
                genotype,
                methylated,
                region_entropy,
                denovo_adj: metrics::DenovoAdjecent::No,
            }
        });

        Ok(current)
    })
}
