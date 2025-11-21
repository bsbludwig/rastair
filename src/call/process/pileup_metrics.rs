use crate::{
    call::{
        methylation::params::ThresholdParams, pileup::Pileup, variant_calling::VariantCallingParams,
    },
    metrics::{self, PileupMetrics},
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
        let span = tracing::trace_span!("pileup_metrics", pos=%pileup.pos);
        let _guard = span.enter();

        let mut current =
            PileupMetrics::new(pileup).wrap_err("Failed to calculate pileup metrics")?;

        current.pos_metrics.extended.region_entropy = segment
            .entropy_around::<100>(current.pileup.idx())
            .wrap_err("Failed to calculate region entropy")?;

        Ok(current)
    })
}
