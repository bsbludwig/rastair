use crate::{
    call::pileup::Pileup,
    metrics::{PileupMetrics, entropy::SlidingEntropy},
    sequence::Segment,
};
use color_eyre::eyre::{Result, WrapErr};
use tracing::instrument;

#[instrument(level = "info", skip_all)]
pub fn calculate_pileup_metrics(
    pileups: impl Iterator<Item = Pileup>,
    segment: &Segment,
) -> impl Iterator<Item = Result<PileupMetrics>> {
    let mut sliding_entropy = SlidingEntropy::new(segment);

    pileups.into_iter().map(move |pileup| {
        let span = tracing::trace_span!("pileup_metrics", pos=%pileup.pos);
        let _guard = span.enter();

        let mut current =
            PileupMetrics::new(pileup).wrap_err("Failed to calculate pileup metrics")?;

        current.pos_metrics.extended.region_entropy =
            sliding_entropy.entropy_at(current.pileup.idx());

        Ok(current)
    })
}
