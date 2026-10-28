use crate::{
    call::tests::pileup::Pileup,
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

        let pos = seqair_types::Pos0::try_from(pileup.pos).wrap_err("position out of range")?;
        let entropy = sliding_entropy.entropy_at(segment.pos_to_idx(pos)?);
        pileup.metrics_with_entropy(segment, entropy).wrap_err("Failed to calculate pileup metrics")
    })
}
