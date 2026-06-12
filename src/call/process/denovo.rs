use crate::metrics::{DenovoAdjecent, FormsDenovo, PileupMetrics};
use color_eyre::eyre::Result;
use seqair_types::{Base::*, Probability};
use tracing::instrument;

/// Set the `denovo_adj` field based on neighboring pileup metrics
#[instrument(level = "debug", skip_all)]
pub fn set_denovo_adj(
    before: Option<&PileupMetrics>,
    current: &mut PileupMetrics,
    after: Option<&PileupMetrics>,
) -> Result<()> {
    if let Some(before) = before
        && let Some(_denovo) =
            before.alts.iter().find(|alt| alt.metrics.denovo == FormsDenovo::ThisBecomesC)
        && current.ref_base() == G
    {
        current.pos_metrics.extended.denovo_adj = DenovoAdjecent::ThisIsTheMatchingG;
    } else if let Some(after) = after
        && let Some(_denovo) =
            after.alts.iter().find(|alt| alt.metrics.denovo == FormsDenovo::ThisBecomesG)
        && current.ref_base() == C
    {
        current.pos_metrics.extended.denovo_adj = DenovoAdjecent::ThisIsTheMatchingC;
    }
    Ok(())
}

/// If one position of a de-novo CpG passes, the entire CpG should pass
#[instrument(level = "debug", skip_all)]
pub fn propagate_denovo_pass_flags(
    before: Option<&PileupMetrics>,
    current: &mut PileupMetrics,
    after: Option<&PileupMetrics>,
    ml_threshold: Option<Probability>,
) -> Result<()> {
    if let Some(before) = before
        && let Some(denovo) =
            before.alts.iter().find(|alt| alt.metrics.denovo == FormsDenovo::ThisBecomesC)
        && current.pos_metrics.extended.denovo_adj == DenovoAdjecent::ThisIsTheMatchingG
        && denovo.filters.pass(ml_threshold)
    {
        current.pos_filters.other_pos_in_denovo_passes = true;
    } else if let Some(after) = after
        && let Some(denovo) =
            after.alts.iter().find(|alt| alt.metrics.denovo == FormsDenovo::ThisBecomesG)
        && current.pos_metrics.extended.denovo_adj == DenovoAdjecent::ThisIsTheMatchingC
        && denovo.filters.pass(ml_threshold)
    {
        current.pos_filters.other_pos_in_denovo_passes = true;
    }
    Ok(())
}
