use color_eyre::eyre::Result;
use rastair_types::Probability;
use tracing::instrument;

use crate::{
    metrics::{DenovoAdjecent, FormsDenovo, PileupMetrics},
    vcf::InCpG,
};

#[instrument(level = "debug", skip_all)]
pub fn set_denovo_adj(
    before: Option<&PileupMetrics>,
    current: &mut PileupMetrics,
    after: Option<&PileupMetrics>,
) -> Result<()> {
    if let Some(before) = before
        && let Some(denovo) = before.alts.iter().filter_map(|alt| alt.metrics.denovo.some()).next()
        && denovo == FormsDenovo::ThisBecomesC
    {
        current.pos_metrics.extended.denovo_adj = DenovoAdjecent::ThisIsTheMatchingG;
    } else if let Some(after) = after
        && let Some(denovo) = after.alts.iter().filter_map(|alt| alt.metrics.denovo.some()).next()
        && denovo == FormsDenovo::ThisBecomesG
    {
        current.pos_metrics.extended.denovo_adj = DenovoAdjecent::ThisIsTheMatchingC;
    }
    Ok(())
}

/// If one position of a (de-novo) CpG passes, the entire CpG should pass
#[instrument(level = "info", skip_all)]
pub fn propagate_cpg_pass_flags(
    before: Option<&PileupMetrics>,
    current: &mut PileupMetrics,
    after: Option<&PileupMetrics>,
    ml_threshold: Option<Probability>,
) -> Result<()> {
    // We'll only change current, so if it already passes, skip it
    if current.pass(ml_threshold) {
        return Ok(());
    }

    // Check if the other position in the CpG passes
    if current.pos_metrics.cpg == InCpG::C
        && let Some(after) = after
        && after.pass(ml_threshold)
    {
        current.pos_filters.other_pos_in_cpg_passes = true;
    } else if current.pos_metrics.cpg == InCpG::G
        && let Some(before) = before
        && before.pass(ml_threshold)
    {
        current.pos_filters.other_pos_in_cpg_passes = true;
    }

    // Check if we're a de-novo CpG candidate and the other position passes
    for alt in &mut current.alts {
        if alt.metrics.denovo == FormsDenovo::ThisBecomesC
            && let Some(after) = after
            && after.pass(ml_threshold)
        {
            alt.filters.filters.other_pos_in_cpg_passes = true;
        } else if alt.metrics.denovo == FormsDenovo::ThisBecomesG
            && let Some(before) = before
            && before.pass(ml_threshold)
        {
            alt.filters.filters.other_pos_in_cpg_passes = true;
        }
    }

    Ok(())
}
