use color_eyre::eyre::Result;
use rastair_types::Probability;
use tracing::instrument;

use crate::{
    metrics::{DenovoAdjecent, FormsDenovo, PileupMetrics},
    utils::{Surrounding, surrounding_pileups},
    vcf::InCpG,
};

#[instrument(level = "info", skip_all)]
pub fn set_denovo_adj(pileups: &mut [PileupMetrics]) {
    let pileups_len = pileups.len();
    for i in 0..pileups_len {
        let Surrounding { before, current, after } = surrounding_pileups(pileups, i);
        if let Some(before) = before
            && let Some(denovo) =
                before.alts.iter().filter_map(|alt| alt.metrics.denovo.some()).next()
            && denovo == FormsDenovo::ThisBecomesC
        {
            current.pos_metrics.extended.denovo_adj = DenovoAdjecent::ThisIsTheMatchingG;
        } else if let Some(after) = after
            && let Some(denovo) =
                after.alts.iter().filter_map(|alt| alt.metrics.denovo.some()).next()
            && denovo == FormsDenovo::ThisBecomesG
        {
            current.pos_metrics.extended.denovo_adj = DenovoAdjecent::ThisIsTheMatchingC;
        }
    }
}

/// If one position of a (de-novo) CpG passes, the entire CpG should pass
#[instrument(level = "info", skip_all)]
pub fn propagate_cpg_pass_flags(
    pileups: &mut [PileupMetrics],
    ml_threshold: Option<Probability>,
) -> Result<()> {
    let pileups_len = pileups.len();
    for i in 0..pileups_len {
        let Surrounding { before, current, after } = surrounding_pileups(pileups, i);

        // We'll only change current, so if it already passes, skip it
        if current.pass(ml_threshold) {
            continue;
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
    }

    Ok(())
}
