use crate::metrics::{AltCall, DenovoAdjecent, FormsDenovo, PileupMetrics, alt_forms_denovo};
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
    // Assigned, not only raised: a re-run after a neighbour's verdict changed must be able to clear it.
    current.pos_filters.other_pos_in_denovo_passes = before.is_some_and(|before| {
        current.pos_metrics.extended.denovo_adj == DenovoAdjecent::ThisIsTheMatchingG
            && before
                .alts
                .iter()
                .find(|alt| alt.metrics.denovo == FormsDenovo::ThisBecomesC)
                .is_some_and(|denovo| denovo.filters.pass(ml_threshold))
    }) || after.is_some_and(|after| {
        current.pos_metrics.extended.denovo_adj == DenovoAdjecent::ThisIsTheMatchingC
            && after
                .alts
                .iter()
                .find(|alt| alt.metrics.denovo == FormsDenovo::ThisBecomesG)
                .is_some_and(|denovo| denovo.filters.pass(ml_threshold))
    });
    Ok(())
}

/// A called `C` alt before a called `G` alt forms a CpG of two variants (`TA>CG`), which the
/// reference context behind [`FormsDenovo`] cannot see: mark both alts as forming it, so that
/// genotyping and methylation weigh their reads as at any de-novo CpG.
///
/// Runs after the alt calls, as it needs the neighbour's, and recomputes each alt from its
/// reference context first, so a re-run after a neighbour's call changed can clear it.
// shortcut: assumes the two alts sit on one chromosome; two hets in trans form no CpG, which
// only phasing could tell.
#[instrument(level = "debug", skip_all)]
pub fn mark_variant_pair_cpgs(
    before: Option<&PileupMetrics>,
    current: &mut PileupMetrics,
    after: Option<&PileupMetrics>,
) -> Result<()> {
    let called = |column: Option<&PileupMetrics>, base| {
        column.is_some_and(|c| {
            c.alts.iter().any(|a| a.base == base && a.call == AltCall::RealVariant)
        })
    };
    let (c_before, g_after) = (called(before, C), called(after, G));
    let reference = current.ref_base();
    for alt in &mut current.alts {
        let from_reference = alt_forms_denovo(alt.base, reference, &current.context);
        let called = alt.call == AltCall::RealVariant;
        alt.metrics.denovo = match alt.base {
            C if called && g_after => FormsDenovo::ThisBecomesC,
            G if called && c_before => FormsDenovo::ThisBecomesG,
            _ => from_reference,
        };
    }
    Ok(())
}
