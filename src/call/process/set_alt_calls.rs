use crate::metrics::{AltCall, PileupMetrics};
use color_eyre::eyre::Result;
use seqair_types::{
    Base::{self, *},
    Probability, SmallVec,
};
use tracing::instrument;

/// Call each alt!
#[instrument(level = "debug", skip_all)]
pub fn set_alt_calls(current: &mut PileupMetrics, ml_threshold: Option<Probability>) -> Result<()> {
    let alleles: SmallVec<Base, 4> =
        std::iter::once(current.ref_base()).chain(current.alts.iter().map(|a| a.base)).collect();
    let cpg_context = *current.pos_metrics.cpg || current.forms_denovo();

    for alt in &mut current.alts {
        if alt.filters.pass(ml_threshold) {
            alt.call = AltCall::RealVariant;
        } else if cpg_context
            && let Some(ref_base) = methylation_evidence(alt.base)
            && alleles.contains(&ref_base)
        {
            alt.call = AltCall::MethylationEvidenceOnly { for_base: ref_base }
        } else {
            alt.call = AltCall::ReadError;
        }
        // TODO: Cover case where this is the other position of a passing de-novo CpG
    }

    Ok(())
}

fn methylation_evidence(base: Base) -> Option<Base> {
    match base {
        A => Some(G),
        T => Some(C),
        _ => None,
    }
}
