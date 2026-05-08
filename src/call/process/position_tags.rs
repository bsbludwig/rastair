use crate::metrics::{AltCall, PileupMetrics, RecordTags};
use tracing::instrument;

/// Add tags to a position for future filtering
#[instrument(level = "debug", skip_all)]
pub fn add_position_tags(current: &mut PileupMetrics) {
    let tags = RecordTags {
        set: true,
        covered: current.pos_metrics.depth > 0,
        cpg: *current.pos_metrics.cpg,
        denovo_cpg: current
            .alts
            .iter()
            .any(|a| a.call == AltCall::RealVariant && *a.metrics.denovo),
        denovo_cpg_partner: *current.pos_metrics.denovo_adj
            && current.pos_filters.other_pos_in_denovo_passes,
        variant: current.alts.iter().any(|a| a.call == AltCall::RealVariant && !*a.metrics.denovo)
            || !current.indel_calls.is_empty(),
    };

    current.tags = tags;
}
