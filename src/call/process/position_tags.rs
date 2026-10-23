use crate::metrics::{AltCall, PileupMetrics, RecordTags};
use tracing::instrument;

/// Add tags to a position for future filtering
#[instrument(level = "debug", skip_all)]
pub fn add_position_tags(current: &mut PileupMetrics) {
    let tags = RecordTags {
        set: true,
        covered: current.pos_metrics.depth > 0,
        denovo_cpg: current
            .alts
            .iter()
            .any(|a| a.call == AltCall::RealVariant && *a.metrics.denovo),
        denovo_cpg_partner: *current.pos_metrics.denovo_adj
            && current.pos_filters.other_pos_in_denovo_passes,
        variant: current.alts.iter().any(|a| a.call == AltCall::RealVariant && !*a.metrics.denovo)
            || current
                .indel_data
                .as_ref()
                .and_then(|d| d.scored.as_ref())
                .is_some_and(|s| s.records.calls().iter().any(|call| call.is_called())),
    };

    current.tags = tags;
}
