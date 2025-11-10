use crate::{
    metrics::{MetricsForAlt, PileupMetrics, ml::types::MachineLearning},
    utils::{Surrounding, logging::ThisIsABug, surrounding_pileups},
    vcf::{low_ml_score, pre_ml},
};
use color_eyre::eyre::{ContextCompat as _, Result};
use tracing::{debug, instrument};

/// Filter out very unlikely alts before running slow ML
fn pre_ml_filter(c: &MetricsForAlt) -> bool {
    c.metrics.pos_metrics.depth > 1 && *c.metrics.pos_metrics.mapq > 5.
}

#[instrument(level = "info", skip_all)]
pub fn add_ml_metrics(pileups: &mut [PileupMetrics], ml: &MachineLearning) -> Result<()> {
    if ml.disabled {
        return Ok(());
    }

    let pileups_len = pileups.len();
    for i in 0..pileups_len {
        let Surrounding { before, current, after } = surrounding_pileups(pileups, i);

        'alts: for alt_base in current.alts() {
            let alt = current
                .alt_metrics(alt_base)
                .wrap_err("Failed to get alt metrics")
                .this_is_a_bug()?;

            if !pre_ml_filter(&alt) {
                let filters = current
                    .alt_filters_mut(alt_base)
                    .wrap_err("Failed to get mutable alt metrics")
                    .this_is_a_bug()?;
                filters.filters.add(pre_ml, || true);

                // Skip expensive ML prediction for this low-quality alt
                continue 'alts;
            }

            if let Some(prediction) = ml.predict(&alt, before, after) {
                let filters = current
                    .alt_filters_mut(alt_base)
                    .wrap_err("Failed to get mutable alt metrics")
                    .this_is_a_bug()?;
                filters.ml.replace(prediction.prediction);
                filters.filters.add(low_ml_score, || !prediction.pass());
            } else {
                debug!(
                    pos=%current.pos(),
                    ref_base=%current.ref_base(),
                    alt_base=%alt_base,
                    "No ML prediction made"
                );
            }
        }
    }

    Ok(())
}
