use crate::{
    metrics::{
        MetricsForAlt, PileupMetrics,
        ml::types::{GpuRastairModel, MachineLearning, MlModel, PlattScaling},
    },
    utils::logging::ThisIsABug,
    vcf::{low_ml_score, pre_ml},
};
use biosphere::gpu::PredictHandle;
use color_eyre::eyre::{ContextCompat as _, Result};
use rastair_types::{Base, Probability};
use tracing::{debug, instrument};

/// Size buffer for reasonably full chunk (10k positions × 4 alts max) per thread.
pub const GPU_BATCH_BUFFER_SIZE: usize = 40_000;

/// Filter out very unlikely alts before running slow ML
fn pre_ml_filter(c: &MetricsForAlt) -> bool {
    c.metrics.pos_metrics.depth > 1 && *c.metrics.pos_metrics.mapq > 5.
}

#[instrument(level = "debug", skip_all)]
pub fn add_ml_metrics(
    before: Option<&PileupMetrics>,
    current: &mut PileupMetrics,
    after: Option<&PileupMetrics>,
    ml: &MachineLearning,
) -> Result<()> {
    if !ml.enabled() {
        return Ok(());
    }

    'alts: for alt_base in current.alts() {
        let alt =
            current.alt_metrics(alt_base).wrap_err("Failed to get alt metrics").this_is_a_bug()?;

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

    Ok(())
}

/// Sequential ML prediction over a Vec of pileups, equivalent to streaming
/// `map_surrounding(add_ml_metrics)`. Used as a CPU fallback when GPU batch
/// prediction is unavailable.
///
/// FIXME: Does this do the same really in regard to matching by position?
pub fn add_ml_metrics_vec(pileups: &mut [PileupMetrics], ml: &MachineLearning) -> Result<()> {
    for i in 0..pileups.len() {
        let (left, rest) = pileups.split_at_mut(i);
        let (current, right) =
            rest.split_first_mut().wrap_err("Failed to split pileups").this_is_a_bug()?;
        let before = left.last().map(|p| p as &_);
        let after = right.first().map(|p| p as &_);
        add_ml_metrics(before, current, after, ml)?;
    }
    Ok(())
}

/// Batch GPU ML prediction over a full chunk of pileups.
///
/// Groups all passing alts by model type (CpG / de-novo CpG / others), runs
/// a single `GpuForest::predict` call per model, and writes the Platt-calibrated
/// predictions back into the pileup filter state. This is 10k-30× fewer GPU
/// dispatches than the per-alt streaming approach.
pub fn batch_add_ml_metrics(
    pileups: &mut [PileupMetrics],
    ml: &MachineLearning,
    gpu: &GpuRastairModel,
) -> Result<()> {
    let Some(rastair_model) = ml.model.as_ref() else {
        return Ok(());
    };
    if pileups.is_empty() {
        return Ok(());
    }
    let positions = pileups.len();

    let calc = &ml.feature_calculator;
    let feature_num = ml.feature_calculator.feature_num();

    let mut pending = PendingGroups {
        cpg: Vec::with_capacity(positions / 2),
        cpg_features: Vec::with_capacity(positions / 2 * feature_num.cpg),
        denovo: Vec::with_capacity(positions / 2),
        denovo_features: Vec::with_capacity(positions / 2 * feature_num.denovo_cpg),
        others: Vec::with_capacity(positions / 2),
        others_features: Vec::with_capacity(positions / 2 * feature_num.others),
    };
    // (pileup_idx, alt_base) pairs that failed pre_ml_filter
    let mut pre_ml_rejected: Vec<(usize, Base)> = Vec::new();

    // Phase 1: read-only pass — compute features for all alts
    {
        let pileups_ref: &[PileupMetrics] = pileups;
        for i in 0..pileups_ref.len() {
            let before = if i > 0 { Some(&pileups_ref[i - 1]) } else { None };
            let after = pileups_ref.get(i + 1);

            // alts() returns a Vec<Base>; collect into local so we can reborrow
            // pileups_ref[i] for alt_metrics inside the loop.
            let alt_bases = pileups_ref[i].alts();

            for alt_base in alt_bases {
                let Some(alt) = pileups_ref[i].alt_metrics(alt_base) else { continue };

                if !pre_ml_filter(&alt) {
                    pre_ml_rejected.push((i, alt_base));
                    continue;
                }

                let (model_type, platt, features_result) = if alt.is_evidence_for_methylation() {
                    (MlModel::Cpg, rastair_model.cpg_platt, calc.calculate_cpg(&alt, before, after))
                } else if *alt.alt.denovo {
                    (
                        MlModel::DenovoCpg,
                        rastair_model.denovo_platt,
                        calc.calculate_denovo_cpg(&alt, before, after),
                    )
                } else {
                    (
                        MlModel::Others,
                        rastair_model.others_platt,
                        calc.calculate_others(&alt, before, after),
                    )
                };

                let features: Vec<f32> = match features_result {
                    Ok(f) => {
                        if f.is_any_nan() {
                            continue;
                        }
                        f.row(0).iter().map(|&v| v as f32).collect()
                    }
                    Err(error) => {
                        debug!(%error, "Failed to calculate features for ML prediction");
                        continue;
                    }
                };

                let pending_item = Pending { pileup_idx: i, alt_base, platt };

                match model_type {
                    MlModel::Cpg => {
                        pending.cpg.push(pending_item);
                        pending.cpg_features.extend(features);
                    }
                    MlModel::DenovoCpg => {
                        pending.denovo.push(pending_item);
                        pending.denovo_features.extend(features);
                    }
                    MlModel::Others => {
                        pending.others.push(pending_item);
                        pending.others_features.extend(features);
                    }
                }
            }
        }
    } // end read-only borrow of pileups

    // Phase 2a: apply pre_ml filter tags
    for (i, alt_base) in pre_ml_rejected {
        if let Some(filters) = pileups[i].alt_filters_mut(alt_base) {
            filters.filters.add(pre_ml, || true);
        }
    }

    // Phase 2b: submit all three GPU dispatches, then collect — GPU work overlaps.
    let h_cpg = gpu.cpg.predict_submit(&pending.cpg_features, pending.cpg.len());
    let h_denovo = gpu.denovo.predict_submit(&pending.denovo_features, pending.denovo.len());
    let h_others = gpu.others.predict_submit(&pending.others_features, pending.others.len());

    let cpg_preds = collect_handle(h_cpg);
    let denovo_preds = collect_handle(h_denovo);
    let others_preds = collect_handle(h_others);

    for (items, preds) in [
        (&pending.cpg, &cpg_preds),
        (&pending.denovo, &denovo_preds),
        (&pending.others, &others_preds),
    ] {
        for (p, &raw_pred) in items.iter().zip(preds.iter()) {
            let calibrated: Probability = p.platt.calibrate_score(raw_pred as f64);
            let threshold = ml.threshold;
            if let Some(filters) = pileups[p.pileup_idx].alt_filters_mut(p.alt_base) {
                filters.ml.replace(calibrated);
                filters.filters.add(low_ml_score, move || calibrated < threshold);
            }
        }
    }

    Ok(())
}

fn collect_handle(handle: Option<PredictHandle<'_>>) -> Vec<f32> {
    handle.map(|h| h.collect()).unwrap_or_default()
}

struct PendingGroups {
    cpg: Vec<Pending>,
    cpg_features: Vec<f32>,
    denovo: Vec<Pending>,
    denovo_features: Vec<f32>,
    others: Vec<Pending>,
    others_features: Vec<f32>,
}

/// Each pending prediction records where to write the result back.
struct Pending {
    pileup_idx: usize,
    alt_base: Base,
    platt: PlattScaling,
}
