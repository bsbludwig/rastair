use crate::{
    metrics::{
        MetricsForAlt, PileupMetrics,
        ml::types::{GpuRastairModel, MachineLearning, MlModel, PlattScaling},
    },
    utils::logging::ThisIsABug,
    vcf::{low_ml_score, pre_ml},
};
use color_eyre::eyre::{ContextCompat as _, Result};
use rastair_types::{Base, Probability};
use tracing::{debug, instrument};

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
pub fn add_ml_metrics_vec(pileups: &mut [PileupMetrics], ml: &MachineLearning) -> Result<()> {
    for i in 0..pileups.len() {
        let (left, rest) = pileups.split_at_mut(i);
        let (current, right) = rest.split_first_mut().expect("index in bounds");
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
    pileups: &mut Vec<PileupMetrics>,
    ml: &MachineLearning,
    gpu: &GpuRastairModel,
) -> Result<()> {
    let Some(rastair_model) = ml.model.as_ref() else {
        return Ok(());
    };
    if pileups.is_empty() {
        return Ok(());
    }

    let calc = ml.feature_calculator.get_calculator();

    // Each pending prediction records where to write the result back.
    struct Pending {
        pileup_idx: usize,
        alt_base: Base,
        model: MlModel,
        features: Vec<f32>, // pre-converted to f32 for direct GPU upload
        platt: PlattScaling,
    }

    let mut pending: Vec<Pending> = Vec::new();
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

                let features = match features_result {
                    Ok(f) => {
                        if f.is_any_nan() {
                            continue;
                        }
                        f.row(0).iter().map(|&v| v as f32).collect()
                    }
                    Err(_) => continue,
                };

                pending.push(Pending {
                    pileup_idx: i,
                    alt_base,
                    model: model_type,
                    features,
                    platt,
                });
            }
        }
    } // end read-only borrow of pileups

    // Phase 2a: apply pre_ml filter tags
    for (i, alt_base) in pre_ml_rejected {
        if let Some(filters) = pileups[i].alt_filters_mut(alt_base) {
            filters.filters.add(pre_ml, || true);
        }
    }

    // Phase 2b: batch GPU predict per model type and assign calibrated predictions
    for (model_type, gpu_forest) in [
        (MlModel::Cpg, &gpu.cpg),
        (MlModel::DenovoCpg, &gpu.denovo),
        (MlModel::Others, &gpu.others),
    ] {
        let batch: Vec<&Pending> = pending.iter().filter(|p| p.model == model_type).collect();
        if batch.is_empty() {
            continue;
        }

        // Stack all feature rows into a flat f32 slice (row-major)
        let features_flat: Vec<f32> =
            batch.iter().flat_map(|p| p.features.iter().copied()).collect();
        let gpu_preds = gpu_forest.predict(&features_flat, batch.len());

        for (p, raw_pred) in batch.iter().zip(gpu_preds.iter()) {
            let calibrated: Probability = p.platt.calibrate_score(*raw_pred as f64);
            let threshold = ml.threshold;
            if let Some(filters) = pileups[p.pileup_idx].alt_filters_mut(p.alt_base) {
                filters.ml.replace(calibrated);
                filters.filters.add(low_ml_score, move || calibrated < threshold);
            }
        }
    }

    Ok(())
}
