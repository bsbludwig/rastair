use crate::{
    call::pileup::indels::IndelAllele,
    metrics::{
        MetricsForAlt, MetricsForIndel, PileupMetrics,
        ml::types::{GpuRastairModel, MachineLearning, MlModel, PlattScaling},
    },
    utils::logging::ThisIsABug,
    vcf::{low_ml_score, pre_ml},
};
use color_eyre::eyre::{ContextCompat as _, Result};
use ndarray::{Array2, s};
use seqair_types::{Base, Probability};
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
    score_indels: bool,
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

    // The hard-filter pathway must not carry an ML score: it did not consult the
    // model, and a score would make `low_ml_score` reachable for calls the model
    // never judged.
    let mut indel_scores: Vec<(usize, Probability)> = Vec::new();
    for (i, call) in current.indel_calls.iter().enumerate().filter(|_| score_indels) {
        let m = MetricsForIndel { metrics: current, indel: call };
        if let Some(pred) = ml.predict_indels(&m) {
            indel_scores.push((i, pred.prediction));
        }
    }
    for (i, score) in indel_scores {
        if let Some(call) = current.indel_calls.get_mut(i) {
            call.ml = Some(score);
        }
    }

    Ok(())
}

/// Sequential ML prediction over a Vec of pileups, equivalent to streaming
/// `map_surrounding(add_ml_metrics)`. Used as a CPU fallback when GPU batch
/// prediction is unavailable.
///
/// FIXME: Does this do the same really in regard to matching by position?
pub fn add_ml_metrics_vec(
    pileups: &mut [PileupMetrics],
    ml: &MachineLearning,
    score_indels: bool,
) -> Result<()> {
    for i in 0..pileups.len() {
        let (left, rest) = pileups.split_at_mut(i);
        let (current, right) =
            rest.split_first_mut().wrap_err("Failed to split pileups").this_is_a_bug()?;
        let before = left.last().map(|p| p as &_);
        let after = right.first().map(|p| p as &_);
        add_ml_metrics(before, current, after, ml, score_indels)?;
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
    score_indels: bool,
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

    // Each alt belongs to exactly one model; positions * 4 is the per-model upper bound.
    let max_alts = positions * 4;
    let mut pending = PendingGroups {
        cpg: Vec::with_capacity(positions),
        cpg_features: Array2::zeros((max_alts, feature_num.cpg)),
        cpg_count: 0,
        denovo: Vec::with_capacity(positions),
        denovo_features: Array2::zeros((max_alts, feature_num.denovo_cpg)),
        denovo_count: 0,
        others: Vec::with_capacity(positions),
        others_features: Array2::zeros((max_alts, feature_num.others)),
        others_count: 0,
        insertion: Vec::with_capacity(positions),
        insertion_features: Array2::zeros((max_alts, feature_num.insertion)),
        insertion_count: 0,
        deletion: Vec::with_capacity(positions),
        deletion_features: Array2::zeros((max_alts, feature_num.deletion)),
        deletion_count: 0,
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

                let f = match features_result {
                    Err(error) => {
                        debug!(%error, "Failed to calculate features for ML prediction");
                        continue;
                    }
                    Ok(f) if f.is_any_nan() => continue,
                    Ok(f) => f,
                };

                let pending_item = Pending::snp(i, alt_base, platt);

                match model_type {
                    MlModel::Cpg => {
                        pending.cpg_features.row_mut(pending.cpg_count).assign(&f.row(0));
                        pending.cpg.push(pending_item);
                        pending.cpg_count += 1;
                    }
                    MlModel::DenovoCpg => {
                        pending.denovo_features.row_mut(pending.denovo_count).assign(&f.row(0));
                        pending.denovo.push(pending_item);
                        pending.denovo_count += 1;
                    }
                    MlModel::Others => {
                        pending.others_features.row_mut(pending.others_count).assign(&f.row(0));
                        pending.others.push(pending_item);
                        pending.others_count += 1;
                    }
                    MlModel::Insertion | MlModel::Deletion => {
                        unreachable!("indels are handled in a separate loop below")
                    }
                }
            }

            for (indel_idx, call) in
                pileups_ref[i].indel_calls.iter().enumerate().filter(|_| score_indels)
            {
                let indel_m = MetricsForIndel { metrics: &pileups_ref[i], indel: call };

                let (model_type, platt, features_result) = match &call.allele {
                    IndelAllele::Insertion(_) => (
                        MlModel::Insertion,
                        rastair_model.insertion_platt,
                        calc.calculate_insertion(&indel_m),
                    ),
                    IndelAllele::Deletion(_) => (
                        MlModel::Deletion,
                        rastair_model.deletion_platt,
                        calc.calculate_deletion(&indel_m),
                    ),
                };

                let f = match features_result {
                    Err(error) => {
                        debug!(%error, "Failed to calculate indel features for ML prediction");
                        continue;
                    }
                    Ok(f) if f.is_any_nan() => continue,
                    Ok(f) => f,
                };

                let pending_item = Pending::indel(i, indel_idx, platt);

                match model_type {
                    MlModel::Insertion => {
                        pending
                            .insertion_features
                            .row_mut(pending.insertion_count)
                            .zip_mut_with(&f.row(0), |d, &s| *d = s);
                        pending.insertion.push(pending_item);
                        pending.insertion_count += 1;
                    }
                    MlModel::Deletion => {
                        pending
                            .deletion_features
                            .row_mut(pending.deletion_count)
                            .zip_mut_with(&f.row(0), |d, &s| *d = s);
                        pending.deletion.push(pending_item);
                        pending.deletion_count += 1;
                    }
                    _ => unreachable!("indel alleles always map to Insertion or Deletion"),
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

    // Phase 2b: predict in sub-batches that fit the GPU buffer.
    // Within each round, submit all models before collecting so GPU work overlaps.
    let cpg_features = pending.cpg_features.slice(s![..pending.cpg_count, ..]);
    let denovo_features = pending.denovo_features.slice(s![..pending.denovo_count, ..]);
    let others_features = pending.others_features.slice(s![..pending.others_count, ..]);
    let insertion_features = pending.insertion_features.slice(s![..pending.insertion_count, ..]);
    let deletion_features = pending.deletion_features.slice(s![..pending.deletion_count, ..]);

    let max_count = pending
        .cpg_count
        .max(pending.denovo_count)
        .max(pending.others_count)
        .max(pending.insertion_count)
        .max(pending.deletion_count);
    let mut cpg_preds = Vec::with_capacity(pending.cpg_count);
    let mut denovo_preds = Vec::with_capacity(pending.denovo_count);
    let mut others_preds = Vec::with_capacity(pending.others_count);
    let mut insertion_preds = Vec::with_capacity(pending.insertion_count);
    let mut deletion_preds = Vec::with_capacity(pending.deletion_count);

    for start in (0..max_count).step_by(GPU_BATCH_BUFFER_SIZE) {
        // Submit all models for this sub-batch round concurrently.
        let h_cpg = {
            if start < pending.cpg_count {
                let end = (start + GPU_BATCH_BUFFER_SIZE).min(pending.cpg_count);
                gpu.cpg.predict_submit(&cpg_features.slice(s![start..end, ..]))?
            } else {
                None
            }
        };
        let h_denovo = {
            if start < pending.denovo_count {
                let end = (start + GPU_BATCH_BUFFER_SIZE).min(pending.denovo_count);
                gpu.denovo.predict_submit(&denovo_features.slice(s![start..end, ..]))?
            } else {
                None
            }
        };
        let h_others = {
            if start < pending.others_count {
                let end = (start + GPU_BATCH_BUFFER_SIZE).min(pending.others_count);
                gpu.others.predict_submit(&others_features.slice(s![start..end, ..]))?
            } else {
                None
            }
        };
        let h_insertion = {
            if start < pending.insertion_count {
                let end = (start + GPU_BATCH_BUFFER_SIZE).min(pending.insertion_count);
                gpu.insertion.predict_submit(&insertion_features.slice(s![start..end, ..]))?
            } else {
                None
            }
        };
        let h_deletion = {
            if start < pending.deletion_count {
                let end = (start + GPU_BATCH_BUFFER_SIZE).min(pending.deletion_count);
                gpu.deletion.predict_submit(&deletion_features.slice(s![start..end, ..]))?
            } else {
                None
            }
        };

        // Collect all before next round.
        if let Some(h) = h_cpg {
            cpg_preds.extend(h.collect()?.iter().copied());
        }
        if let Some(h) = h_denovo {
            denovo_preds.extend(h.collect()?.iter().copied());
        }
        if let Some(h) = h_others {
            others_preds.extend(h.collect()?.iter().copied());
        }
        if let Some(h) = h_insertion {
            insertion_preds.extend(h.collect()?.iter().copied());
        }
        if let Some(h) = h_deletion {
            deletion_preds.extend(h.collect()?.iter().copied());
        }
    }

    // Phase 2c: write calibrated predictions back.
    for (items, preds) in [
        (&pending.cpg, &cpg_preds),
        (&pending.denovo, &denovo_preds),
        (&pending.others, &others_preds),
        (&pending.insertion, &insertion_preds),
        (&pending.deletion, &deletion_preds),
    ] {
        for (p, &raw_pred) in items.iter().zip(preds.iter()) {
            let calibrated: Probability = p.platt.calibrate_score(f64::from(raw_pred));
            let threshold = ml.threshold;
            if let Some(indel_idx) = p.indel_idx {
                if let Some(call) = pileups[p.pileup_idx].indel_calls.get_mut(indel_idx) {
                    call.ml = Some(calibrated);
                }
            } else if let Some(filters) = pileups[p.pileup_idx].alt_filters_mut(p.alt_base) {
                filters.ml.replace(calibrated);
                filters.filters.add(low_ml_score, move || calibrated < threshold);
            }
        }
    }

    Ok(())
}

struct PendingGroups {
    cpg: Vec<Pending>,
    cpg_features: Array2<f32>,
    cpg_count: usize,
    denovo: Vec<Pending>,
    denovo_features: Array2<f32>,
    denovo_count: usize,
    others: Vec<Pending>,
    others_features: Array2<f32>,
    others_count: usize,
    insertion: Vec<Pending>,
    insertion_features: Array2<f32>,
    insertion_count: usize,
    deletion: Vec<Pending>,
    deletion_features: Array2<f32>,
    deletion_count: usize,
}

/// Each pending prediction records where to write the result back.
struct Pending {
    pileup_idx: usize,
    /// For SNPs: the alt base. For indels: unused.
    alt_base: Base,
    /// For indels: index into `indel_calls`. For SNPs: `None`.
    indel_idx: Option<usize>,
    platt: PlattScaling,
}

impl Pending {
    fn snp(pileup_idx: usize, alt_base: Base, platt: PlattScaling) -> Self {
        Self { pileup_idx, alt_base, indel_idx: None, platt }
    }

    fn indel(pileup_idx: usize, indel_idx: usize, platt: PlattScaling) -> Self {
        Self { pileup_idx, alt_base: Base::Unknown, indel_idx: Some(indel_idx), platt }
    }
}
