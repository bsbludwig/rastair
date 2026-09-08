use crate::utils::logging::ThisIsABug;
use crate::{
    call::pileup::indels::IndelAllele,
    metrics::{
        MetricsForAlt, MetricsForIndel, PileupMetrics,
        ml::types::{
            ByModel, GpuRastairModel, MachineLearning, MlModel, PlattScaling, RastairFlatModel,
        },
    },
    vcf::{low_ml_score, pre_ml},
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr as _, ensure};
use ndarray::{Array2, s};
use seqair_types::{Base, Probability};
use tracing::debug;
#[cfg(test)]
use tracing::instrument;

/// Rows a single dispatch can carry, and so the size of the per-model GPU
/// buffers: `max_samples × n_trees × 4 B` of intermediates, 26 MB per model at
/// this value against the 800-tree bundled forests.
///
/// A larger round is not refused, it is split — [`submit_and_collect`] chunks
/// by this — so the only cost of a low value is more rounds, and the only cost
/// of a high one is a reservation that on a discrete GPU is real VRAM. A 100 kb
/// segment produces ~4,900 rows across all five models, so a region is never
/// split, and there is room left should the inference queue ever back up
/// enough for regions to coalesce (measured so far: it does not).
pub const GPU_BATCH_BUFFER_SIZE: usize = 8_192;

/// Filter out very unlikely alts before running slow ML
fn pre_ml_filter(c: &MetricsForAlt) -> bool {
    c.metrics.pos_metrics.depth > 1 && *c.metrics.pos_metrics.mapq > 5.
}

/// One candidate at a time, the way scoring worked before regions were batched.
/// Kept as the reference the batched paths are tested against.
#[cfg(test)]
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

/// Build one feature row per ML candidate, and tag the alts not worth scoring.
///
/// Reads `pileups[i - 1]` and `pileups[i + 1]`, so this has to run where the
/// region's `Vec` lives. Candidates whose features fail to compute or come out
/// `NaN` are dropped here and simply never receive a score.
pub fn extract_ml_rows(
    pileups: &mut [PileupMetrics],
    ml: &MachineLearning,
    model: &RastairFlatModel,
    score_indels: bool,
) -> Result<MlBatch> {
    let calc = &ml.feature_calculator;
    let feature_num = calc.feature_num();

    let mut batch = MlBatch::from_fn(|m| ModelBatch::new(feature_num.get(m)));
    let mut pre_ml_rejected: Vec<(usize, Base)> = Vec::new();

    for i in 0..pileups.len() {
        let before = if i > 0 { pileups.get(i - 1) } else { None };
        let after = pileups.get(i + 1);
        let Some(current) = pileups.get(i) else { continue };

        for alt_base in current.alts() {
            let Some(alt) = current.alt_metrics(alt_base) else { continue };

            if !pre_ml_filter(&alt) {
                pre_ml_rejected.push((i, alt_base));
                continue;
            }

            let (which, features) = if alt.is_evidence_for_methylation() {
                (MlModel::Cpg, calc.calculate_cpg(&alt, before, after))
            } else if *alt.alt.denovo {
                (MlModel::DenovoCpg, calc.calculate_denovo_cpg(&alt, before, after))
            } else {
                (MlModel::Others, calc.calculate_others(&alt, before, after))
            };

            let Some(f) = usable_features(features, "alt") else { continue };
            batch[which].push(Pending::snp(i, alt_base, model.platt(which)), &f)?;
        }

        for (indel_idx, call) in current.indel_calls.iter().enumerate().filter(|_| score_indels) {
            let indel = MetricsForIndel { metrics: current, indel: call };

            let (which, features) = match &call.allele {
                IndelAllele::Insertion(_) => (MlModel::Insertion, calc.calculate_insertion(&indel)),
                IndelAllele::Deletion(_) => (MlModel::Deletion, calc.calculate_deletion(&indel)),
            };

            let Some(f) = usable_features(features, "indel") else { continue };
            batch[which].push(Pending::indel(i, indel_idx, model.platt(which)), &f)?;
        }
    }

    for (i, alt_base) in pre_ml_rejected {
        if let Some(filters) = pileups.get_mut(i).and_then(|p| p.alt_filters_mut(alt_base)) {
            filters.filters.add(pre_ml, || true);
        }
    }

    Ok(batch)
}

/// The one feature row a calculator produced, unless it failed or is `NaN`.
///
/// A `NaN` anywhere in the row would propagate through the forest into the
/// score, so such a candidate is left unscored rather than scored wrongly.
fn usable_features(features: Result<Array2<f32>>, kind: &'static str) -> Option<Array2<f32>> {
    match features {
        Err(error) => {
            debug!(%error, kind, "Failed to calculate features for ML prediction");
            None
        }
        Ok(f) if f.is_any_nan() => None,
        Ok(f) => Some(f),
    }
}

/// Score `pileups` on the calling thread with the flat forests: the CPU twin of
/// the inference thread's scoring.
///
/// Every candidate in the region goes into one `predict` call per model. That
/// batch is what lets biosphere walk 16 samples per tree in lockstep, which is
/// ~5x faster than scoring candidates one at a time as [`add_ml_metrics`]
/// does.
pub fn score_on_cpu(
    pileups: &mut [PileupMetrics],
    ml: &MachineLearning,
    score_indels: bool,
) -> Result<()> {
    let Some(model) = ml.model.as_ref() else {
        return Ok(());
    };
    if pileups.is_empty() {
        return Ok(());
    }

    let (rows, targets) = extract_ml_rows(pileups, ml, model, score_indels)?.split()?;
    let scores = MlScores::from_fn(|m| {
        let rows = rows[m].view();
        if rows.nrows() == 0 {
            return Vec::new();
        }
        model.forest(m).predict(&rows).to_vec()
    });
    apply_ml_scores(pileups, &targets, &scores, ml.threshold);
    Ok(())
}

/// Score a batch on the GPU: one dispatch per model per [`GPU_BATCH_BUFFER_SIZE`] rows.
///
/// Every model is submitted before any is collected. The five forests sit on
/// five separate devices, so this is what lets their GPU work overlap.
pub fn submit_and_collect(rows: &MlRows, gpu: &GpuRastairModel) -> Result<MlScores> {
    let mut scores = MlScores::from_fn(|m| Vec::with_capacity(rows[m].nrows()));
    let longest = MlModel::ALL.into_iter().map(|m| rows[m].nrows()).max().unwrap_or(0);

    for start in (0..longest).step_by(GPU_BATCH_BUFFER_SIZE) {
        let mut handles = ByModel::from_fn(|_| None);

        for model in MlModel::ALL {
            let rows = rows[model].view();
            if start >= rows.nrows() {
                continue;
            }
            let end = (start + GPU_BATCH_BUFFER_SIZE).min(rows.nrows());
            handles[model] = gpu.forest(model).predict_submit(&rows.slice(s![start..end, ..]))?;
        }

        for model in MlModel::ALL {
            if let Some(handle) = handles[model].take() {
                scores[model].extend(handle.collect()?.iter().map(|&score| f64::from(score)));
            }
        }
    }

    Ok(scores)
}

/// Platt-calibrate the raw scores and write them back into the region.
///
/// A candidate whose score is missing — the batch was truncated, or the model
/// returned fewer rows than were submitted — is skipped, not defaulted: `zip`
/// stops at the shorter side.
pub fn apply_ml_scores(
    pileups: &mut [PileupMetrics],
    targets: &MlTargets,
    scores: &MlScores,
    threshold: Probability,
) {
    for model in MlModel::ALL {
        for (p, &raw) in targets[model].iter().zip(scores[model].iter()) {
            let Some(pileup) = pileups.get_mut(p.pileup_idx) else { continue };
            let calibrated: Probability = p.platt.calibrate_score(raw);

            if let Some(indel_idx) = p.indel_idx {
                if let Some(call) = pileup.indel_calls.get_mut(indel_idx) {
                    call.ml = Some(calibrated);
                }
            } else if let Some(filters) = pileup.alt_filters_mut(p.alt_base) {
                filters.ml.replace(calibrated);
                filters.filters.add(low_ml_score, move || calibrated < threshold);
            }
        }
    }
}

/// One region's ML work: a feature row per candidate, grouped by the model that
/// scores it.
///
/// Holds no borrow of the region, so it can be handed to another thread.
pub type MlBatch = ByModel<ModelBatch>;

/// Raw, uncalibrated forest scores for an [`MlBatch`], in the same row order.
///
/// `f64` because that is what [`FlatForest::predict`] returns and what
/// [`PlattScaling::calibrate_score`] takes; the CPU path then calibrates the
/// very same value the one-candidate-at-a-time path did, so batching does not
/// change a single output. The GPU's `f32` scores widen losslessly.
///
/// [`FlatForest::predict`]: biosphere::FlatForest::predict
pub type MlScores = ByModel<Vec<f64>>;

/// The half of an [`MlBatch`] the GPU needs: a filled feature matrix per model.
pub type MlRows = ByModel<Array2<f32>>;

/// The half of an [`MlBatch`] the GPU does not need: where each row's score
/// goes. Row `i` of `MlRows[m]` scores `MlTargets[m][i]`.
pub type MlTargets = ByModel<Vec<Pending>>;

impl MlBatch {
    /// Split off the rows, so they can be scored somewhere the region is not.
    ///
    /// The row-to-target alignment [`ModelBatch`] maintains is what survives
    /// the split, so this is the only place allowed to take the two apart.
    pub fn split(self) -> Result<(MlRows, MlTargets)> {
        let mut rows = MlRows::from_fn(|_| Array2::zeros((0, 0)));
        let mut targets = MlTargets::from_fn(|_| Vec::new());
        for (model, batch) in self {
            (rows[model], targets[model]) = batch.split()?;
        }
        Ok((rows, targets))
    }
}

/// Feature rows for one model, and where each row's score has to be written back.
///
/// Rows are appended to a flat `Vec` that grows with the candidates actually
/// found, rather than reserved at `positions × 4` up front: a 100 kb region
/// has ~5,000 rows across all five models, and the bound was two orders of
/// magnitude above that. `pending.len() × n_features` is always the length of
/// `features`, which [`Self::push`] enforces, so the two cannot drift apart.
pub struct ModelBatch {
    pending: Vec<Pending>,
    features: Vec<f32>,
    n_features: usize,
}

impl ModelBatch {
    fn new(n_features: usize) -> Self {
        Self { pending: Vec::new(), features: Vec::new(), n_features }
    }

    fn push(&mut self, item: Pending, features: &Array2<f32>) -> Result<()> {
        let row = features.as_slice().wrap_err("Feature row is not contiguous").this_is_a_bug()?;
        ensure!(
            row.len() == self.n_features,
            "Feature row has {} values, model expects {}",
            row.len(),
            self.n_features
        );
        self.features.extend_from_slice(row);
        self.pending.push(item);
        Ok(())
    }

    fn split(self) -> Result<(Array2<f32>, Vec<Pending>)> {
        let rows = Array2::from_shape_vec((self.pending.len(), self.n_features), self.features)
            .wrap_err("Feature rows do not match the number of candidates")
            .this_is_a_bug()?;
        Ok((rows, self.pending))
    }
}

pub struct Pending {
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
