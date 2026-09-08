use crate::{
    call::pileup::indels::IndelAllele,
    metrics::{
        MetricsForAlt, MetricsForIndel, PileupMetrics,
        ml::types::{
            ByModel, GpuRastairModel, MachineLearning, MlModel, PlattScaling, RastairFlatModel,
        },
    },
    utils::logging::ThisIsABug,
    vcf::RastairFilter,
};
use color_eyre::eyre::{ContextCompat as _, Result};
use ndarray::{Array2, ArrayView2, s};
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
            filters.filters.add(RastairFilter::PreMl, || true);

            // Skip expensive ML prediction for this low-quality alt
            continue 'alts;
        }

        if let Some(prediction) = ml.predict(&alt, before, after) {
            let filters = current
                .alt_filters_mut(alt_base)
                .wrap_err("Failed to get mutable alt metrics")
                .this_is_a_bug()?;
            filters.ml.replace(prediction.prediction);
            filters.filters.add(RastairFilter::LowMlScore, || !prediction.pass());
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
    if let Some(ref d) = current.indel_data {
        for (i, call) in d.calls.iter().enumerate().filter(|_| score_indels) {
            let m = MetricsForIndel { metrics: current, indel: call };
            if let Some(pred) = ml.predict_indels(&m) {
                indel_scores.push((i, pred.prediction));
            }
        }
    }
    for (i, score) in indel_scores {
        if let Some(ref mut d) = current.indel_data {
            if let Some(call) = d.calls.get_mut(i) {
                call.ml = Some(score);
            }
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
/// Three stages, deliberately separate: [`extract_ml_rows`] needs the region's
/// pileups and their neighbours, [`submit_and_collect`] needs nothing but `f32`
/// rows, and [`apply_ml_scores`] needs the region again. Only the middle one
/// touches the GPU.
pub fn batch_add_ml_metrics(
    pileups: &mut [PileupMetrics],
    ml: &MachineLearning,
    gpu: &GpuRastairModel,
    score_indels: bool,
) -> Result<()> {
    let Some(model) = ml.model.as_ref() else {
        return Ok(());
    };
    if pileups.is_empty() {
        return Ok(());
    }

    let batch = extract_ml_rows(pileups, ml, model, score_indels);
    let scores = submit_and_collect(&batch, gpu)?;
    apply_ml_scores(pileups, &batch, &scores, ml.threshold);
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
) -> MlBatch {
    let calc = &ml.feature_calculator;
    let feature_num = calc.feature_num();

    // Each candidate belongs to exactly one model, so `positions * 4` bounds
    // every per-model row count at once.
    let max_rows = pileups.len() * 4;
    let mut batch = MlBatch::from_fn(|m| ModelBatch::new(max_rows, feature_num.get(m)));
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
            batch[which].push(Pending::snp(i, alt_base, model.platt(which)), &f);
        }

        let indel_calls = current.indel_data.as_deref().map_or(&[][..], |d| &d.calls);
        for (indel_idx, call) in indel_calls.iter().enumerate().filter(|_| score_indels) {
            let indel = MetricsForIndel { metrics: current, indel: call };

            let (which, features) = match &call.allele {
                IndelAllele::Insertion(_) => (MlModel::Insertion, calc.calculate_insertion(&indel)),
                IndelAllele::Deletion(_) => (MlModel::Deletion, calc.calculate_deletion(&indel)),
            };

            let Some(f) = usable_features(features, "indel") else { continue };
            batch[which].push(Pending::indel(i, indel_idx, model.platt(which)), &f);
        }
    }

    for (i, alt_base) in pre_ml_rejected {
        if let Some(filters) = pileups.get_mut(i).and_then(|p| p.alt_filters_mut(alt_base)) {
            filters.filters.add(RastairFilter::PreMl, || true);
        }
    }

    batch
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

/// Score a batch on the GPU: one dispatch per model per [`GPU_BATCH_BUFFER_SIZE`] rows.
///
/// Every model is submitted before any is collected. The five forests sit on
/// five separate devices, so this is what lets their GPU work overlap.
pub fn submit_and_collect(batch: &MlBatch, gpu: &GpuRastairModel) -> Result<MlScores> {
    let mut scores = MlScores::from_fn(|m| Vec::with_capacity(batch[m].len()));
    let longest = MlModel::ALL.into_iter().map(|m| batch[m].len()).max().unwrap_or(0);

    for start in (0..longest).step_by(GPU_BATCH_BUFFER_SIZE) {
        let mut handles = ByModel::from_fn(|_| None);

        for model in MlModel::ALL {
            let rows = batch[model].rows();
            if start >= rows.nrows() {
                continue;
            }
            let end = (start + GPU_BATCH_BUFFER_SIZE).min(rows.nrows());
            handles[model] = gpu.forest(model).predict_submit(&rows.slice(s![start..end, ..]))?;
        }

        for model in MlModel::ALL {
            if let Some(handle) = handles[model].take() {
                scores[model].extend(handle.collect()?.iter().copied());
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
    batch: &MlBatch,
    scores: &MlScores,
    threshold: Probability,
) {
    for model in MlModel::ALL {
        for (p, &raw) in batch[model].pending.iter().zip(scores[model].iter()) {
            let Some(pileup) = pileups.get_mut(p.pileup_idx) else { continue };
            let calibrated: Probability = p.platt.calibrate_score(f64::from(raw));

            if let Some(indel_idx) = p.indel_idx {
                if let Some(d) = pileup.indel_data.as_mut()
                    && let Some(call) = d.calls.get_mut(indel_idx)
                {
                    call.ml = Some(calibrated);
                }
            } else if let Some(filters) = pileup.alt_filters_mut(p.alt_base) {
                filters.ml.replace(calibrated);
                filters.filters.add(RastairFilter::LowMlScore, move || calibrated < threshold);
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
pub type MlScores = ByModel<Vec<f32>>;

/// Feature rows for one model, and where each row's score has to be written back.
///
/// `features` is allocated to an upper bound; `pending.len()` is how many of its
/// rows are filled, so the two cannot drift apart.
pub struct ModelBatch {
    pending: Vec<Pending>,
    features: Array2<f32>,
}

impl ModelBatch {
    fn new(max_rows: usize, n_features: usize) -> Self {
        Self { pending: Vec::new(), features: Array2::zeros((max_rows, n_features)) }
    }

    fn push(&mut self, item: Pending, features: &Array2<f32>) {
        self.features.row_mut(self.pending.len()).assign(&features.row(0));
        self.pending.push(item);
    }

    fn len(&self) -> usize {
        self.pending.len()
    }

    fn rows(&self) -> ArrayView2<'_, f32> {
        self.features.slice(s![..self.pending.len(), ..])
    }
}

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
