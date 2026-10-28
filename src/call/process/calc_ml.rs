use crate::metrics::{
    AltVerdict, MetricsForAlt, PileupMetrics, PositionMetrics,
    ml::{
        features,
        types::{ByModel, GpuRastairModel, MachineLearning, MlModel, PlattScaling},
    },
};
use crate::utils::{adjacent, logging::ThisIsABug};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr as _, ensure};
use ndarray::{Array2, s};
use seqair_types::{Base, Probability};
use tracing::debug;

/// Rows a single dispatch can carry, and so the size of the per-model GPU
/// buffers: `max_samples × n_trees × 4 B` of intermediates, 26 MB per model at
/// this value against the 800-tree bundled forests.
///
/// A larger round is not refused, it is split — [`submit_and_collect`] chunks
/// by this — so the only cost of a low value is more rounds, and the only cost
/// of a high one is a reservation that on a discrete GPU is real VRAM. A 100 kb
/// segment produces ~4,900 rows across all models, so a region is never
/// split, and there is room left should the inference queue ever back up
/// enough for regions to coalesce (measured so far: it does not).
pub const GPU_BATCH_BUFFER_SIZE: usize = 8_192;

/// Why an alt gets no row, and so no score.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Unscored {
    /// Too unlikely to be worth the forests: at an `N`, an `N` itself, or at
    /// a column of depth 1 or mean MAPQ 5 or less. `call` marks it
    /// [`crate::vcf::RastairFilter::PreMl`] ([`worth_scoring`]).
    PreMl,
    /// Its features failed to compute or are not all finite: a forest
    /// splits on thresholds, so an infinity would land in a leaf no training
    /// row reached.
    NoFeatures(MlModel),
}

/// The row a forest scores `alt` with: the one gate `call` and `ml train`
/// share, so the forests are fitted on exactly the rows they score.
///
/// `before` and `after` count only when they are the genome's neighbours of
/// `alt`'s column; a filtered region's neighbour in memory may not be.
pub fn ml_row(
    alt: &MetricsForAlt,
    before: Option<&PileupMetrics>,
    after: Option<&PileupMetrics>,
) -> Result<(MlModel, Array2<f32>), Unscored> {
    let column = alt.metrics;
    if !worth_scoring(&column.pos_metrics, column.reference_base, alt.alt.base) {
        return Err(Unscored::PreMl);
    }
    let before = before.filter(|before| adjacent(before, column));
    let after = after.filter(|after| adjacent(column, after));
    let (model, row) = features::features(alt, before, after);
    usable(model, row)
}

/// Whether the forests score `alt` over `reference` at a column with `pos`
/// at all ([`Unscored::PreMl`]).
pub fn worth_scoring(pos: &PositionMetrics, reference: Base, alt: Base) -> bool {
    let likely = pos.depth > 1 && *pos.mapq > 5.;
    likely && reference != Base::Unknown && alt != Base::Unknown
}

fn usable(model: MlModel, row: Result<Array2<f32>>) -> Result<(MlModel, Array2<f32>), Unscored> {
    match row {
        Ok(row) if row.iter().all(|value| value.is_finite()) => Ok((model, row)),
        Ok(_) => Err(Unscored::NoFeatures(model)),
        Err(error) => {
            debug!(%error, model = model.name(), "Failed to calculate features for ML prediction");
            Err(Unscored::NoFeatures(model))
        }
    }
}

/// Build one feature row per ML candidate `only` takes.
///
/// Every alt `only` takes loses any earlier score here, so a rescoring (the
/// SNV cascade) cannot leave an alt it fails to score with the earlier
/// stage's score. Alts [`ml_row`] gives no row stay unscored, which
/// `set_alt_calls` rejects as `low_ml_score`, or as `pre_ml` when they were
/// not [`worth_scoring`].
pub fn extract_ml_rows(
    pileups: &mut [PileupMetrics],
    ml: &MachineLearning,
    only: impl Fn(&PileupMetrics, Base) -> bool,
) -> Result<MlBatch> {
    let model = &ml.model;
    let mut batch = MlBatch::from_fn(|m| ModelBatch::new(features::feature_count(m)));
    let mut taken: Vec<(usize, Base)> = Vec::new();

    for i in 0..pileups.len() {
        let before = i.checked_sub(1).and_then(|before| pileups.get(before));
        let after = pileups.get(i + 1);
        let Some(current) = pileups.get(i) else { continue };

        for alt_base in current.alts.iter().map(|alt| alt.base) {
            if !only(current, alt_base) {
                continue;
            }
            taken.push((i, alt_base));
            let Some(alt) = current.alt_metrics(alt_base) else { continue };
            match ml_row(&alt, before, after) {
                Ok((which, row)) => batch[which]
                    .push(Pending { pileup_idx: i, alt_base, platt: model.platt(which) }, &row)?,
                Err(Unscored::PreMl | Unscored::NoFeatures(_)) => {}
            }
        }
    }

    for (i, alt_base) in taken {
        if let Some(alt) = pileups.get_mut(i).and_then(|p| p.find_alt_mut(alt_base)) {
            alt.verdict = AltVerdict::Unscored;
        }
    }

    Ok(batch)
}

/// Score `pileups` on the calling thread with the flat forests: the CPU twin of
/// the inference thread's scoring.
///
/// Every candidate in the region goes into one `predict` call per model. That
/// batch is what lets biosphere walk 16 samples per tree in lockstep, which is
/// ~5x faster than scoring candidates one at a time. Only the alts `only`
/// takes are scored.
pub fn score_on_cpu(
    pileups: &mut [PileupMetrics],
    ml: &MachineLearning,
    only: impl Fn(&PileupMetrics, Base) -> bool,
) -> Result<()> {
    let model = &ml.model;
    if pileups.is_empty() {
        return Ok(());
    }

    let (rows, targets) = extract_ml_rows(pileups, ml, only)?.split()?;
    let scores = MlScores::from_fn(|m| {
        let rows = rows[m].view();
        if rows.nrows() == 0 {
            return Vec::new();
        }
        model.forest(m).predict(&rows).to_vec()
    });
    apply_ml_scores(pileups, &targets, &scores)
}

/// Score a batch on the GPU: one dispatch per model per [`GPU_BATCH_BUFFER_SIZE`] rows.
///
/// Every model is submitted before any is collected, so their GPU work can
/// overlap.
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
/// A model that returned a different number of scores than it was sent rows
/// is an error, checked before anything is written.
pub fn apply_ml_scores(
    pileups: &mut [PileupMetrics],
    targets: &MlTargets,
    scores: &MlScores,
) -> Result<()> {
    for model in MlModel::ALL {
        ensure!(
            scores[model].len() == targets[model].len(),
            "The {model:?} forest returned {} scores for {} rows",
            scores[model].len(),
            targets[model].len()
        );
    }
    for model in MlModel::ALL {
        for (p, &raw) in targets[model].iter().zip(scores[model].iter()) {
            let Some(pileup) = pileups.get_mut(p.pileup_idx) else { continue };
            let calibrated: Probability = p.platt.calibrate_score(raw);

            if let Some(alt) = pileup.find_alt_mut(p.alt_base) {
                alt.verdict = AltVerdict::Scored(calibrated);
            }
        }
    }
    Ok(())
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
/// has ~5,000 rows across all models, and the bound was two orders of
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
    alt_base: Base,
    platt: PlattScaling,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::{RecordFilters, ml::MachineLearningParams, tests::utils::test_call};
    use crate::pileups;
    use seqair_types::Base::{G, T};

    /// One rule for `call` and `ml train`: a row with an infinity is as
    /// unscored as one with a `NaN` or none at all.
    #[test]
    fn only_finite_rows_are_scored() {
        let row = |values: [f32; 2]| Ok(Array2::from_shape_vec((1, 2), values.to_vec())?);
        assert!(usable(MlModel::Cpg, row([1.0, 2.0])).is_ok());
        for unusable in [
            row([1.0, f32::INFINITY]),
            row([f32::NEG_INFINITY, 2.0]),
            row([f32::NAN, 2.0]),
            Err(color_eyre::eyre::eyre!("no neighbour")),
        ] {
            assert_eq!(usable(MlModel::Cpg, unusable), Err(Unscored::NoFeatures(MlModel::Cpg)));
        }
    }

    #[test]
    fn only_the_alts_the_predicate_takes_are_scored() -> Result<()> {
        let (segment, pileups) = pileups!(
            [ A ] Ref,
            [ T ] OT, [ T ] OB, [ T ] OT,
            [ G ] OT, [ G ] OB, [ G ] OT,
            [ A ] OT, [ A ] OB,
        );
        let mut columns = test_call(segment, pileups, RecordFilters::all())?;
        let column = columns.first_mut().wrap_err("one column")?;
        for base in [T, G] {
            column.find_alt_mut(base).wrap_err("both alts")?.verdict =
                crate::metrics::AltVerdict::Unscored;
        }
        let ml = MachineLearningParams::default().init(1)?.forests;
        score_on_cpu(&mut columns, &ml, |_, alt| alt == T)?;
        let score = |base| columns.first()?.find_alt(base)?.verdict.ml();
        assert!(score(T).is_some(), "the alt the predicate takes is scored");
        assert!(score(G).is_none(), "the alt it refuses keeps its score");
        Ok(())
    }

    #[test]
    fn a_rescored_alt_loses_its_earlier_score_and_short_scores_are_an_error() -> Result<()> {
        let (segment, pileups) = pileups!(
            [ A ] Ref,
            [ T ] OT, [ T ] OB, [ T ] OT,
            [ G ] OT, [ G ] OB, [ G ] OT,
            [ A ] OT, [ A ] OB,
        );
        let mut columns = test_call(segment, pileups, RecordFilters::all())?;
        let earlier = Probability::new(1.0)?;
        let column = columns.first_mut().wrap_err("one column")?;
        for base in [T, G] {
            column.find_alt_mut(base).wrap_err("both alts")?.verdict =
                crate::metrics::AltVerdict::Scored(earlier);
        }
        let ml = MachineLearningParams::default().init(1)?.forests;
        let (rows, targets) = extract_ml_rows(&mut columns, &ml, |_, alt| alt == T)?.split()?;
        let score = |columns: &[PileupMetrics], base| columns.first()?.find_alt(base)?.verdict.ml();
        assert!(score(&columns, T).is_none(), "the alt being rescored is unscored until scored");
        assert_eq!(score(&columns, G), Some(earlier), "the alt not taken keeps its score");

        let sent = MlModel::ALL.into_iter().map(|m| rows[m].nrows()).sum::<usize>();
        assert_eq!(sent, 1, "one row, for T");
        let short = MlScores::from_fn(|_| Vec::new());
        assert!(apply_ml_scores(&mut columns, &targets, &short).is_err());
        assert!(score(&columns, T).is_none(), "nothing is written from a short reply");
        Ok(())
    }
}
