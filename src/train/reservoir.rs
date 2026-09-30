//! Bounded collection of labelled candidates, and the train/holdout draw.
//!
//! A genome offers far more candidates than training uses, so each model keeps
//! a uniform sample per class, sized to what the draw will ask for plus a full
//! holdout. Which candidates survive is decided by a seeded random key, so the
//! sample depends on `--seed` alone and not on how segments were scheduled.

use crate::metrics::ml::types::MlModel;
use color_eyre::eyre::{Context as _, ContextCompat as _, Result, ensure};
use ndarray::{Array1, Array2};
use rand::prelude::*;
use seqair_types::SmolStr;
use std::cmp::Ordering;
use tracing::warn;

/// Most examples held back from training for Platt calibration, and the
/// headroom per class a reservoir keeps beyond the training draw.
const MAX_HOLDOUT: usize = 100_000;

/// Whether the truth set claims a candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Label {
    Positive,
    Negative,
}

impl Label {
    pub(super) const fn of(in_truth: bool) -> Self {
        if in_truth { Self::Positive } else { Self::Negative }
    }

    const fn is_positive(self) -> bool {
        matches!(self, Self::Positive)
    }

    /// The target value biosphere fits on.
    pub(super) const fn weight(self) -> f64 {
        match self {
            Self::Positive => 1.0,
            Self::Negative => 0.0,
        }
    }
}

/// One value per [`Label`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(super) struct ByLabel<T> {
    pub(super) positive: T,
    pub(super) negative: T,
}

impl<T> ByLabel<T> {
    fn get_mut(&mut self, label: Label) -> &mut T {
        match label {
            Label::Positive => &mut self.positive,
            Label::Negative => &mut self.negative,
        }
    }

    fn map<U>(self, mut f: impl FnMut(T) -> U) -> ByLabel<U> {
        ByLabel { positive: f(self.positive), negative: f(self.negative) }
    }
}

/// How many examples of each class a training draw asks for.
pub(super) type SamplingRequest = ByLabel<usize>;

impl SamplingRequest {
    /// What the draw takes from `available` examples: the request, but never
    /// more than [`training_share`] of a class, so a holdout always remains.
    pub(super) fn draw_from(self, available: ByLabel<usize>) -> Self {
        ByLabel {
            positive: self.positive.min(training_share(available.positive)),
            negative: self.negative.min(training_share(available.negative)),
        }
    }
}

/// How many of `available` examples the training draw may take.
fn training_share(available: usize) -> usize {
    if available <= 1 {
        return available;
    }
    (available / 5 * 4).max(1)
}

/// Uniform keys for the reservoir, one stream per segment and model.
///
/// Seeding per segment rather than per run makes the keys a segment draws
/// independent of how rayon interleaved it with the others, and per model
/// keeps one model's sample independent of what the others were offered.
pub(super) struct KeySource(StdRng);

impl KeySource {
    pub(super) fn for_segment(seed: u64, segment: usize, model: MlModel) -> Self {
        let segment = u64::try_from(segment).unwrap_or(u64::MAX);
        let model = u64::try_from(model.index()).unwrap_or(u64::MAX);
        Self(StdRng::from_seed(bytemuck::cast([seed, segment, model, 0].map(u64::to_le))))
    }

    fn next_key(&mut self) -> u64 {
        self.0.random()
    }
}

/// One collected candidate; its features are the row at `offset` in
/// [`TrainingData::features`].
#[derive(Debug, Clone)]
pub(super) struct Example {
    pub(super) label: Label,
    pub(super) chrom: SmolStr,
    pub(super) pos: u64,
    /// Decides which examples a trimmed reservoir keeps: the smallest.
    key: u64,
    offset: usize,
}

impl Example {
    /// The order a reservoir keeps the smallest of, and leaves its survivors in.
    /// Position breaks key ties so the order is total.
    fn by_key(&self, other: &Self) -> Ordering {
        self.key
            .cmp(&other.key)
            .then_with(|| self.chrom.cmp(&other.chrom))
            .then_with(|| self.pos.cmp(&other.pos))
    }
}

/// One model's candidates: a uniform sample per class of everything offered,
/// with their feature rows in one flat buffer.
pub(super) struct TrainingData {
    /// Row-major, `stride` values per example.
    features: Vec<f32>,
    stride: usize,
    examples: Vec<Example>,
    caps: ByLabel<usize>,
    /// Everything ever offered, including what trimming discarded; the holdout
    /// keeps this ratio rather than the reservoir's.
    seen: ByLabel<u64>,
    /// Candidates offered without a usable feature row.
    rejected: u64,
}

impl TrainingData {
    /// A reservoir for rows of `stride` features, with room for `request`
    /// plus a full holdout of each class.
    pub(super) fn for_request(stride: usize, request: SamplingRequest) -> Self {
        Self::bounded(stride, request.map(|n| n.saturating_add(MAX_HOLDOUT)))
    }

    fn bounded(stride: usize, caps: ByLabel<usize>) -> Self {
        Self {
            features: Vec::new(),
            stride,
            examples: Vec::new(),
            caps,
            seen: ByLabel::default(),
            rejected: 0,
        }
    }

    pub(super) fn add_example(
        &mut self,
        row: &[f32],
        label: Label,
        chrom: SmolStr,
        pos: u64,
        keys: &mut KeySource,
    ) -> Result<()> {
        ensure!(
            row.len() == self.stride,
            "Feature row has {} values, this model takes {}",
            row.len(),
            self.stride
        );
        let offset = self.features.len();
        self.features.extend_from_slice(row);
        self.examples.push(Example { label, chrom, pos, key: keys.next_key(), offset });
        *self.seen.get_mut(label) += 1;
        self.trim_if_oversized();
        Ok(())
    }

    pub(super) fn merge(&mut self, other: TrainingData) -> Result<()> {
        ensure!(
            self.stride == other.stride,
            "Cannot merge rows of {} values into rows of {}",
            other.stride,
            self.stride
        );
        let base = self.features.len();
        self.features.extend(other.features);
        self.examples.extend(other.examples.into_iter().map(|mut example| {
            example.offset = example.offset.saturating_add(base);
            example
        }));
        self.seen.positive += other.seen.positive;
        self.seen.negative += other.seen.negative;
        self.rejected += other.rejected;
        self.trim_if_oversized();
        Ok(())
    }

    /// Trim only once the reservoir holds twice its caps, so trimming costs
    /// amortised O(1) per example.
    fn trim_if_oversized(&mut self) {
        let total = self.caps.positive.saturating_add(self.caps.negative);
        if self.examples.len() >= total.saturating_mul(2) {
            self.trim();
        }
    }

    /// Keep, per class, the examples with the smallest keys.
    ///
    /// The keys are uniform, so the kept set is a uniform sample of the class
    /// whatever order the examples arrived in, which is what lets segments be
    /// trimmed separately and merged.
    fn trim(&mut self) {
        let positives = self.positives();
        if positives <= self.caps.positive
            && self.examples.len().saturating_sub(positives) <= self.caps.negative
        {
            return;
        }

        let (mut positive, mut negative): (Vec<Example>, Vec<Example>) =
            std::mem::take(&mut self.examples).into_iter().partition(|e| e.label.is_positive());
        let mut keep = Vec::with_capacity(
            self.caps.positive.min(positive.len()) + self.caps.negative.min(negative.len()),
        );
        for (class, cap) in
            [(&mut positive, self.caps.positive), (&mut negative, self.caps.negative)]
        {
            if let Some(last) = cap.checked_sub(1)
                && class.len() > cap
            {
                class.select_nth_unstable_by(last, Example::by_key);
            }
            class.truncate(cap);
            keep.append(class);
        }
        self.rebuild(keep);
    }

    /// Trim to the caps and sort the survivors by key.
    ///
    /// The kept set does not depend on merge order but its order does, and the
    /// draw shuffles positions within that order.
    pub(super) fn finish(&mut self) {
        self.trim();
        let mut examples = std::mem::take(&mut self.examples);
        examples.sort_unstable_by(Example::by_key);
        self.rebuild(examples);
    }

    /// Rewrite the feature buffer to hold exactly these examples' rows, in
    /// their order.
    fn rebuild(&mut self, examples: Vec<Example>) {
        let old = std::mem::take(&mut self.features);
        let stride = self.stride;
        let mut features = Vec::with_capacity(examples.len().saturating_mul(stride));
        let mut kept = Vec::with_capacity(examples.len());
        for mut example in examples {
            let row =
                example.offset.checked_add(stride).and_then(|end| old.get(example.offset..end));
            let Some(row) = row else {
                warn!(offset = example.offset, stride, "training example has no feature row");
                continue;
            };
            example.offset = features.len();
            features.extend_from_slice(row);
            kept.push(example);
        }
        self.features = features;
        self.examples = kept;
    }

    pub(super) fn row_of(&self, example: &Example) -> Option<&[f32]> {
        self.features.get(example.offset..example.offset.checked_add(self.stride)?)
    }

    pub(super) fn examples(&self) -> &[Example] {
        &self.examples
    }

    pub(super) fn seen(&self) -> ByLabel<u64> {
        self.seen
    }

    pub(super) fn reject(&mut self) {
        self.rejected += 1;
    }

    pub(super) fn rejected(&self) -> u64 {
        self.rejected
    }

    pub(super) fn len(&self) -> usize {
        self.examples.len()
    }

    pub(super) fn is_empty(&self) -> bool {
        self.examples.is_empty()
    }

    pub(super) fn positives(&self) -> usize {
        self.examples.iter().filter(|e| e.label.is_positive()).count()
    }

    /// Retained examples per class.
    pub(super) fn kept(&self) -> ByLabel<usize> {
        let positive = self.positives();
        ByLabel { positive, negative: self.len().saturating_sub(positive) }
    }
}

/// Feature rows widened to `f64`, with their labels' weights.
pub(super) struct Matrix {
    pub(super) features: Array2<f64>,
    pub(super) labels: Array1<f64>,
}

/// The rows a forest is fit on, and the held-out rows its Platt scaling is.
pub(super) struct Split {
    pub(super) train: Matrix,
    pub(super) holdout: Matrix,
}

/// Draw `request` examples for training and a holdout from what is left.
///
/// The holdout keeps the class ratio of everything collection saw: the
/// reservoir is far more balanced than the population, and calibrating on it
/// would shift what an ML threshold means.
pub(super) fn split(data: &TrainingData, request: SamplingRequest, seed: u64) -> Result<Split> {
    let mut rng = StdRng::seed_from_u64(seed);

    let (mut positives, mut negatives): (Vec<usize>, Vec<usize>) =
        (0..data.len()).partition(|&i| data.examples.get(i).is_some_and(|e| e.label.is_positive()));

    let draw = request.draw_from(ByLabel { positive: positives.len(), negative: negatives.len() });
    ensure!(draw.positive > 0, "No positive examples available for training");
    ensure!(draw.negative > 0, "No negative examples available for training");

    positives.shuffle(&mut rng);
    negatives.shuffle(&mut rng);
    let (train_pos, left_pos) = positives.split_at(draw.positive);
    let (train_neg, left_neg) = negatives.split_at(draw.negative);

    let holdout = holdout_shape(data.seen, left_pos.len(), left_neg.len());
    let holdout_indices =
        left_pos.iter().take(holdout.positive).chain(left_neg.iter().take(holdout.negative));

    Ok(Split {
        train: build_matrix(data, train_pos.iter().chain(train_neg).copied().collect())?,
        holdout: build_matrix(data, holdout_indices.copied().collect())?,
    })
}

/// Split a holdout of at most [`MAX_HOLDOUT`] examples between the classes in
/// the proportion collection observed, without exceeding either leftover.
/// A class with leftovers always gets at least one slot: Platt calibration
/// cannot be fit without both.
///
/// Falls back to [`capped_leftovers`] when there is no population to follow or
/// the proportional split rounds to nothing: an empty holdout aborts training.
#[expect(
    clippy::cast_possible_truncation,
    reason = "the f64 to usize casts saturate, and each is clamped by a `min` against a count"
)]
fn holdout_shape(seen: ByLabel<u64>, leftover_pos: usize, leftover_neg: usize) -> ByLabel<usize> {
    let total_seen = seen.positive.saturating_add(seen.negative);
    if total_seen == 0 {
        return capped_leftovers(leftover_pos, leftover_neg);
    }
    let positive_share = seen.positive as f64 / total_seen as f64;
    let negative_share = 1.0 - positive_share;

    // Largest holdout whose class split fits inside both leftovers.
    let mut size = MAX_HOLDOUT.min(leftover_pos.saturating_add(leftover_neg));
    if positive_share > 0.0 {
        size = size.min((leftover_pos as f64 / positive_share).floor() as usize);
    }
    if negative_share > 0.0 {
        size = size.min((leftover_neg as f64 / negative_share).floor() as usize);
    }

    let want_pos = ((size as f64) * positive_share).round() as usize;
    let want_neg = size.saturating_sub(want_pos);
    if want_pos.saturating_add(want_neg) == 0 {
        return capped_leftovers(leftover_pos, leftover_neg);
    }
    ByLabel {
        positive: want_pos.max(1).min(leftover_pos),
        negative: want_neg.max(1).min(leftover_neg),
    }
}

/// Every leftover, together bounded by [`MAX_HOLDOUT`] and split in the
/// leftovers' own proportion.
fn capped_leftovers(leftover_pos: usize, leftover_neg: usize) -> ByLabel<usize> {
    let total = leftover_pos.saturating_add(leftover_neg);
    if total <= MAX_HOLDOUT {
        return ByLabel { positive: leftover_pos, negative: leftover_neg };
    }
    let share = MAX_HOLDOUT as u128 * leftover_pos as u128 / total as u128;
    let want_pos = usize::try_from(share).unwrap_or(MAX_HOLDOUT);
    ByLabel {
        positive: want_pos.min(leftover_pos),
        negative: MAX_HOLDOUT.saturating_sub(want_pos).min(leftover_neg),
    }
}

/// The selected examples' rows, in stored (key) order, widened to `f64`.
fn build_matrix(data: &TrainingData, mut indices: Vec<usize>) -> Result<Matrix> {
    ensure!(
        !indices.is_empty(),
        "No examples left for this matrix: training needs a holdout for Platt calibration. \
         Provide more training data."
    );
    ensure!(data.stride > 0, "Training data carries no features");
    indices.sort_unstable();

    let mut values = Vec::with_capacity(indices.len().saturating_mul(data.stride));
    let mut labels = Vec::with_capacity(indices.len());
    for &index in &indices {
        let example =
            data.examples.get(index).with_context(|| format!("No training example {index}"))?;
        let row = data
            .row_of(example)
            .with_context(|| format!("Training example {index} has no feature row"))?;
        values.extend(row.iter().copied().map(f64::from));
        labels.push(example.label.weight());
    }

    Ok(Matrix {
        features: Array2::from_shape_vec((indices.len(), data.stride), values)
            .wrap_err("Failed to build feature matrix")?,
        labels: Array1::from_vec(labels),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const POS: Label = Label::Positive;
    const NEG: Label = Label::Negative;
    const UNBOUNDED: ByLabel<usize> = ByLabel { positive: usize::MAX, negative: usize::MAX };

    fn caps(positive: usize, negative: usize) -> ByLabel<usize> {
        ByLabel { positive, negative }
    }

    /// Row `i` is `[i, i + 0.5, i + 0.25]`, at position `i`.
    fn data_of(labels: &[Label]) -> TrainingData {
        let mut keys = KeySource::for_segment(0, 0, MlModel::Cpg);
        let mut data = TrainingData::bounded(3, UNBOUNDED);
        for (i, &label) in (0_u16..).zip(labels) {
            let base = f32::from(i);
            data.add_example(
                &[base, base + 0.5, base + 0.25],
                label,
                "chr1".into(),
                i.into(),
                &mut keys,
            )
            .unwrap();
        }
        data
    }

    fn row(data: &TrainingData, index: usize) -> Option<&[f32]> {
        data.row_of(data.examples.get(index)?)
    }

    #[test]
    fn rows_keep_their_boundaries() {
        let data = data_of(&[POS, NEG, NEG]);
        assert_eq!(data.kept(), caps(1, 2));
        assert_eq!(row(&data, 0), Some(&[0.0, 0.5, 0.25][..]));
        assert_eq!(row(&data, 2), Some(&[2.0, 2.5, 2.25][..]));
        assert_eq!(row(&data, 3), None);
    }

    /// A reservoir capped below the draw would discard the very examples the
    /// draw asks for, and one without holdout headroom would starve Platt.
    #[test]
    fn a_reservoir_has_room_for_its_draw_and_a_full_holdout() {
        let data = TrainingData::for_request(3, caps(8_000, 20_000));
        assert_eq!(data.caps, caps(8_000 + MAX_HOLDOUT, 20_000 + MAX_HOLDOUT));
    }

    #[test]
    fn a_row_of_the_wrong_width_is_rejected() {
        let mut data = TrainingData::bounded(3, UNBOUNDED);
        let mut keys = KeySource::for_segment(0, 0, MlModel::Cpg);
        assert!(data.add_example(&[1.0, 2.0], POS, "chr1".into(), 0, &mut keys).is_err());
        assert!(data.is_empty());
    }

    #[test]
    fn rows_of_another_width_are_not_merged() {
        let mut data = TrainingData::bounded(2, UNBOUNDED);
        assert!(data.merge(data_of(&[POS])).is_err());
        assert!(data.is_empty());
    }

    #[test]
    fn merging_appends_rows_without_shifting_them() {
        let mut left = data_of(&[POS, NEG]);
        left.merge(data_of(&[NEG])).unwrap();
        assert_eq!(left.len(), 3);
        assert_eq!(row(&left, 2), Some(&[0.0, 0.5, 0.25][..]));
    }

    #[test]
    fn a_matrix_widens_only_the_selected_rows_in_stored_order() {
        let data = data_of(&[POS, NEG, NEG, POS]);
        let matrix = build_matrix(&data, vec![3, 0]).unwrap();
        assert_eq!(matrix.features.shape(), &[2, 3]);
        assert_eq!(matrix.features.row(0).to_vec(), vec![0.0_f64, 0.5, 0.25]);
        assert_eq!(matrix.features.row(1).to_vec(), vec![3.0_f64, 3.5, 3.25]);
        assert_eq!(matrix.labels.to_vec(), vec![1.0, 1.0]);
    }

    #[test]
    fn an_empty_selection_is_an_error() {
        assert!(build_matrix(&data_of(&[POS]), Vec::new()).is_err());
    }

    #[test]
    fn a_draw_is_capped_by_both_the_request_and_the_pool() {
        let request = SamplingRequest { positive: 8_000, negative: 20_000 };
        assert_eq!(request.draw_from(caps(100_000, 100_000)), request);
        assert_eq!(request.draw_from(caps(10, 5_000)), caps(8, 4_000));
    }

    #[test]
    fn the_training_draw_never_takes_a_whole_class() {
        assert_eq!(training_share(0), 0);
        assert_eq!(training_share(1), 1);
        assert_eq!(training_share(10), 8);
        assert!(training_share(1_000_000) > 20_000);
    }

    #[test]
    fn a_reservoir_never_exceeds_its_class_caps_but_counts_everything() {
        let mut keys = KeySource::for_segment(0, 0, MlModel::Cpg);
        let mut data = TrainingData::bounded(3, caps(3, 4));
        for i in 0..200_u16 {
            let row = [f32::from(i), 0.0, 0.0];
            data.add_example(&row, Label::of(i % 2 == 0), "chr1".into(), i.into(), &mut keys)
                .unwrap();
        }
        data.trim();
        assert!(data.kept().positive <= 3 && data.kept().negative <= 4);
        for example in data.examples() {
            assert_eq!(data.row_of(example).map(<[f32]>::len), Some(3));
        }
        assert_eq!(data.seen(), ByLabel { positive: 100, negative: 100 });
    }

    /// A row must stay with its own example through trimming: each row here
    /// states its position and class, so a row paired with a neighbour's
    /// label shows up as a mismatch.
    #[test]
    fn trimming_keeps_every_row_with_its_own_example() {
        let marker = |label: Label| if label.is_positive() { 1.0_f32 } else { -1.0 };
        let mut keys = KeySource::for_segment(0, 0, MlModel::Cpg);
        let mut data = TrainingData::bounded(3, caps(7, 9));
        for pos in 0..500_u16 {
            let label = Label::of(pos % 3 == 0);
            let row = [f32::from(pos), marker(label), 0.0];
            data.add_example(&row, label, "chr1".into(), pos.into(), &mut keys).unwrap();
        }
        data.trim();

        assert!(data.kept().positive <= 7 && data.kept().negative <= 9);
        for example in data.examples() {
            let row = data.row_of(example).unwrap();
            assert_eq!(row.first().copied().map(f64::from), Some(example.pos as f64));
            assert_eq!(row.get(1).copied(), Some(marker(example.label)));
        }
    }

    #[test]
    fn merged_reservoirs_stay_bounded() {
        let mut left = TrainingData::bounded(3, caps(2, 2));
        let mut right = TrainingData::bounded(3, caps(2, 2));
        for (segment, target) in [&mut left, &mut right].into_iter().enumerate() {
            let mut keys = KeySource::for_segment(0, segment, MlModel::Cpg);
            for i in 0..50_u16 {
                let row = [f32::from(i), 1.0, 2.0];
                target
                    .add_example(&row, Label::of(i % 2 == 0), "chr1".into(), i.into(), &mut keys)
                    .unwrap();
            }
        }
        left.merge(right).unwrap();
        left.trim();
        assert!(left.kept().positive <= 2 && left.kept().negative <= 2);
        assert_eq!(left.seen(), ByLabel { positive: 50, negative: 50 });
    }

    /// `--seed` only reproduces a run if the keys are seeded and the merged
    /// order does not depend on which segment finished first.
    #[test]
    fn the_same_seed_collects_the_same_examples_in_the_same_order() {
        fn segment_of(seed: u64, segment: u16, caps: ByLabel<usize>) -> TrainingData {
            let mut keys = KeySource::for_segment(seed, segment.into(), MlModel::Cpg);
            let mut data = TrainingData::bounded(3, caps);
            for i in 0..40_u16 {
                let pos = segment * 100 + i;
                let row = [f32::from(pos), 0.0, 0.0];
                data.add_example(&row, Label::of(i % 2 == 0), "chr1".into(), pos.into(), &mut keys)
                    .unwrap();
            }
            data
        }

        fn collect(seed: u64, merge_order: [u16; 4], caps: ByLabel<usize>) -> Vec<(u64, u64)> {
            let mut merged = TrainingData::bounded(3, caps);
            for segment in merge_order {
                merged.merge(segment_of(seed, segment, caps)).unwrap();
            }
            merged.finish();
            merged.examples().iter().map(|e| (e.pos, e.key)).collect()
        }

        // Below the caps, where merging is concatenation, and tight enough
        // that the reservoir does the choosing.
        for caps in [caps(1_000, 1_000), caps(5, 5)] {
            let in_order = collect(7, [0, 1, 2, 3], caps);
            assert!(!in_order.is_empty());
            assert_eq!(
                in_order,
                collect(7, [3, 1, 0, 2], caps),
                "merge order mattered at {caps:?}"
            );
            assert_ne!(in_order, collect(8, [0, 1, 2, 3], caps), "seed ignored at {caps:?}");
        }
    }

    #[test]
    fn every_segment_and_model_draws_its_own_keys() {
        let keys = |segment, model| KeySource::for_segment(7, segment, model).next_key();
        let first = keys(0, MlModel::Cpg);
        assert_eq!(first, keys(0, MlModel::Cpg));
        assert_ne!(first, keys(1, MlModel::Cpg));
        assert_ne!(first, keys(0, MlModel::Deletion));
    }

    #[test]
    fn the_holdout_follows_the_population_ratio_not_the_kept_one() {
        let seen = ByLabel { positive: 1_000, negative: 99_000 };
        let holdout = holdout_shape(seen, 100_000, 100_000);
        let total = holdout.positive + holdout.negative;
        let ratio = holdout.positive as f64 / total as f64;
        assert!((ratio - 0.01).abs() < 0.005, "holdout ratio was {ratio}");
        assert!(total <= MAX_HOLDOUT);
    }

    #[test]
    fn the_holdout_never_asks_for_more_than_a_class_has_left() {
        let holdout = holdout_shape(ByLabel { positive: 500, negative: 500 }, 4, 900);
        assert!(holdout.positive <= 4 && holdout.negative <= 900);
        assert!(holdout.positive + holdout.negative > 0);
    }

    #[test]
    fn even_the_holdout_fallback_stays_within_the_bound() {
        let holdout = holdout_shape(ByLabel { positive: 500, negative: 500 }, 0, MAX_HOLDOUT * 5);
        assert_eq!(holdout, caps(0, MAX_HOLDOUT));
    }

    /// A class rarer than one in twice the holdout size would round to no
    /// slot at all, and the calibration would then be fit on one class.
    #[test]
    fn a_rare_class_keeps_a_holdout_slot() {
        let holdout = holdout_shape(ByLabel { positive: 3, negative: 10_000_000 }, 1, 500_000);
        assert_eq!(holdout.positive, 1);
        assert!(holdout.negative > 0);
    }

    #[test]
    fn a_holdout_with_no_population_counts_takes_what_is_there() {
        assert_eq!(holdout_shape(ByLabel::default(), 3, 7), caps(3, 7));
    }

    #[test]
    fn a_small_model_still_leaves_a_holdout() {
        let mut labels = vec![POS; 5];
        labels.extend([NEG; 5]);
        let data = data_of(&labels);
        let split =
            split(&data, SamplingRequest { positive: 8_000, negative: 20_000 }, 42).unwrap();
        assert!(split.train.features.nrows() > 0);
        assert_eq!(split.train.features.nrows(), split.train.labels.len());
        assert!(split.holdout.features.nrows() > 0);
        assert_eq!(split.holdout.features.nrows(), split.holdout.labels.len());
        assert!(split.train.features.nrows() + split.holdout.features.nrows() <= data.len());
    }
}
