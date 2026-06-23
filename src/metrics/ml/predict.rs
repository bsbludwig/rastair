use super::types::{MachineLearning, MlModel, Prediction};
use crate::metrics::{MetricsForAlt, MetricsForIndel, PileupMetrics};
use color_eyre::eyre::ensure;
use tracing::{debug, instrument, warn};

impl MachineLearning {
    #[instrument(level = "debug", skip_all)]
    #[allow(clippy::unwrap_in_result, reason = "it's fine")]
    pub fn predict(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Option<Prediction> {
        if !self.enabled() {
            return None;
        }

        let Some(model) = self.model.as_ref() else {
            warn!("No model available");
            return None;
        };

        let calc = &self.feature_calculator;

        let (name, flat_forest, platt, features) = if current.is_evidence_for_methylation() {
            (MlModel::Cpg, &model.cpg, &model.cpg_platt, calc.calculate_cpg(current, before, after))
        } else if *current.alt.denovo {
            (
                MlModel::DenovoCpg,
                &model.denovo,
                &model.denovo_platt,
                calc.calculate_denovo_cpg(current, before, after),
            )
        } else {
            (
                MlModel::Others,
                &model.others,
                &model.others_platt,
                calc.calculate_others(current, before, after),
            )
        };
        let features = features.and_then(|x| {
            ensure!(!x.is_any_nan(), "Failed to calculate features (one metric was NaN)");
            Ok(x)
        });
        let features = match features {
            Err(error) => {
                debug!(%error, "Failed to generate features for ML prediction");
                return None;
            }
            Ok(x) => x,
        };
        let prediction = flat_forest.predict(&features.view());

        match prediction.get(0).copied() {
            Some(p) => Some(Prediction {
                prediction: platt.calibrate_score(p),
                threshold: self.threshold,
                allele: current.alt.base.into(),
                features: features.row(0).to_owned(),
                model: name,
            }),
            None => {
                warn!(model=?name, "No predictions");
                None
            }
        }
    }

    #[instrument(level = "debug", skip_all)]
    #[allow(clippy::unwrap_in_result, reason = "it's fine")]
    pub fn predict_indels(&self, current: &MetricsForIndel) -> Option<Prediction> {
        if !self.enabled() {
            return None;
        }

        let Some(model) = self.model.as_ref() else {
            warn!("No model available");
            return None;
        };

        let calc = &self.feature_calculator;

        let (name, flat_forest, platt, features) = if current.indel.allele.is_insertion() {
            (
                MlModel::Insertion,
                &model.insertion,
                &model.insertion_platt,
                calc.calculate_insertion(current),
            )
        } else {
            (
                MlModel::Deletion,
                &model.deletion,
                &model.deletion_platt,
                calc.calculate_deletion(current),
            )
        };

        let features = features.and_then(|x| {
            ensure!(!x.is_any_nan(), "Failed to calculate features (one metric was NaN)");
            Ok(x)
        });
        let features = match features {
            Err(error) => {
                debug!(%error, "Failed to generate features for ML prediction");
                return None;
            }
            Ok(x) => x,
        };
        let prediction = flat_forest.predict(&features.view());

        match prediction.get(0).copied() {
            Some(p) => Some(Prediction {
                prediction: platt.calibrate_score(p),
                threshold: self.threshold,
                allele: current.indel.allele.bases().iter().map(|b| b.as_str()).collect(),
                features: features.row(0).to_owned(),
                model: name,
            }),
            None => {
                warn!(model=?name, "No predictions");
                None
            }
        }
    }
}
