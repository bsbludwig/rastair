use super::types::{MachineLearning, MlModel, Prediction};
use crate::metrics::{MetricsForAlt, PileupMetrics};
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

        let Some(rastair_model) = self.model.as_ref() else {
            warn!("No model available");
            return None;
        };

        let feature_calculator = self.feature_calculator.get_calculator();

        let (name, model, platt, features) = if current.is_evidence_for_methylation() {
            (
                MlModel::Cpg,
                &rastair_model.cpg,
                &rastair_model.cpg_platt,
                feature_calculator.calculate_cpg(current, before, after),
            )
        } else if *current.alt.denovo {
            (
                MlModel::DenovoCpg,
                &rastair_model.denovo,
                &rastair_model.denovo_platt,
                feature_calculator.calculate_denovo_cpg(current, before, after),
            )
        } else {
            (
                MlModel::Others,
                &rastair_model.others,
                &rastair_model.others_platt,
                feature_calculator.calculate_others(current, before, after),
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
        let prediction = model.predict(&features.view());

        match prediction.get(0).copied() {
            Some(p) => Some(Prediction {
                prediction: platt.calibrate_score(p),
                threshold: self.threshold,
                allele: current.alt.base,
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
