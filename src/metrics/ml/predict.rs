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
        let Some(flat_model) = self.flat_model.as_ref() else {
            warn!("No flat model available");
            return None;
        };

        let calc = &self.feature_calculator;

        let (name, flat_forest, platt, features) = if current.is_evidence_for_methylation() {
            (
                MlModel::Cpg,
                &flat_model.cpg,
                &rastair_model.cpg_platt,
                calc.calculate_cpg(current, before, after),
            )
        } else if *current.alt.denovo {
            (
                MlModel::DenovoCpg,
                &flat_model.denovo,
                &rastair_model.denovo_platt,
                calc.calculate_denovo_cpg(current, before, after),
            )
        } else {
            (
                MlModel::Others,
                &flat_model.others,
                &rastair_model.others_platt,
                calc.calculate_others(current, before, after),
            )
        };
        let features = features.and_then(|x| {
            ensure!(!x.is_any_nan(), "Failed to calculate features (one metric was NaN)");
            Ok(x)
        });
        let features_f64 = match features {
            Err(error) => {
                debug!(%error, "Failed to generate features for ML prediction");
                return None;
            }
            Ok(x) => x,
        };
        let features_f32 = features_f64.mapv(|x| x as f32);
        let prediction = flat_forest.predict(&features_f32.view());

        match prediction.get(0).copied() {
            Some(p) => Some(Prediction {
                prediction: platt.calibrate_score(p as f64),
                threshold: self.threshold,
                allele: current.alt.base,
                features: features_f32.row(0).to_owned(),
                model: name,
            }),
            None => {
                warn!(model=?name, "No predictions");
                None
            }
        }
    }
}
