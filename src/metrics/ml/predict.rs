use super::types::{MachineLearning, MlModel, Prediction};
use crate::metrics::{MetricsForAlt, PileupMetrics};
use color_eyre::eyre::ensure;
use rastair_types::Probability;
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
        if self.disabled {
            return None;
        }

        let (name, model, features) = if current.is_evidence_for_methylation() {
            (MlModel::Cpg, self.cpg.as_ref(), super::cpg(current, before, after))
        } else if *current.alt.denovo {
            (
                MlModel::DenovoCpg,
                self.denovo_cpg.as_ref(),
                super::denovo_cpg(current, before, after),
            )
        } else {
            (MlModel::Others, self.others.as_ref(), super::others(current, before, after))
        };

        let Some(model) = model else {
            warn!(model=?name, "No model found");
            return None;
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
                prediction: Probability::new(p)
                    .expect("Got invalid probability value in prediction"),
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
