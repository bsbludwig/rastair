use super::types::{MachineLearning, MlModel, Prediction};
use crate::{
    metrics::{MetricsForAlt, PileupMetrics},
    vcf::InCpG,
};
use rastair_types::{Base, Probability};
use tracing::{instrument, warn};

impl MachineLearning {
    #[instrument(level = "debug", skip_all)]
    #[allow(clippy::unwrap_in_result)] // it's fine
    pub fn predict2(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Option<Prediction> {
        if self.disabled {
            return None;
        }

        let pos_metrics = &current.metrics.pos_metrics;
        let alt = current.alt;

        let (name, model, features) = if (pos_metrics.cpg == InCpG::C && alt.base == Base::T)
            || (pos_metrics.cpg == InCpG::G && current.alt.base == Base::A)
        {
            (MlModel::Cpg, self.cpg.as_ref(), super::cpg(current, before, after))
        } else if *alt.denovo {
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
        let features = match features {
            Err(error) => {
                warn!(%error, "Failed to generate features for ML prediction");
                return None;
            }
            Ok(x) => x,
        };
        let prediction = model.predict(&features.view());

        match prediction.get(0).copied() {
            Some(p) => Some(Prediction {
                prediction: Probability::new(p).expect("Probability should be valid"),
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
