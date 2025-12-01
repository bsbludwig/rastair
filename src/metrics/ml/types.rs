use super::features::FeatureCalculator;
use biosphere::RandomForest;
use ndarray::Array1;
use rastair_types::{Base, Probability};

/// Combined model file containing all three random forest models
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct RastairModel {
    pub cpg: RandomForest,
    pub denovo: RandomForest,
    pub others: RandomForest,
}

pub struct MachineLearning {
    pub threshold: Probability,
    pub model: Option<Box<RastairModel>>,
    pub feature_calculator: Box<dyn FeatureCalculator>,
}

impl MachineLearning {
    /// Create a disabled ML instance
    pub fn disabled() -> Self {
        Self {
            threshold: Probability::new(1.).expect("1 is a valid probability"),
            model: None,
            feature_calculator: Box::new(super::StandardFeatures),
        }
    }

    pub fn enabled(&self) -> bool {
        self.model.is_some()
    }
}

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct Prediction {
    /// The model used for this prediction
    pub model: MlModel,
    /// The alt base this prediction is for
    pub allele: Base,
    /// Probability of the alt being a variant
    pub prediction: Probability,
    /// Threshold for calling a variant
    pub threshold: Probability,
    /// Features used for this prediction
    #[serde(skip)]
    pub features: Array1<f64>,
}

impl Prediction {
    pub fn pass(&self) -> bool {
        self.prediction >= self.threshold
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum MlModel {
    Cpg,
    DenovoCpg,
    Others,
}
