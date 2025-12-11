use std::fmt;

use biosphere::RandomForest;
use ndarray::Array1;
use rastair_types::{Base, Probability};

/// Combined model file containing all three random forest models
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct RastairModel {
    #[serde(default)]
    pub feature_set: MlFeatureSet,
    pub cpg: RandomForest,
    pub denovo: RandomForest,
    pub others: RandomForest,
}

#[derive(Debug, Clone, Copy, Default, clap::ValueEnum, serde::Serialize, serde::Deserialize)]
pub enum MlFeatureSet {
    #[default]
    Standard,
    Simple,
}

impl fmt::Display for MlFeatureSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MlFeatureSet::Standard => write!(f, "standard"),
            MlFeatureSet::Simple => write!(f, "simple"),
        }
    }
}

/// Instance of machine learning model and parameters
pub struct MachineLearning {
    pub threshold: Probability,
    pub model: Option<Box<RastairModel>>,
    pub feature_calculator: MlFeatureSet,
}

impl MachineLearning {
    /// Create a disabled ML instance
    pub fn disabled() -> Self {
        Self {
            threshold: Probability::ZERO,
            model: None,
            feature_calculator: MlFeatureSet::Standard,
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
