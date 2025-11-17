use biosphere::RandomForest;
use ndarray::Array1;
use rastair_types::{Base, Probability};

pub struct MachineLearning {
    pub disabled: bool,
    pub threshold: Probability,
    pub cpg: Option<Box<RandomForest>>,
    pub denovo_cpg: Option<Box<RandomForest>>,
    pub others: Option<Box<RandomForest>>,
}

impl MachineLearning {
    /// Create a disabled ML instance
    pub fn disabled() -> Self {
        Self {
            disabled: true,
            threshold: Probability::new(1.).expect("1 is a valid probability"),
            cpg: None,
            denovo_cpg: None,
            others: None,
        }
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
