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
    pub model: MlModel,
    pub allele: Base,
    pub prediction: Probability,
    pub threshold: Probability,
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
