use crate::metrics::ml::features::FeatureCalculatorBox;
use biosphere::FlatForest;
use biosphere::gpu::GpuForest;
use ndarray::Array1;
use rastair_types::{Base, Probability, SmolStr};
use std::fmt;

#[derive(Clone, Copy, Debug, serde::Serialize, serde::Deserialize)]
pub struct PlattScaling {
    pub a: f64,
    pub b: f64,
}

impl Default for PlattScaling {
    fn default() -> Self {
        Self { a: 1.0, b: 0.0 }
    }
}

impl PlattScaling {
    pub fn calibrate_score(&self, score: f64) -> Probability {
        let z = self.a * score + self.b;
        let p = if z >= 0.0 {
            let ez = (-z).exp();
            ez / (1.0 + ez)
        } else {
            let ez = z.exp();
            1.0 / (1.0 + ez)
        };
        Probability::new(p).unwrap_or(Probability::ZERO)
    }
}

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct RastairFlatModel {
    pub cpg: FlatForest,
    #[serde(default)]
    pub cpg_platt: PlattScaling,
    pub denovo: FlatForest,
    #[serde(default)]
    pub denovo_platt: PlattScaling,
    pub others: FlatForest,
    #[serde(default)]
    pub others_platt: PlattScaling,
    pub insertion: FlatForest,
    #[serde(default)]
    pub insertion_platt: PlattScaling,
    pub deletion: FlatForest,
    #[serde(default)]
    pub deletion_platt: PlattScaling,
    #[serde(default)]
    pub feature_set: MlFeatureSet,
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

/// Flat (BFS-linearised, f32) forests for CPU inference.
///
/// Built once from the loaded [`RastairModel`] and stored in [`MachineLearning`].
/// Uses the same f32 representation as the GPU path so both paths are numerically consistent.
pub struct FlatRastairModel {
    pub cpg: FlatForest,
    pub denovo: FlatForest,
    pub others: FlatForest,
    pub insertion: FlatForest,
    pub deletion: FlatForest,
}

/// GPU-accelerated forests for each model type, used as per-thread prototypes.
///
/// Create once via [`MachineLearningParams::init`], then call [`GpuRastairModel::fork`]
/// inside each worker thread to get a thread-local handle that shares compiled
/// pipelines and uploaded node data without re-uploading.
pub struct GpuRastairModel {
    pub cpg: GpuForest,
    pub denovo: GpuForest,
    pub others: GpuForest,
    pub insertion: GpuForest,
    pub deletion: GpuForest,
}

impl GpuRastairModel {
    /// Create per-thread handles that share GPU pipelines and node data with `self`.
    pub fn fork(&self, max_samples: usize) -> Self {
        Self {
            cpg: self.cpg.fork(max_samples),
            denovo: self.denovo.fork(max_samples),
            others: self.others.fork(max_samples),
            insertion: self.insertion.fork(max_samples),
            deletion: self.deletion.fork(max_samples),
        }
    }
}

/// Instance of machine learning model and parameters
pub struct MachineLearning {
    pub threshold: Probability,
    pub model: Option<Box<RastairFlatModel>>,
    pub feature_set: MlFeatureSet,
    pub feature_calculator: FeatureCalculatorBox,
    /// Prototype GPU forests. Worker threads call [`GpuRastairModel::fork`] on
    /// first use to obtain thread-local handles without recompiling shaders.
    pub gpu_prototype: Option<GpuRastairModel>,
}

impl MachineLearning {
    /// Create a disabled ML instance
    pub fn disabled() -> Self {
        let feature_set = MlFeatureSet::Standard;
        Self {
            threshold: Probability::ZERO,
            model: None,
            feature_set,
            feature_calculator: feature_set.get_calculator(),
            gpu_prototype: None,
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
    pub allele: SmolStr,
    /// Probability of the alt being a variant
    pub prediction: Probability,
    /// Threshold for calling a variant
    pub threshold: Probability,
    /// Features used for this prediction
    #[serde(skip)]
    pub features: Array1<f32>,
}

impl Prediction {
    pub fn pass(&self) -> bool {
        self.prediction >= self.threshold
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum MlModel {
    Others,
    Cpg,
    DenovoCpg,
    Insertion,
    Deletion,
}
