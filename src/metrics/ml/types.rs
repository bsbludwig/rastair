use crate::metrics::ml::features::FeatureCalculatorBox;
use biosphere::FlatForest;
use biosphere::gpu::GpuForest;
use ndarray::Array1;
use seqair_types::{Probability, SmolStr};
use std::fmt;
use std::ops::{Index, IndexMut};

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

/// One value per [`MlModel`], addressed by the model itself.
///
/// The five models are always handled together — a feature block, a dispatch
/// and a score vector each — so this replaces the five parallel fields that
/// shape would otherwise grow.
#[derive(Debug, Clone)]
pub struct ByModel<T>([T; MlModel::COUNT]);

impl<T> ByModel<T> {
    /// One entry per model, built in [`MlModel::ALL`] order.
    pub fn from_fn(f: impl FnMut(MlModel) -> T) -> Self {
        Self(MlModel::ALL.map(f))
    }

    pub fn iter(&self) -> impl Iterator<Item = (MlModel, &T)> {
        MlModel::ALL.into_iter().zip(self.0.iter())
    }
}

impl<T> IntoIterator for ByModel<T> {
    type Item = (MlModel, T);
    type IntoIter = std::iter::Zip<
        std::array::IntoIter<MlModel, { MlModel::COUNT }>,
        std::array::IntoIter<T, { MlModel::COUNT }>,
    >;

    fn into_iter(self) -> Self::IntoIter {
        MlModel::ALL.into_iter().zip(self.0)
    }
}

// `MlModel::index` is a total map onto `0..COUNT` and the array has exactly
// that many slots, so neither impl can be out of bounds.
impl<T> Index<MlModel> for ByModel<T> {
    type Output = T;

    #[expect(clippy::indexing_slicing, reason = "MlModel::index is total over the array")]
    fn index(&self, model: MlModel) -> &T {
        &self.0[model.index()]
    }
}

impl<T> IndexMut<MlModel> for ByModel<T> {
    #[expect(clippy::indexing_slicing, reason = "MlModel::index is total over the array")]
    fn index_mut(&mut self, model: MlModel) -> &mut T {
        &mut self.0[model.index()]
    }
}

impl RastairFlatModel {
    pub fn forest(&self, model: MlModel) -> &FlatForest {
        match model {
            MlModel::Others => &self.others,
            MlModel::Cpg => &self.cpg,
            MlModel::DenovoCpg => &self.denovo,
            MlModel::Insertion => &self.insertion,
            MlModel::Deletion => &self.deletion,
        }
    }

    pub fn platt(&self, model: MlModel) -> PlattScaling {
        match model {
            MlModel::Others => self.others_platt,
            MlModel::Cpg => self.cpg_platt,
            MlModel::DenovoCpg => self.denovo_platt,
            MlModel::Insertion => self.insertion_platt,
            MlModel::Deletion => self.deletion_platt,
        }
    }
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

/// The five forests on one shared [`biosphere::gpu::GpuContext`], owned by the
/// inference thread.
///
/// They used to sit on five devices, because biosphere's UMA upload waited on
/// the device's most recent submission from any source, which on a shared
/// device was the previous model's dispatch. Biosphere now uploads with
/// `write_buffer` everywhere, and with that fixed one device measures the same
/// as five on Metal (chr12:20-30Mb, 1.87 s vs 1.90 s) while being what a
/// discrete GPU wants: one queue pipelines the five dispatches that separate
/// devices would time-slice.
pub struct GpuRastairModel {
    pub cpg: GpuForest,
    pub denovo: GpuForest,
    pub others: GpuForest,
    pub insertion: GpuForest,
    pub deletion: GpuForest,
}

impl GpuRastairModel {
    pub fn forest(&self, model: MlModel) -> &GpuForest {
        match model {
            MlModel::Others => &self.others,
            MlModel::Cpg => &self.cpg,
            MlModel::DenovoCpg => &self.denovo,
            MlModel::Insertion => &self.insertion,
            MlModel::Deletion => &self.deletion,
        }
    }
}

/// Instance of machine learning model and parameters
pub struct MachineLearning {
    pub threshold: Probability,
    pub model: Option<Box<RastairFlatModel>>,
    pub feature_set: MlFeatureSet,
    pub feature_calculator: FeatureCalculatorBox,
    /// The thread that owns the GPU forests, when `--gpu` is on. Workers hand
    /// it feature rows rather than each forking a copy of the forests.
    pub inference: Option<crate::call::process::InferenceStage>,
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
            inference: None,
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

impl MlModel {
    /// Every model, in the order [`MlModel::index`] assigns.
    pub const ALL: [Self; 5] =
        [Self::Others, Self::Cpg, Self::DenovoCpg, Self::Insertion, Self::Deletion];

    pub const COUNT: usize = Self::ALL.len();

    /// Position in a per-model array. Only meaningful within one run — nothing
    /// on disk is keyed by it.
    pub const fn index(self) -> usize {
        match self {
            Self::Others => 0,
            Self::Cpg => 1,
            Self::DenovoCpg => 2,
            Self::Insertion => 3,
            Self::Deletion => 4,
        }
    }
}
