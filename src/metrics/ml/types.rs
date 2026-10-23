use crate::metrics::ml::features::FeatureCalculatorBox;
use biosphere::FlatForest;
use biosphere::gpu::GpuForest;
use color_eyre::eyre::WrapErr as _;
use seqair_types::Probability;
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

/// The field order is the on-disk layout: msgpack writes a struct as an
/// array, so it is the order of the bundled `rastair_default.rff.mpk.lz4`.
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct RastairFlatModel {
    pub cpg: FlatForest,
    pub denovo: FlatForest,
    pub others: FlatForest,
    #[serde(default)]
    pub cpg_platt: PlattScaling,
    #[serde(default)]
    pub denovo_platt: PlattScaling,
    #[serde(default)]
    pub others_platt: PlattScaling,
    #[serde(default)]
    pub feature_set: MlFeatureSet,
    /// `None` in a model written before layouts were recorded.
    #[serde(default)]
    pub feature_layout: Option<FeatureLayout>,
    /// The forests `--snv-engine` re-scores the engine's alts with, fitted on
    /// the alts the forests above keep, with the engine's features from their
    /// clusters. `None` in a model trained without the engine.
    #[serde(default)]
    pub snv_engine: Option<Box<RastairFlatModel>>,
}

/// A fingerprint of each forest's feature names, in layout order.
///
/// A forest reads a row by position, so one fitted on another layout would
/// score the wrong features without failing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct FeatureLayout {
    pub cpg: u64,
    pub denovo: u64,
    pub others: u64,
}

impl FeatureLayout {
    pub fn of(feature_set: MlFeatureSet) -> Self {
        let names = feature_set.get_calculator().feature_names();
        Self {
            cpg: fingerprint(&names.cpg),
            denovo: fingerprint(&names.denovo_cpg),
            others: fingerprint(&names.others),
        }
    }

    fn get(self, model: MlModel) -> u64 {
        match model {
            MlModel::Others => self.others,
            MlModel::Cpg => self.cpg,
            MlModel::DenovoCpg => self.denovo,
        }
    }
}

/// 64-bit FNV-1a over the names, each terminated by a NUL: stable across
/// builds and platforms, unlike `std`'s hashers.
fn fingerprint(names: &[&str]) -> u64 {
    names
        .iter()
        .flat_map(|name| name.bytes().chain([0]))
        .fold(0xcbf2_9ce4_8422_2325, |hash, byte| {
            (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3)
        })
}

/// One value per [`MlModel`], addressed by the model itself.
///
/// The three models are always handled together — a feature block, a dispatch
/// and a score vector each — so they are one array rather than three parallel
/// fields.
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

    /// Map every entry in [`MlModel::ALL`] order, stopping at the first error.
    pub(crate) fn try_map<U, E>(
        self,
        mut f: impl FnMut(MlModel, T) -> Result<U, E>,
    ) -> Result<ByModel<U>, E> {
        let [a, b, c] = self.0;
        let [ma, mb, mc] = MlModel::ALL;
        Ok(ByModel([f(ma, a)?, f(mb, b)?, f(mc, c)?]))
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
    pub(crate) fn from_trained(
        feature_set: MlFeatureSet,
        trained: ByModel<(FlatForest, PlattScaling)>,
    ) -> Self {
        const {
            assert!(matches!(MlModel::ALL, [MlModel::Others, MlModel::Cpg, MlModel::DenovoCpg]));
        }
        let [(others, others_platt), (cpg, cpg_platt), (denovo, denovo_platt)] = trained.0;
        Self {
            cpg,
            denovo,
            others,
            cpg_platt,
            denovo_platt,
            others_platt,
            feature_set,
            feature_layout: Some(FeatureLayout::of(feature_set)),
            snv_engine: None,
        }
    }

    /// Refuse a model whose forests were not fitted on the features this
    /// build computes.
    pub fn check_feature_layout(&self) -> color_eyre::Result<()> {
        let Some(recorded) = self.feature_layout else {
            color_eyre::eyre::bail!(
                "The model records no feature layout, so it predates the current features; \
                 retrain it with `rastair ml train`"
            );
        };
        let current = FeatureLayout::of(self.feature_set);
        let differing: Vec<&str> = MlModel::ALL
            .into_iter()
            .filter(|&model| recorded.get(model) != current.get(model))
            .map(MlModel::name)
            .collect();
        color_eyre::eyre::ensure!(
            differing.is_empty(),
            "The model's {} forest(s) were trained on a different feature layout than this \
             rastair computes; retrain it with `rastair ml train`",
            differing.join(", ")
        );
        if let Some(engine) = &self.snv_engine {
            engine.check_feature_layout().wrap_err("The model's SNV engine forests")?;
        }
        Ok(())
    }

    pub fn forest(&self, model: MlModel) -> &FlatForest {
        match model {
            MlModel::Others => &self.others,
            MlModel::Cpg => &self.cpg,
            MlModel::DenovoCpg => &self.denovo,
        }
    }

    pub fn platt(&self, model: MlModel) -> PlattScaling {
        match model {
            MlModel::Others => self.others_platt,
            MlModel::Cpg => self.cpg_platt,
            MlModel::DenovoCpg => self.denovo_platt,
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

/// The three forests on one shared [`biosphere::gpu::GpuContext`], owned by the
/// inference thread.
///
/// One queue pipelines the dispatches, which separate devices would
/// time-slice, and on Metal one device measures the same as five
/// (chr12:20-30Mb, 1.87 s vs 1.90 s). Sharing relies on biosphere uploading
/// with `write_buffer`: an upload that waited on the device's most recent
/// submission would wait on the previous model's dispatch.
pub struct GpuRastairModel {
    pub cpg: GpuForest,
    pub denovo: GpuForest,
    pub others: GpuForest,
}

impl GpuRastairModel {
    pub fn forest(&self, model: MlModel) -> &GpuForest {
        match model {
            MlModel::Others => &self.others,
            MlModel::Cpg => &self.cpg,
            MlModel::DenovoCpg => &self.denovo,
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
    /// The forests that re-score the alts the SNV engine scored, under
    /// `--snv-engine`.
    pub snv_engine: Option<Box<MachineLearning>>,
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
            snv_engine: None,
        }
    }

    pub fn enabled(&self) -> bool {
        self.model.is_some()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum MlModel {
    Others,
    Cpg,
    DenovoCpg,
}

impl MlModel {
    /// Every model, in the order [`MlModel::index`] assigns.
    pub const ALL: [Self; 3] = [Self::Others, Self::Cpg, Self::DenovoCpg];

    pub const COUNT: usize = Self::ALL.len();

    /// Position in a per-model array. Only meaningful within one run — nothing
    /// on disk is keyed by it.
    pub const fn index(self) -> usize {
        match self {
            Self::Others => 0,
            Self::Cpg => 1,
            Self::DenovoCpg => 2,
        }
    }

    /// What the model's exports and messages call it.
    pub const fn name(self) -> &'static str {
        match self {
            Self::Others => "other",
            Self::Cpg => "cpg",
            Self::DenovoCpg => "denovo",
        }
    }
}
