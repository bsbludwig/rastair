use crate::metrics::ml::features::feature_names;
use biosphere::FlatForest;
use biosphere::gpu::GpuForest;
use seqair_types::Probability;
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

/// A model file. msgpack writes a struct as an array, so the field order is
/// the on-disk layout, here and in [`ForestSet`].
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct RastairFlatModel {
    /// The layout both forest sets were fitted on: one feature layout serves
    /// both stages of the cascade.
    pub feature_layout: FeatureLayout,
    /// The forests that score every alt.
    pub forests: ForestSet,
    /// The forests `--snv-engine` re-scores the engine's alts with, fitted on
    /// the alts [`Self::forests`] keep, with the engine's features from their
    /// clusters. `None` in a model trained without `--snv-engine`.
    pub snv_engine: Option<ForestSet>,
}

/// One forest per [`MlModel`], each with the Platt scaling fitted on its
/// holdout.
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct ForestSet {
    pub cpg: FlatForest,
    pub denovo: FlatForest,
    pub others: FlatForest,
    pub cpg_platt: PlattScaling,
    pub denovo_platt: PlattScaling,
    pub others_platt: PlattScaling,
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
    /// This build's layout.
    pub fn current() -> Self {
        Self {
            cpg: fingerprint(&feature_names(MlModel::Cpg)),
            denovo: fingerprint(&feature_names(MlModel::DenovoCpg)),
            others: fingerprint(&feature_names(MlModel::Others)),
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

impl ForestSet {
    pub(crate) fn from_trained(trained: ByModel<(FlatForest, PlattScaling)>) -> Self {
        const {
            assert!(matches!(MlModel::ALL, [MlModel::Others, MlModel::Cpg, MlModel::DenovoCpg]));
        }
        let [(others, others_platt), (cpg, cpg_platt), (denovo, denovo_platt)] = trained.0;
        Self { cpg, denovo, others, cpg_platt, denovo_platt, others_platt }
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

impl RastairFlatModel {
    /// Refuse a model whose forests were not fitted on the features this
    /// build computes.
    pub fn check_feature_layout(&self) -> color_eyre::Result<()> {
        let current = FeatureLayout::current();
        let differing: Vec<&str> = MlModel::ALL
            .into_iter()
            .filter(|&model| self.feature_layout.get(model) != current.get(model))
            .map(MlModel::name)
            .collect();
        color_eyre::eyre::ensure!(
            differing.is_empty(),
            "The model's {} forest(s) were trained on a different feature layout than this \
             rastair computes; retrain it with `rastair ml train`",
            differing.join(", ")
        );
        Ok(())
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
    pub model: Box<ForestSet>,
    /// The thread that owns the GPU forests, when `--gpu` is on. Workers hand
    /// it feature rows rather than each forking a copy of the forests.
    pub inference: Option<crate::call::process::InferenceStage>,
}

/// The forests `call` scores with: `forests` score every alt, and under
/// `--snv-engine` `engine` re-scores the alts the SNV engine scored.
pub struct Cascade {
    pub forests: MachineLearning,
    pub engine: Option<MachineLearning>,
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
