//! The forests `call` scores alts with: the bundled model or `--model`, on
//! the CPU or the GPU, and under `--snv-engine` the second set. Their features
//! are `metrics::ml::features`.
//!
//! The default `--ml` threshold is 0.4, not 0.5: the engine-stage score is
//! under-confident in the middle (at 10x, alts scored .3–.4 are about 52 %
//! true), so 0.4 gains recall at low depth and is neutral at 26x and 73x.

use crate::{
    call::process::{GPU_BATCH_BUFFER_SIZE, InferenceStage},
    metrics::ml::types::{Cascade, ForestSet, GpuRastairModel, MachineLearning, RastairFlatModel},
    utils::cli,
};
use better_default::Default;
use biosphere::gpu::{GpuContext, GpuForest, GpuInitError};
use clap::value_parser;
use clio::ClioPath;
use color_eyre::{
    Result, Section,
    eyre::{Context, ensure, eyre},
};
use seqair_types::Probability;
use std::{fs, io::Read, path::Path, sync::Arc};
use tracing::{debug, instrument};

pub const DEFAULT_ML_THRESHOLD: Probability = Probability::new_panicky(0.4);
/// clap takes the value of a bare `--ml` as a string, which cannot be derived
/// from the constant at compile time; a test keeps the two equal.
const DEFAULT_ML_THRESHOLD_STR: &str = "0.4";

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct MachineLearningParams {
    /// The ML score an alt needs to be called a variant
    ///
    /// The forests score every alt; one at or above this threshold is PASS.
    /// `--ml` without a value uses the default threshold.
    #[arg(long = "ml", default_value_t = DEFAULT_ML_THRESHOLD, default_missing_value = DEFAULT_ML_THRESHOLD_STR, num_args = 0..=1)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(DEFAULT_ML_THRESHOLD)]
    pub ml: Probability,
    /// Use GPU-accelerated ML predictions which might speed up large datasets,
    /// but can be slower for small ones due to overhead. The
    /// pair-HMM, indel and SNV clusters alike, runs on the GPU too.
    ///
    /// Requires a Metal/Vulkan/DX12-capable GPU.
    #[arg(long, help_heading = cli::sections::PROCESSING)]
    #[default(false)]
    pub gpu: bool,
    /// Path to the combined model file containing CpG, denovo, and others models
    ///
    /// Default is the bundled model in the Rastair binary.
    #[arg(long, value_parser=value_parser!(ClioPath).exists().is_file())]
    #[arg(help_heading = cli::sections::FILTER)]
    #[serde(skip)]
    model: Option<ClioPath>,
    /// Re-score the forest's SNV candidates with the pair-HMM engine that
    /// calls indels
    ///
    /// Alts the forest scores at 0.005 or more, the two it scores highest at
    /// a position, are scored as events of the engine's haplotypes, and a
    /// second set of forests in the model scores them again with that
    /// evidence. On by default; `--snv-engine=false`
    /// scores with the first forests alone. The engine genotypes the called
    /// SNVs either way.
    #[arg(long, default_value_t = true)]
    #[arg(action = clap::ArgAction::Set, num_args = 0..=1, require_equals = true, default_missing_value = "true")]
    #[arg(help_heading = cli::sections::FILTER)]
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    #[default(true)]
    pub snv_engine: bool,
}

impl MachineLearningParams {
    /// The forests that score every alt, and under `--snv-engine` the
    /// model's second set, which scores the alts the engine has scored.
    #[instrument(name = "init_ml", skip(self))]
    pub fn init(&self, workers: usize) -> Result<Cascade> {
        let floor = crate::call::pileup::from_seqair::CASCADE_MIN_SCORE;
        ensure!(
            *self.ml >= floor,
            "--ml must be at least {floor}: an SNV is called only with the pair-HMM's genotype, \
             and the pair-HMM scores only alts the model scores that high"
        );
        let RastairFlatModel { forests, snv_engine: engine_forests, .. } = load_combined_model(
            self.model.as_ref().map(|x| x.path()),
            &include_bytes!("../../models/rastair_default.rff.mpk.lz4")[..],
        )
        .wrap_err("Failed to load combined RF model")?;
        // One context for both sets of forests: one device is as fast as several.
        let gpu = self
            .gpu
            .then(|| {
                GpuContext::new()
                    .wrap_err("Failed to initialise GPU context")
                    .note(GpuInitError::hints())
            })
            .transpose()?;
        let forests = Self::with_forests(forests, gpu.as_ref(), workers)?;
        let engine = if self.snv_engine {
            let engine_forests = engine_forests.ok_or_else(|| {
                eyre!(
                    "The model has no SNV engine forests; train one with \
                     `rastair ml train --snv-engine`, or pass `--snv-engine=false`"
                )
            })?;
            Some(Self::with_forests(engine_forests, gpu.as_ref(), workers)?)
        } else {
            None
        };
        Ok(Cascade { forests, engine })
    }

    /// One set of forests, on the GPU of `gpu` when given.
    pub(crate) fn with_forests(
        model: ForestSet,
        gpu: Option<&Arc<GpuContext>>,
        workers: usize,
    ) -> Result<MachineLearning> {
        let gpu = gpu.map(|ctx| {
            let gpu_forest =
                |forest| GpuForest::with_context(Arc::clone(ctx), forest, GPU_BATCH_BUFFER_SIZE);
            GpuRastairModel {
                cpg: gpu_forest(&model.cpg),
                denovo: gpu_forest(&model.denovo),
                others: gpu_forest(&model.others),
            }
        });

        let inference = gpu
            .map(|gpu| InferenceStage::spawn(gpu, workers))
            .transpose()
            .wrap_err("Failed to start the GPU inference thread")?;

        Ok(MachineLearning { model: Box::new(model), inference })
    }
}

#[instrument(level = "debug", skip_all)]
pub fn load_combined_model(path: Option<&Path>, built_in: &[u8]) -> Result<RastairFlatModel> {
    fn load_model_from_file(path: &Path) -> Result<RastairFlatModel> {
        ensure!(path.exists(), "Model file does not exist");
        let file = fs::read(path).wrap_err("Failed to read model file")?;
        load_combined(&file[..]).wrap_err("Failed to load combined model")
    }

    fn load_combined(reader: impl Read) -> Result<RastairFlatModel> {
        let decompress = lz4::Decoder::new(reader).wrap_err("Failed to create LZ4 decoder")?;
        rmp_serde::from_read(decompress).wrap_err(
            "Failed to deserialize combined model; a model file of another rastair version \
             cannot be read, retrain it with `rastair ml train`",
        )
    }

    let model = if let Some(path) = path {
        let model = load_model_from_file(path)
            .wrap_err_with(|| format!("Failed to load combined model from {path:?}"))?;
        debug!(?path, "Loaded combined model from file");
        model.check_feature_layout().wrap_err_with(|| format!("Cannot use model {path:?}"))?;
        model
    } else {
        let model = load_combined(built_in)?;
        debug!("Loaded built-in combined model");
        model.check_feature_layout().wrap_err("Cannot use the bundled model")?;
        model
    };
    Ok(model)
}

#[cfg(test)]
mod tests {
    use super::*;

    const BUNDLED: &[u8] = include_bytes!("../../models/rastair_default.rff.mpk.lz4");

    #[test]
    fn a_bare_ml_flag_means_the_default_threshold() -> Result<()> {
        assert_eq!(DEFAULT_ML_THRESHOLD_STR.parse::<Probability>()?, DEFAULT_ML_THRESHOLD);
        Ok(())
    }

    #[test]
    fn the_bundled_model_matches_the_feature_layout() -> Result<()> {
        load_combined_model(None, BUNDLED)?;
        Ok(())
    }

    #[test]
    fn a_model_fitted_on_another_layout_is_refused() -> Result<()> {
        let mut model = load_combined_model(None, BUNDLED)?;
        model.feature_layout.denovo ^= 1;
        let error = model.check_feature_layout().expect_err("a changed layout is refused");
        let message = format!("{error:#}");
        assert!(message.contains("denovo forest"), "{message}");
        assert!(message.contains("rastair ml train"), "{message}");
        Ok(())
    }
}
