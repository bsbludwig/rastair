//! Random Forest based classification of positions
//!
//! This module implements feature extraction for machine learning classification
//! of methylation sites in CpG contexts. The feature extraction logic replicates
//! Ben's R notebook analysis for training random forest models on VCF data.
//!
//! ## CpG Feature Extraction
//!
//! The `cpg` module contains functions to extract features from VCF records
//! including:
//! - Basic variant information (ref, alt, mapping quality, etc.)
//! - Sequence context one-hot encoding
//! - Normalized allele depths and strand bias counts
//! - Base and mapping quality metrics by strand
//! - Adjacent position features for methylation evidence
//!
//! The extracted features are returned as an ndarray Array1<f64> suitable
//! for input to a random forest classifier.

use crate::{
    call::process::{GPU_BATCH_BUFFER_SIZE, InferenceStage},
    metrics::ml::types::{GpuRastairModel, MachineLearning, RastairFlatModel},
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

pub const DEFAULT_ML_THRESHOLD: Probability = Probability::new_panicky(0.5);

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct MachineLearningParams {
    /// Only use hard thresholds to call variants and methylation events.
    ///
    /// This disables using the machine learning models. This will make rastair
    /// much faster, but at the cost of accuracy.
    #[arg(long, help_heading = cli::sections::FILTER)]
    #[default(false)]
    pub no_ml: bool,
    /// Use machine learning model with this threshold value to call variants
    /// and methylation events
    ///
    /// When specified, a ML model will classify positions with a prediction
    /// score. Anything above this threshold is considered PASS.
    ///
    /// For consistency with `--no-ml`, this option can be also be specified as
    /// `--ml` without a value, which will use the default threshold.
    #[arg(long = "ml", default_value_t = DEFAULT_ML_THRESHOLD, default_missing_value = "0.5", num_args = 0..=1)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(DEFAULT_ML_THRESHOLD)]
    pub ml: Probability,
    /// Use GPU-accelerated ML predictions which might speed up large datasets,
    /// but can be slower for small ones due to overhead. The
    /// pair-HMM indel scoring runs on the GPU too.
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
    /// Alts the forest scores at 0.005 or more are scored as events of the
    /// engine's haplotypes, and a second set of forests in the model scores
    /// them again with that evidence. On by default; `--snv-engine=false`
    /// scores with the first forests alone.
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
    pub fn init(&self, workers: usize) -> Result<MachineLearning> {
        if self.no_ml {
            return Ok(MachineLearning::disabled());
        };

        let mut model = load_combined_model(
            self.model.as_ref().map(|x| x.path()),
            &include_bytes!("../../models/rastair_default.rff.mpk.lz4")[..],
        )
        .wrap_err("Failed to load combined RF model")?;
        let engine_forests = model.snv_engine.take();
        let mut ml = self.with_forests(model, workers)?;
        if self.snv_engine {
            let forests = engine_forests.ok_or_else(|| {
                eyre!(
                    "The model has no SNV engine forests; train one with \
                     `rastair ml train --snv-engine`, or pass `--snv-engine=false`"
                )
            })?;
            ml.snv_engine = Some(Box::new(self.with_forests(*forests, workers)?));
        }
        Ok(ml)
    }

    /// One set of forests, ignoring the model's SNV engine forests.
    pub(crate) fn with_forests(
        &self,
        model: RastairFlatModel,
        workers: usize,
    ) -> Result<MachineLearning> {
        let gpu = if self.gpu {
            let ctx = GpuContext::new()
                .wrap_err("Failed to initialise GPU context")
                .note(GpuInitError::hints())?;
            let gpu_forest =
                |forest| GpuForest::with_context(Arc::clone(&ctx), forest, GPU_BATCH_BUFFER_SIZE);

            Some(GpuRastairModel {
                cpg: gpu_forest(&model.cpg),
                denovo: gpu_forest(&model.denovo),
                others: gpu_forest(&model.others),
            })
        } else {
            None
        };

        let inference = gpu
            .map(|gpu| InferenceStage::spawn(gpu, workers))
            .transpose()
            .wrap_err("Failed to start the GPU inference thread")?;

        Ok(MachineLearning {
            threshold: self.ml,
            feature_set: model.feature_set,
            feature_calculator: model.feature_set.get_calculator(),
            model: Some(Box::new(model)),
            inference,
            snv_engine: None,
        })
    }

    pub fn threshold(&self) -> Option<Probability> {
        if self.no_ml { None } else { Some(self.ml) }
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
        rmp_serde::from_read(decompress).wrap_err("Failed to deserialize combined model")
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
    use color_eyre::eyre::ContextCompat as _;

    const BUNDLED: &[u8] = include_bytes!("../../models/rastair_default.rff.mpk.lz4");

    #[test]
    fn the_bundled_model_matches_the_feature_layout() -> Result<()> {
        load_combined_model(None, BUNDLED)?;
        Ok(())
    }

    #[test]
    fn a_model_fitted_on_another_layout_is_refused() -> Result<()> {
        let mut model = load_combined_model(None, BUNDLED)?;
        let layout = model.feature_layout.as_mut().wrap_err("bundled model has a layout")?;
        layout.denovo ^= 1;
        let error = model.check_feature_layout().expect_err("a changed layout is refused");
        let message = format!("{error:#}");
        assert!(message.contains("denovo forest"), "{message}");
        assert!(message.contains("rastair ml train"), "{message}");

        model.feature_layout = None;
        let error = model.check_feature_layout().expect_err("a missing layout is refused");
        assert!(format!("{error:#}").contains("rastair ml train"));
        Ok(())
    }
}
