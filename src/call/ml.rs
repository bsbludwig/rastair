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
    metrics::ml::types::{GpuRastairModel, MachineLearning, RastairModel},
    utils::cli,
};
use better_default::Default;
use biosphere::{FlatForest, gpu::GpuForest};
use clap::value_parser;
use clio::ClioPath;
use color_eyre::{
    Result,
    eyre::{Context, ensure},
};
use rastair_types::Probability;
use std::{fs, io::Read, path::Path};
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
    /// Use GPU-accelerated batch inference for ML predictions.
    ///
    /// Instead of running one random forest prediction per alt allele (CPU,
    /// sequential), batches all alts in a chunk into a single GPU dispatch per
    /// model type. Requires a Metal/Vulkan/DX12-capable GPU.
    ///
    /// Off by default; use this flag to benchmark GPU vs CPU throughput.
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
}

impl MachineLearningParams {
    #[instrument(name = "init_ml", skip(self))]
    pub fn init(&self) -> Result<MachineLearning> {
        if self.no_ml {
            return Ok(MachineLearning::disabled());
        };

        let combined = load_combined_model(
            self.model.as_ref().map(|x| x.path()),
            &include_bytes!("../../models/rastair_default.rf.mpk.lz4")[..],
        )
        .wrap_err("Failed to load combined RF model")?;

        let feature_nums = combined.feature_set.get_calculator().feature_num();
        // Size buffers for a full chunk (10k positions × 4 alts max) per thread.
        let max_samples = 40_000;
        let gpu_prototype = self.gpu.then(|| GpuRastairModel {
            cpg: GpuForest::from_flat_forest(
                &FlatForest::from_forest(&combined.cpg, feature_nums.cpg),
                max_samples,
            ),
            denovo: GpuForest::from_flat_forest(
                &FlatForest::from_forest(&combined.denovo, feature_nums.denovo_cpg),
                max_samples,
            ),
            others: GpuForest::from_flat_forest(
                &FlatForest::from_forest(&combined.others, feature_nums.others),
                max_samples,
            ),
        });

        Ok(MachineLearning {
            threshold: self.ml,
            feature_calculator: combined.feature_set,
            model: Some(Box::new(combined)),
            gpu_prototype,
        })
    }

    pub fn threshold(&self) -> Option<Probability> {
        if self.no_ml { None } else { Some(self.ml) }
    }
}

#[instrument(level = "debug", skip_all)]
pub fn load_combined_model(path: Option<&Path>, built_in: &[u8]) -> Result<RastairModel> {
    fn load_model_from_file(path: &Path) -> Result<RastairModel> {
        ensure!(path.exists(), "Model file does not exist");
        let file = fs::read(path).wrap_err("Failed to read model file")?;
        load_combined(&file[..]).wrap_err("Failed to load combined model")
    }

    fn load_combined(reader: impl Read) -> Result<RastairModel> {
        let decompress = lz4::Decoder::new(reader).wrap_err("Failed to create LZ4 decoder")?;
        rmp_serde::from_read(decompress).wrap_err("Failed to deserialize combined model")
    }

    if let Some(path) = path {
        let model = load_model_from_file(path)
            .wrap_err_with(|| format!("Failed to load combined model from {path:?}"));
        debug!(?path, "Loaded combined model from file");
        model
    } else {
        let model = load_combined(built_in);
        debug!("Loaded built-in combined model");
        model
    }
}
