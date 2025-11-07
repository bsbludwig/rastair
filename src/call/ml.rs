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

use crate::{metrics::ml::types::MachineLearning, utils::cli};
use biosphere::RandomForest;
use clap::value_parser;
use clio::ClioPath;
use color_eyre::{
    Result,
    eyre::{Context, ensure},
};
use rastair_types::Probability;
use std::{fs, io::Read, path::Path};
use tracing::instrument;

pub const DEFAULT_ML_THRESHOLD: Probability = Probability::new_panicky(0.8);

#[derive(Debug, Clone, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct MachineLearningParams {
    /// Only use hard thresholds to call variants and methylation events.
    ///
    /// This disables using the machine learning models. This will make rastair
    /// much faster, but at the cost of accuracy.
    #[arg(long = "thresholds")]
    #[arg(help_heading = cli::sections::FILTER)]
    pub no_ml: bool,
    /// Use machine learning model with this threshold value to call variants
    /// and methylation events
    ///
    /// When specified, a ML model will classify positions with a prediction
    /// score. Anything above this threshold is considered PASS.
    ///
    /// For consistency with `--thresholds`, this option can be also be
    /// specified as `--ml` without a value, which will use the default
    /// threshold.
    #[arg(long = "ml", default_value_t = DEFAULT_ML_THRESHOLD, default_missing_value = "0.8", num_args = 0..=1)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub ml: Probability,
    /// Path to the model for CpG positions
    ///
    /// Default is the bundled model in the Rastair binary.
    #[arg(long, value_parser=value_parser!(ClioPath).exists().is_file())]
    #[arg(help_heading = cli::sections::FILTER)]
    #[serde(skip)]
    model_cpg: Option<ClioPath>,
    /// Path to the model for de novo CpG positions
    ///
    /// Default is the bundled model in the Rastair binary.
    #[arg(long, value_parser=value_parser!(ClioPath).exists().is_file())]
    #[arg(help_heading = cli::sections::FILTER)]
    #[serde(skip)]
    model_denovo_cpg: Option<ClioPath>,
    /// Path to the model for other positions
    ///
    /// Default is the bundled model in the Rastair binary.
    #[arg(long, value_parser=value_parser!(ClioPath).exists().is_file())]
    #[arg(help_heading = cli::sections::FILTER)]
    #[serde(skip)]
    model_others: Option<ClioPath>,
}

impl MachineLearningParams {
    #[instrument(name = "init_ml", skip(self))]
    pub fn init(&self) -> Result<MachineLearning> {
        if self.no_ml {
            return Ok(MachineLearning::disabled());
        };

        Ok(MachineLearning {
            disabled: false,
            threshold: self.ml,
            cpg: Some(Box::new(
                load_model(
                    self.model_cpg.as_ref().map(|x| x.path()),
                    &include_bytes!("../../models/BS_RF_800-2_CpG.rf.mpk.lz4")[..],
                )
                .wrap_err("Failed to load CpG RF model")?,
            )),
            denovo_cpg: Some(Box::new(
                load_model(
                    self.model_cpg.as_ref().map(|x| x.path()),
                    &include_bytes!("../../models/BS_RF_800-2_denovo.rf.mpk.lz4")[..],
                )
                .wrap_err("Failed to load DeNovo CpG RF model")?,
            )),
            others: Some(Box::new(
                load_model(
                    self.model_cpg.as_ref().map(|x| x.path()),
                    &include_bytes!("../../models/BS_RF_800-2_other.rf.mpk.lz4")[..],
                )
                .wrap_err("Failed to load Others RF model")?,
            )),
        })
    }
}

#[instrument(level = "debug", skip_all)]
pub fn load_model(path: Option<&Path>, built_in: &[u8]) -> Result<RandomForest> {
    fn load_model_from_file(path: &Path) -> Result<RandomForest> {
        ensure!(path.exists(), "Model file does not exist");
        let file = fs::read(path).wrap_err("Failed to read model file")?;
        load_rf(&file[..]).wrap_err("Failed to load model")
    }

    fn load_rf(reader: impl Read) -> Result<RandomForest> {
        let decompress = lz4::Decoder::new(reader).wrap_err("Failed to create LZ4 decoder")?;
        rmp_serde::from_read(decompress).wrap_err("Failed to deserialize random forest")
    }

    let Some(path) = path else {
        return load_rf(built_in);
    };
    load_model_from_file(path).wrap_err_with(|| format!("Failed to load model from {path:?}"))
}
