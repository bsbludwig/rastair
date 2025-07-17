use crate::call::ml::{MachineLearning, models::load_model};
use clap::value_parser;
use clio::ClioPath;
use color_eyre::{Result, eyre::Context};

#[derive(Debug, Clone, clap::Args)]
pub struct MachineLearningParams {
    /// Use machine learning model with this threshold value to call variants
    /// and methylation events
    ///
    /// When specified, a ML model will classify positions with a prediction
    /// score. Anything above this threshold is considered PASS.
    #[arg(long = "ml", required = false, default_missing_value = "0.8", num_args = 0..=1)]
    pub ml: Option<f64>,
    /// Path to the model for CpG positions
    ///
    /// Default is the bundled model in the Rastair binary.
    #[arg(long, value_parser=value_parser!(ClioPath).exists().is_file())]
    model_cpg: Option<ClioPath>,
    /// Path to the model for de novo CpG positions
    ///
    /// Default is the bundled model in the Rastair binary.
    #[arg(long, value_parser=value_parser!(ClioPath).exists().is_file())]
    model_denovo_cpg: Option<ClioPath>,
    /// Path to the model for other positions
    ///
    /// Default is the bundled model in the Rastair binary.
    #[arg(long, value_parser=value_parser!(ClioPath).exists().is_file())]
    model_others: Option<ClioPath>,
}

impl MachineLearningParams {
    pub fn init(&self) -> Result<MachineLearning> {
        let Some(threshold) = self.ml else {
            return Ok(MachineLearning::disabled());
        };

        Ok(MachineLearning {
            threshold,
            cpg: Some(Box::new(
                load_model(
                    self.model_cpg.as_ref().map(|x| x.path()),
                    &include_bytes!("../../../models/BS_RF_800-2_CpG.rf.mpk.lz4")[..],
                )
                .wrap_err("Failed to load CpG RF model")?,
            )),
            denovo_cpg: Some(Box::new(
                load_model(
                    self.model_cpg.as_ref().map(|x| x.path()),
                    &include_bytes!("../../../models/BS_RF_800-2_denovo.rf.mpk.lz4")[..],
                )
                .wrap_err("Failed to load DeNovo CpG RF model")?,
            )),
            others: Some(Box::new(
                load_model(
                    self.model_cpg.as_ref().map(|x| x.path()),
                    &include_bytes!("../../../models/BS_RF_800-2_other.rf.mpk.lz4")[..],
                )
                .wrap_err("Failed to load Others RF model")?,
            )),
        })
    }
}
