use crate::call::vcf;
use color_eyre::{Result, eyre::Context};
use std::fmt;
use tracing::instrument;

#[derive(Debug, Clone, clap::Args)]
pub struct MethylationCallingParams {
    /// The methylation calling mode
    #[arg(long, default_value_t = MethylationCallingMode::None)]
    pub calling: MethylationCallingMode,

    #[command(flatten)]
    thresholds: super::threshold::ThresholdConfig,
}

impl MethylationCallingParams {
    /// Call methylation events based on the configured mode
    #[instrument(level = "trace", skip_all)]
    pub fn call(&self, record: vcf::Record) -> Result<vcf::Record> {
        match self.calling {
            MethylationCallingMode::None => {
                // No methylation calling, return the record as is
                Ok(record)
            }
            MethylationCallingMode::Thresholds => super::threshold::call(record, &self.thresholds)
                .wrap_err("Failed to call methylation based on thresholds"),
            // MethylationCallingMode::ML => {
            //     // Placeholder for ML-based calling logic
            //     unimplemented!("ML-based methylation calling is not implemented yet")
            // }
        }
    }
}

#[derive(Debug, Clone, Copy, clap::ValueEnum)]
pub enum MethylationCallingMode {
    /// Don't perform methylation calling
    None,
    /// Call methylation events based on thresholds
    Thresholds,
    // /// Call methylation events based on ML model
    // ML,
}

impl fmt::Display for MethylationCallingMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MethylationCallingMode::None => write!(f, "none"),
            MethylationCallingMode::Thresholds => write!(f, "thresholds"),
            // MethylationCallingMode::ML => write!(f, "ml"),
        }
    }
}
