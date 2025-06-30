use crate::{call::process::IncludeAllCpGs, vcf};
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
    // todo: make this a generic filter that can be called directly when looking at a pileup
    pub fn should_include_all_cpgs(&self) -> IncludeAllCpGs {
        if matches!(self.calling, MethylationCallingMode::None) {
            IncludeAllCpGs::No
        } else {
            IncludeAllCpGs::Yes
        }
    }

    /// Call methylation events based on the configured mode
    #[instrument(level = "trace", skip_all)]
    pub fn call(
        &self,
        record: &mut vcf::Record,
        before: Option<&vcf::Record>,
        after: Option<&vcf::Record>,
    ) -> Result<()> {
        match self.calling {
            MethylationCallingMode::None => {
                // No methylation calling, return the record as is
                Ok(())
            }
            MethylationCallingMode::Thresholds => {
                super::threshold::call(&self.thresholds, record, before, after)
                    .wrap_err("Failed to call methylation based on thresholds")
            } // MethylationCallingMode::ML => {
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
