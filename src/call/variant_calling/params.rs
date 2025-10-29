use super::ErrorModel;
use crate::{
    call::variant_calling::{
        QualityFilterParams, read_flags::ReadFlags, read_masking::ReadMaskParams,
    },
    utils::cli,
};
use better_default::Default;

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct VariantCallingParams {
    /// The error model to use
    ///
    /// This should match the sequencing platform used to generate the data
    #[arg(long, default_value = "novaseq6000")]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub error_model: ErrorModel,

    /// Whether to keep overlapping reads
    #[arg(long, default_value_t = false)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub keep_overlapping_reads: bool,

    #[command(flatten)]
    pub quality: QualityFilterParams,

    #[command(flatten)]
    pub read_masking: ReadMaskParams,

    #[command(flatten)]
    pub read_flags: ReadFlags,
}
