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
    /// Accepts platform names or a custom error rate (e.g., 0.005)
    #[arg(long, default_value = "novaseq6000", value_parser = ErrorModel::value_parser())]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub error_model: ErrorModel,

    /// Whether to keep overlapping reads
    #[arg(long, default_value_t = false)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub keep_overlapping_reads: bool,

    // The minimum number of reads to call a position as a variant
    #[arg(long, default_value_t = 3)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(3)]
    pub v_min_depth: u32,

    // The minimum number of reads to call a position as a variant
    #[arg(long, default_value_t = 1000)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(1000)]
    pub max_coverage: u32,

    #[command(flatten)]
    pub quality: QualityFilterParams,

    #[command(flatten)]
    pub read_masking: ReadMaskParams,

    #[command(flatten)]
    pub read_flags: ReadFlags,
}
