use super::ErrorModel;
use crate::call::variant_calling::{
    QualityFilterParams, read_flags::ReadFlags, read_masking::ReadMaskParams,
};

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct VariantCallingParams {
    /// The error model to use
    ///
    /// This should match the sequencing platform used to generate the data
    #[arg(long, default_value = "novaseq6000")]
    pub error_model: ErrorModel,

    /// Whether to keep overlapping reads
    #[arg(long, default_value_t = false)]
    pub keep_overlapping_reads: bool,

    /// Only report positions that are CpGs in the reference or variants that would result in a new CpG
    #[arg(long, default_value_t = false)]
    pub cpgs_only: bool,

    #[command(flatten)]
    pub quality: QualityFilterParams,

    #[command(flatten)]
    pub read_masking: ReadMaskParams,

    #[command(flatten)]
    pub read_flags: ReadFlags,
}
