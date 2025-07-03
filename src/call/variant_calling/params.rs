use super::ErrorModel;
use crate::call::variant_calling::read_masking::ReadMaskParams;

#[derive(Debug, Clone, clap::Args)]
pub struct VariantCallingParams {
    /// The error model to use
    ///
    /// This should match the sequencing platform used to generate the data
    #[arg(long, default_value = "novaseq6000")]
    pub error_model: ErrorModel,

    /// Whether to keep overlapping reads
    #[arg(long, default_value_t = false)]
    pub keep_overlapping_reads: bool,

    /// Only look at sites that are CpG in the reference
    #[arg(long, default_value_t = false)]
    pub cpgs_only: bool,

    #[command(flatten)]
    read_masking: ReadMaskParams,
}

impl Default for VariantCallingParams {
    fn default() -> Self {
        Self {
            error_model: ErrorModel::Novaseq6000,
            keep_overlapping_reads: false,
            cpgs_only: false,
            read_masking: ReadMaskParams::default(),
        }
    }
}
