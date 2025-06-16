use super::ErrorModel;

#[derive(Debug, Clone, Default, clap::Args)]
pub struct VariantCallingParams {
    /// The error model to use
    ///
    /// This should match the sequencing platform used to generate the data
    #[arg(long, default_value = "novaseq6000")]
    pub error_model: ErrorModel,
}
