use crate::{call::ml::DEFAULT_ML_THRESHOLD, utils::cli};
use clap::value_parser;
use clio::ClioPath;
use color_eyre::eyre::Result;
use rastair_types::Probability;
use tracing::instrument;

#[derive(Debug, clap::Args)]
pub struct TrainModelParams {
    /// Path to sorted and indexed BAM file
    #[arg(help_heading = cli::sections::INPUT, value_hint=clap::ValueHint::FilePath)]
    pub bam: ClioPath,

    /// Path to sorted and indexed (via samtools faidx) FASTA file. Can be bgzip
    /// compressed, but requires both a gzi index and a fai index
    #[arg(short='r', long, value_parser=value_parser!(ClioPath).exists().is_file(), value_hint=clap::ValueHint::FilePath)]
    #[arg(help_heading = cli::sections::INPUT)]
    pub fasta_file: ClioPath,

    /// Path to the ground truth file (VCF) to train with
    #[arg(help_heading = cli::sections::INPUT, value_hint=clap::ValueHint::FilePath)]
    pub truth: ClioPath,

    /// Use machine learning model with this threshold value to call variants
    /// and methylation events
    #[arg(long = "ml", default_value_t = DEFAULT_ML_THRESHOLD, default_missing_value = "0.8", num_args = 0..=1)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub ml: Probability,
}

#[instrument(level = "info", skip_all)]
pub fn train_model(params: &TrainModelParams) -> Result<()> {
    unimplemented!();
}
