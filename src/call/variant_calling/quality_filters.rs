use crate::{call::pileup::SimpleRead, utils::cli};
use better_default::Default;

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct QualityFilterParams {
    /// Minimum mapping quality to consider a read
    #[arg(short = 'q', long, default_value_t = 1)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(1)]
    pub min_mapq: u8,
    /// Minimum base quality to consider a base
    #[arg(short = 'Q', long, default_value_t = 10)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(10)]
    pub min_baseq: u8,
}

impl QualityFilterParams {
    /// Check if the read and base pass the quality filters
    pub fn filter(&self, read: &SimpleRead) -> bool {
        read.mapq >= self.min_mapq && read.qual >= self.min_baseq
    }
}
