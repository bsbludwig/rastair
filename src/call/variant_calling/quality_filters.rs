use crate::call::variants::SeenBase;

#[derive(Debug, Clone, clap::Args)]
pub struct QualityFilterParams {
    /// Minimum mapping quality to consider a read
    #[arg(short = 'q', long, default_value_t = 1)]
    pub min_mapq: u8,
    /// Minimum base quality to consider a base
    #[arg(short = 'Q', long, default_value_t = 10)]
    pub min_baseq: u8,
}

impl Default for QualityFilterParams {
    fn default() -> Self {
        // keep in sync with the default values in the CLI
        Self { min_mapq: 1, min_baseq: 10 }
    }
}

impl QualityFilterParams {
    /// Check if the read and base pass the quality filters
    pub fn filter(&self, read: &SeenBase) -> bool {
        read.mapq >= self.min_mapq && read.qual >= self.min_baseq
    }
}
