#[derive(Debug, Clone, clap::Args)]
pub struct ThresholdConfig {
    /// The minimum variant allele frequency
    #[clap(long, default_value_t = 0.2)]
    pub m_vaf_min: f64,

    /// The minimum number of reads to call a position as methylated
    #[clap(long, default_value_t = 3)]
    pub m_min_depth: usize,
}

mod call;
pub use call::call;

mod filters;

#[cfg(test)]
mod tests;
