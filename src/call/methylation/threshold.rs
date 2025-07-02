#[derive(Debug, Clone, clap::Args)]
pub struct ThresholdConfig {
    /// The minimum variant allele frequency
    #[clap(long, default_value_t = 0.2)]
    pub m_vaf_min: f64,

    /// The minimum number of reads to call a position as methylated
    #[clap(long, default_value_t = 3)]
    pub m_min_depth: usize,

    /// The minimum quality ratio `(ad_alt*bq_alt + 1) / (ad_ref*bq_ref + 1)`
    #[clap(long, default_value_t = 0.27)]
    pub m_bq_ratio_min: f64,

    /// The minimum relative position in read for alt allele evidence
    #[clap(long, default_value_t = 0.2)]
    pub m_read_position_min: f64,

    /// The maximum relative position in read for alt allele evidence
    #[clap(long, default_value_t = 0.8)]
    pub m_read_position_max: f64,

    /// The maximum coverage depth for methylation calling
    #[clap(long, default_value_t = 1000)]
    pub m_max_coverage: usize,
}

mod call;
pub use call::call;
mod filters;
mod utils;

#[cfg(test)]
mod tests;
