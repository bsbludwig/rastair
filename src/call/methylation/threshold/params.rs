#[derive(Debug, Clone, clap::Args)]
pub struct ThresholdParams {
    /// The minimum variant allele frequency
    #[clap(long, default_value_t = 0.2)]
    pub m_vaf_min: f64,

    /// The minimum number of reads to call a position as methylated
    #[clap(long, default_value_t = 3)]
    pub m_min_depth: usize,

    /// The minimum number of reads required as evidence for a de novo CpG
    #[clap(long, default_value_t = 2)]
    pub m_min_denovo_depth: u32,

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

impl Default for ThresholdParams {
    fn default() -> Self {
        Self {
            m_vaf_min: 0.2,
            m_min_depth: 3,
            m_min_denovo_depth: 2,
            m_bq_ratio_min: 0.27,
            m_read_position_min: 0.2,
            m_read_position_max: 0.8,
            m_max_coverage: 1000,
        }
    }
}
