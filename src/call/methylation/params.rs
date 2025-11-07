use crate::utils::cli;
use better_default::Default;

#[derive(Debug, Clone, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct MethylationCallingParams {
    #[command(flatten)]
    pub thresholds: ThresholdParams,
}

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct ThresholdParams {
    /// The minimum variant allele frequency
    #[arg(long, default_value_t = 0.2)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(0.2)]
    pub m_vaf_min: f64,

    /// The minimum number of reads to call a position as methylated
    #[arg(long, default_value_t = 3)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(3)]
    pub m_min_depth: usize,

    /// The minimum number of reads required as evidence for a de novo CpG
    #[arg(long, default_value_t = 2)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(2)]
    pub m_min_denovo_depth: u32,

    /// The minimum quality ratio `(ad_alt*bq_alt + 1) / (ad_ref*bq_ref + 1)`
    #[arg(long, default_value_t = 0.27)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(0.27)]
    pub m_bq_ratio_min: f64,

    /// The minimum relative position in read for alt allele evidence
    #[arg(long, default_value_t = 0.2)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(0.2)]
    pub m_read_position_min: f64,

    /// The maximum relative position in read for alt allele evidence
    #[arg(long, default_value_t = 0.8)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(0.8)]
    pub m_read_position_max: f64,

    /// The maximum coverage depth for methylation calling
    #[arg(long, default_value_t = 1000)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(1000)]
    pub m_max_coverage: usize,
}
