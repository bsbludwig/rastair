use crate::{call::process::IncludeAllCpGs, utils::cli, vcf};
use better_default::Default;
use color_eyre::{Result, eyre::Context};
use tracing::instrument;

#[derive(Debug, Clone, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct MethylationCallingParams {
    /// Calculate threshold values and filters for methylation
    #[arg(long, default_value_t = false)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    skip_methylation_calling: bool,

    #[command(flatten)]
    thresholds: ThresholdParams,
}

impl MethylationCallingParams {
    pub fn methylation_calling(&self) -> bool {
        !self.skip_methylation_calling
    }

    pub fn should_include_all_cpgs(&self) -> IncludeAllCpGs {
        if self.methylation_calling() { IncludeAllCpGs::Yes } else { IncludeAllCpGs::No }
    }

    /// Call methylation events based on the configured mode
    #[instrument(level = "trace", skip_all)]
    pub fn call(
        &self,
        record: &mut vcf::Record,
        before: Option<&vcf::Record>,
        after: Option<&vcf::Record>,
    ) -> Result<()> {
        if self.methylation_calling() {
            super::call(&self.thresholds, record, before, after)
                .wrap_err("Failed to call methylation based on thresholds")
        } else {
            // If no methylation calling is configured, just return the record as is
            Ok(())
        }
    }
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
