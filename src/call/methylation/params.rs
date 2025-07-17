use crate::{call::process::IncludeAllCpGs, vcf};
use color_eyre::{Result, eyre::Context};
use tracing::instrument;

#[derive(Debug, Clone, clap::Args)]
pub struct MethylationCallingParams {
    /// Calculate threshold values and filters for methylation
    #[arg(long, default_value_t = false)]
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
