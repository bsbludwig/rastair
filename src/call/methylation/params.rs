use crate::utils::cli;
use better_default::Default;

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct MethylationCallingParams {
    #[command(flatten)]
    pub thresholds: ThresholdParams,

    /// Rescue soft-clipped CpG-partner bases
    ///
    /// Aligners often soft-clip a read's fringe base(s) when that mismatches
    /// the reference. But in TAPS, a methylated C reads as T, so a leading C→T
    /// mismatch might get soft-clipped, discarding real methylation evidence
    /// (and same for a trailing G->A). With this flag, a single soft-clipped
    /// base immediately next to an aligned base is "rescued" when it is the
    /// missing partner of a reference CpG (ref C on OT, ref G on OB) and
    /// counted as a normal observation.
    ///
    /// Recovered bases bypass read-end masking (`--nOT`/`--nOB`) since they are
    /// fringe bases by definition, but still go through base-quality,
    /// mapping-quality and methylation read-position filters.
    ///
    /// (Only takes effect on the seqair backend.)
    #[arg(long)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub rescue_soft_clip_cpg: bool,
}

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct ThresholdParams {
    /// The minimum variant allele frequency
    #[arg(long, default_value_t = 0.2)]
    #[arg(help_heading = cli::sections::METHYLATION)]
    #[default(0.2)]
    pub m_vaf_min: f64,

    /// The minimum number of reads to call a position as methylated
    #[arg(long, default_value_t = 3)]
    #[arg(help_heading = cli::sections::METHYLATION)]
    #[default(3)]
    pub m_min_depth: u32,

    /// The minimum quality ratio `(ad_alt*bq_alt + 1) / (ad_ref*bq_ref + 1)`
    #[arg(long, default_value_t = 0.27)]
    #[arg(help_heading = cli::sections::METHYLATION)]
    #[default(0.27)]
    pub m_bq_ratio_min: f64,

    /// The minimum relative position in read for alt allele evidence
    #[arg(long, default_value_t = 0.2)]
    #[arg(help_heading = cli::sections::METHYLATION)]
    #[default(0.2)]
    pub m_read_position_min: f64,

    /// The maximum relative position in read for alt allele evidence
    #[arg(long, default_value_t = 0.8)]
    #[arg(help_heading = cli::sections::METHYLATION)]
    #[default(0.8)]
    pub m_read_position_max: f64,

    /// The maximum coverage depth for methylation calling
    #[arg(long, default_value_t = 1000)]
    #[arg(help_heading = cli::sections::METHYLATION)]
    #[default(1000)]
    pub m_max_coverage: u32,
}
