use super::ErrorModel;
use crate::{
    call::variant_calling::{
        QualityFilterParams, read_flags::ReadFlags, read_masking::ReadMaskParams,
    },
    utils::cli,
};
use better_default::Default;
use seqair::reader::DepthLimit;
use std::num::NonZeroU32;

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct VariantCallingParams {
    /// Enable unpaired mode
    ///
    /// In this mode, unpaired reads are accepted and strand assignment uses
    /// only the read's alignment direction (forward=OT, reverse=OB).
    #[arg(long, default_value_t = false)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub unpaired: bool,

    /// The error model to use
    ///
    /// Accepts platform names or a custom error rate (e.g., 0.005)
    #[arg(long, default_value = "novaseq6000", value_parser = ErrorModel::value_parser())]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub error_model: ErrorModel,

    /// Whether to keep overlapping paired-end reads
    ///
    /// In unpaired (`--single-strand`) mode this is ignored because read-pair
    /// overlap deduplication is disabled.
    #[arg(long, default_value_t = false)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub keep_overlapping_reads: bool,

    // The minimum number of reads to call a position as a variant
    #[arg(long, default_value_t = 3)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(3)]
    pub v_min_depth: u32,

    /// Maximum number of reads to consider at one position
    ///
    /// Counts reads that pass the flag, mapping-quality and tag filters, so a
    /// pileup of MAPQ 0 reads cannot crowd out the usable ones; where it binds,
    /// the reads starting first are kept. Base quality and read-end masking
    /// apply afterwards, so a capped position counts somewhat fewer. `0`
    /// means no limit.
    #[arg(long, default_value_t = MaxCoverage::DEFAULT)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(MaxCoverage::DEFAULT)]
    pub max_coverage: MaxCoverage,

    #[command(flatten)]
    pub quality: QualityFilterParams,

    #[command(flatten)]
    pub read_masking: ReadMaskParams,

    #[command(flatten)]
    pub read_flags: ReadFlags,
}

/// How many reads to consider at one position, or no limit at all.
///
/// Not a `u32`, because `0` on the command line means *no limit*, and a plain
/// number invites `depth.min(0)`: zero reads at every position. This type has
/// no number to misread, only [`depth_limit`](MaxCoverage::depth_limit).
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(transparent)]
pub struct MaxCoverage(Option<NonZeroU32>);

impl MaxCoverage {
    /// Deep enough that no real pileup reaches it, shallow enough to bound a
    /// collapsed repeat.
    pub const DEFAULT: Self = Self(NonZeroU32::new(1000));

    /// A cap of `n` reads, or no cap for `0` — the same reading the command
    /// line gets.
    #[must_use]
    pub const fn new(n: u32) -> Self {
        Self(NonZeroU32::new(n))
    }

    /// The cap seqair applies while loading, after the record filter.
    #[must_use]
    pub fn depth_limit(self) -> DepthLimit {
        self.0.map_or(DepthLimit::Unlimited, DepthLimit::PerColumn)
    }
}

impl Default for MaxCoverage {
    fn default() -> Self {
        Self::DEFAULT
    }
}

impl std::fmt::Display for MaxCoverage {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0.map_or(0, NonZeroU32::get))
    }
}

impl std::str::FromStr for MaxCoverage {
    type Err = std::num::ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self(NonZeroU32::new(s.parse()?)))
    }
}
