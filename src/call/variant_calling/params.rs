use super::ErrorModel;
use crate::{
    call::variant_calling::{
        QualityFilterParams, read_flags::ReadFlags, read_masking::ReadMaskParams,
    },
    utils::cli,
};
use better_default::Default;
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

    /// Guess OT/OB read orientation from mismatch motifs instead of SAM flags
    ///
    /// Scans read mismatches against the reference and counts `TG` versus `CA`
    /// motifs in a 2 bp window anchored at each mismatch (current+next and
    /// previous+current), using the htslib/reference-oriented read sequence.
    /// `TG > CA` means OT, `CA > TG` means OB, and ties / evidence-free reads
    /// are split pseudo-randomly but reproducibly per read.
    ///
    /// Useful for non-directional libraries, commonly from tagmentation prep.
    ///
    /// This option currently affects only the pileup-based `call` workflow.
    #[arg(long, default_value_t = false)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub guess_read_orientation: bool,

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

    /// Depth threshold below which linear name dedup is used instead of a hashmap
    ///
    /// At pileup positions with depth ≤ this value, read name deduplication
    /// uses a linear scan through parallel suffix/name arrays rather than an
    /// `FxHashMap`. Set to 0 to always use the hashmap.
    ///
    /// Only affects the htslib pileup path. The seqair path deduplicates
    /// overlapping mates from seqair's per-store mate links and never matches
    /// read names, so it ignores this.
    #[arg(long, default_value_t = 30)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    #[default(30)]
    pub linear_dedup_threshold: usize,

    // The minimum number of reads to call a position as a variant
    #[arg(long, default_value_t = 3)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(3)]
    pub v_min_depth: u32,

    /// Maximum number of reads to consider at one position
    ///
    /// Counts reads that pass the mapping- and base-quality filters, so a
    /// pileup of MAPQ 0 reads cannot crowd out the usable ones. To bound
    /// memory, at most ten times this many reads are loaded per position at
    /// all. `0` means no limit.
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
/// Not a `u32`, because `0` on the command line means *no limit* and a plain
/// number cannot say so: `depth.min(0)` is zero reads at every position, which
/// is what rastair did, while the depth cap handed to seqair read the same `0`
/// as "unlimited" two lines away and loaded every read. The run produced an
/// empty VCF at unbounded memory. Through this type neither reading is
/// available — there is no number to `min` against, only [`clamp`] and
/// [`per_column`], and "unlimited" is [`None`] in both.
///
/// [`clamp`]: MaxCoverage::clamp
/// [`per_column`]: MaxCoverage::per_column
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(transparent)]
pub struct MaxCoverage(Option<NonZeroU32>);

impl MaxCoverage {
    /// Deep enough that no real pileup reaches it, shallow enough to bound a
    /// collapsed repeat.
    pub const DEFAULT: Self = Self(NonZeroU32::new(1000));

    /// No cap at all.
    pub const UNLIMITED: Self = Self(None);

    /// A cap of `n` reads, or [`UNLIMITED`](Self::UNLIMITED) for `0` — the
    /// same reading the command line gets.
    #[must_use]
    pub const fn new(n: u32) -> Self {
        Self(NonZeroU32::new(n))
    }

    /// The cap as seqair's [`DepthLimit`] wants it: `None` is unlimited.
    ///
    /// [`DepthLimit`]: seqair::reader::DepthLimit
    #[must_use]
    pub fn per_column(self) -> Option<NonZeroU32> {
        self.0
    }

    /// `depth`, capped. Unlimited returns it unchanged — which is the whole
    /// point, since `depth.min(cap)` cannot express that.
    #[must_use]
    pub fn clamp(self, depth: usize) -> usize {
        match self.0 {
            Some(cap) => depth.min(cap.get() as usize),
            None => depth,
        }
    }

    /// How many reads a column may *count*: those that pass the read filters,
    /// up to the cap. Unlimited is `usize::MAX`.
    #[must_use]
    pub fn kept(self) -> usize {
        self.0.map_or(usize::MAX, |cap| cap.get() as usize)
    }

    /// How many reads a column may *load*, before any filter: ten times the
    /// cap. The cap itself is applied to filtered reads, so it must not be
    /// the engine's truncation limit — a pileup that is 85 % MAPQ 0 would
    /// otherwise have its usable reads evicted by the useless ones, and
    /// which ones survive would depend on the engine's truncation order.
    /// The ceiling only bounds memory on collapsed repeats; the two engines
    /// still truncate differently above it, which is why it is generous.
    #[must_use]
    pub fn load_ceiling(self) -> Option<NonZeroU32> {
        const LOAD_FACTOR: NonZeroU32 = NonZeroU32::new(10).unwrap();
        self.0.map(|cap| cap.saturating_mul(LOAD_FACTOR))
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
