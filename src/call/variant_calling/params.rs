use super::ErrorModel;
use crate::{
    call::variant_calling::{
        QualityFilterParams, read_flags::ReadFlags, read_masking::ReadMaskParams,
    },
    utils::cli,
};
use better_default::Default;

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
    #[arg(long, default_value_t = 30)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    #[default(30)]
    pub linear_dedup_threshold: usize,

    // The minimum number of reads to call a position as a variant
    #[arg(long, default_value_t = 3)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(3)]
    pub v_min_depth: u32,

    // The maximum number of reads to consider for a position (for performance reasons)
    #[arg(long, default_value_t = 1000)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(1000)]
    pub max_coverage: u32,

    #[command(flatten)]
    pub quality: QualityFilterParams,

    #[command(flatten)]
    pub read_masking: ReadMaskParams,

    #[command(flatten)]
    pub read_flags: ReadFlags,
}
