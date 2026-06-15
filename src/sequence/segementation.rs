use crate::{sequence::ChunkRegion, utils::cli};
use better_default::Default;
use color_eyre::eyre::ContextCompat as _;
use color_eyre::eyre::Result;
use seqair_types::{Base, SmallVec};
use std::ops::{Deref, Range};
use tracing::warn;

#[derive(Debug, clap::Args, Clone, Default)]
pub struct SegmentationParams {
    /// Maximum length of a segment in bases
    ///
    /// Used for splitting work between threads. Tweak this to adjust memory
    /// usage.
    #[arg(long, default_value_t = 100_000)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    #[default(100_000)]
    pub segment_max_length: u64,

    /// Number of bases to overlap between segments
    ///
    /// Helpful to avoid missing variants at the edges of segments.
    #[arg(long, default_value_t = 200)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    #[default(200)]
    pub segment_overlap: u64,

    /// Maximum estimated compressed bytes to load per segment (memory budget)
    ///
    /// seqair backend only. A segment whose index-estimated compressed size
    /// exceeds this is subdivided into smaller sub-segments, so peak memory
    /// per worker stays bounded regardless of local coverage. The decoded
    /// in-memory size is a few times larger than this compressed budget.
    #[arg(long, default_value_t = 256 * 1024 * 1024)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    #[default(256 * 1024 * 1024)]
    pub segment_max_bytes: u64,
}

impl SegmentationParams {
    pub fn sanitize(&mut self) {
        let defaults = Self::default();

        if self.segment_max_length == 0 {
            let default = defaults.segment_max_length;
            warn!(
                default,
                "Segment max length is 0. This is invalid. Overwriting to default value."
            );
            self.segment_max_length = default;
        };

        let low_max_length = 1_000;
        if self.segment_max_length < low_max_length {
            warn!(
                max_length = self.segment_max_length,
                "Segment max length is set to a very low value (<{low_max_length}). This may cause performance degradation."
            );
        };

        let high_max_length = 1_000_000;
        if self.segment_max_length > high_max_length {
            warn!(
                max_length = self.segment_max_length,
                "Segment max length is set to a very high value (>{high_max_length}). This may cause high memory usage."
            );
        };

        if self.segment_overlap == 0 {
            let default = defaults.segment_overlap;
            warn!(
                default,
                "Segment max length is 0. This is invalid. Overwriting to default value."
            );
            self.segment_overlap = default;
        };
        let overlap_max_limit = 1000;
        if self.segment_overlap > overlap_max_limit {
            warn!(
                overlap = self.segment_overlap,
                "Segment overlap is set to a very high value (>{overlap_max_limit}). This may cause excessive redundant processing."
            );
        };
        if self.segment_overlap * 2 >= self.segment_max_length {
            warn!(
                overlap = self.segment_overlap,
                max_length = self.segment_max_length,
                "Segment overlap is more than half of segment max length. This may lead to inefficient processing."
            );
        };

        // `0` intentionally disables the budget; only warn for a tiny non-zero
        // value, which would subdivide regions far more than necessary.
        let low_max_bytes = 1024 * 1024;
        if self.segment_max_bytes != 0 && self.segment_max_bytes < low_max_bytes {
            warn!(
                max_bytes = self.segment_max_bytes,
                "Segment max bytes is set very low (<1 MiB); this may subdivide regions excessively."
            );
        };
    }
}

#[derive(Debug, Clone)]
pub struct Segment {
    pub range: ChunkRegion,
    pub sequence: Vec<u8>,
    /// Number of bases of overlap at the start of this segment
    pub overlap_start: u64,
    /// Number of bases of overlap at the end of this segment
    pub overlap_end: u64,
}

impl Deref for Segment {
    type Target = ChunkRegion;

    fn deref(&self) -> &Self::Target {
        &self.range
    }
}

impl Segment {
    /// Get a slice of the sequence
    pub fn sequence_slice<const N: usize>(
        &self,
        start: usize,
        end: usize,
    ) -> Result<SmallVec<Base, N>> {
        Ok(self.get(start, end)?.iter().map(Base::from).collect())
    }

    /// Get a slice of the sequence
    pub fn get(&self, start: usize, end: usize) -> Result<&[u8]> {
        let start = start.min(self.sequence.len());
        let end = end.min(self.sequence.len());

        self.sequence.get(start..end).wrap_err_with(|| FailedToReadSequenceSlice(start..end))
    }
}

#[derive(Debug, thiserror::Error)]
#[error("Failed to read sequence slice {0:?}")]
struct FailedToReadSequenceSlice(Range<usize>);
