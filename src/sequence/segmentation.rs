use crate::{sequence::ChunkRegion, utils::cli};
use better_default::Default;
use color_eyre::eyre::ContextCompat as _;
use color_eyre::eyre::Result;
use seqair_types::{Base, Pos0, SmallVec};
use std::{
    num::NonZeroU64,
    ops::{Deref, Range},
};
use tracing::warn;

/// The range `--segment-max-length` accepts: a shorter segment only adds index queries and
/// shrinks ML batches, and a length of 0 never finishes segmenting.
pub fn segment_length() -> clap::builder::RangedU64ValueParser {
    clap::value_parser!(u64).range(1_000..=1_000_000)
}

#[derive(Debug, clap::Args, Clone, Default)]
pub struct SegmentationParams {
    /// Maximum length of a segment in bases
    ///
    /// Used for splitting work between threads. Tweak this to adjust memory
    /// usage: peak memory scales with this times `--threads`, and so does the
    /// batch each ML dispatch gets, which is why smaller is not free. On
    /// NA12878 chr12 at 26x, dropping this to 10 000 costs ~1.45x wall time —
    /// the same reads, but a tenth of the rows per GPU round trip and ten
    /// times as many BAM index queries. Between 1,000 and 1,000,000.
    #[arg(long, default_value_t = 100_000, value_parser = segment_length())]
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
    /// A segment whose index-estimated compressed size exceeds this is
    /// subdivided into smaller sub-segments, read one after the other. The
    /// decoded in-memory size is a few times larger than this compressed
    /// budget. `call` keeps every sub-segment's reads until the SNV engine has
    /// scored the segment, so there it bounds how much is read at a time, not
    /// peak memory.
    #[arg(long, default_value_t = ByteBudget::DEFAULT)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    #[default(ByteBudget::DEFAULT)]
    pub segment_max_bytes: ByteBudget,
}

/// The compressed bytes a segment may load, or no budget at all: `0` on the
/// command line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ByteBudget(Option<NonZeroU64>);

impl ByteBudget {
    pub const DEFAULT: Self = Self(NonZeroU64::new(256 * 1024 * 1024));

    /// A budget of `bytes`, or none for `0`, as the command line reads it.
    #[must_use]
    pub const fn new(bytes: u64) -> Self {
        Self(NonZeroU64::new(bytes))
    }

    #[must_use]
    pub const fn get(self) -> Option<NonZeroU64> {
        self.0
    }
}

impl std::fmt::Display for ByteBudget {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0.map_or(0, NonZeroU64::get))
    }
}

impl std::str::FromStr for ByteBudget {
    type Err = std::num::ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self::new(s.parse()?))
    }
}

impl SegmentationParams {
    pub fn sanitize(&mut self) {
        // `segment_length()` already holds `segment_max_length` to its range.
        if self.segment_overlap == 0 {
            let default = Self::default().segment_overlap;
            warn!(default, "Segment overlap is 0. This is invalid. Overwriting to default value.");
            self.segment_overlap = default;
        };
        let overlap_max_limit = 1000;
        if self.segment_overlap > overlap_max_limit {
            warn!(
                overlap = self.segment_overlap,
                "Segment overlap is set to a very high value (>{overlap_max_limit}). This may cause excessive redundant processing."
            );
        };
        if self.segment_overlap.saturating_mul(2) >= self.segment_max_length {
            warn!(
                overlap = self.segment_overlap,
                max_length = self.segment_max_length,
                "Segment overlap is more than half of segment max length. This may lead to inefficient processing."
            );
        };

        // Only a tiny budget is suspicious: it would subdivide regions far
        // more than necessary.
        let low_max_bytes = 1024 * 1024;
        if let Some(max_bytes) = self.segment_max_bytes.get()
            && max_bytes.get() < low_max_bytes
        {
            warn!(
                max_bytes,
                "Segment max bytes is set very low (<1 MiB); this may subdivide regions excessively."
            );
        };
    }
}

#[derive(Debug, Clone)]
pub struct Segment {
    pub range: std::sync::Arc<ChunkRegion>,
    pub sequence: Vec<u8>,
}

impl Deref for Segment {
    type Target = ChunkRegion;

    fn deref(&self) -> &Self::Target {
        &self.range
    }
}

impl Segment {
    /// Whether `pos` is one this segment emits, rather than one of the overlap
    /// positions a neighbouring segment owns.
    pub(crate) fn is_core(&self, pos: Pos0) -> bool {
        let core_start = self.region.start + self.overlap_start;
        let core_end = self.region.end.saturating_sub(self.overlap_end);
        (core_start..core_end).contains(&pos.as_u64())
    }

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sanitize_takes_any_overlap() {
        let mut params = SegmentationParams { segment_overlap: u64::MAX, ..Default::default() };
        params.sanitize();
        assert_eq!(params.segment_overlap, u64::MAX);
    }
}
