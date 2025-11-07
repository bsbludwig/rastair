use crate::utils::logging::ThisIsABug as _;
use color_eyre::{
    Result,
    eyre::{Context as _, ContextCompat},
};
use smol_str::SmolStr;
use std::fmt;

/// A genomic region with chromosome and coordinates
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Region {
    pub contig: SmolStr,
    /// 1-based start position (inclusive)
    pub start: u64,
    /// 1-based end position (inclusive)
    pub end: u64,
}

impl Region {
    pub fn range(&self) -> std::ops::Range<u64> {
        self.start..self.end
    }

    /// Returns true if the given position falls within this region's bounds
    pub fn contains(&self, pos: u64) -> bool {
        self.range().contains(&pos)
    }

    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }
}

#[test]
fn test_region_contains() {
    let region = Region { contig: "chr1".into(), start: 100, end: 200 };
    assert!(region.contains(150));
    assert!(!region.contains(50));
    assert!(!region.contains(250));
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Display for Region {
    /// Write 1-based region in the format `contig:start-end`
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}-{}", self.contig, self.start + 1, self.end)
    }
}

impl TryFrom<Region> for noodles::core::Region {
    type Error = color_eyre::Report;

    fn try_from(value: Region) -> Result<Self, Self::Error> {
        let start = usize::try_from(value.start).wrap_err("start position too large for usize")?;
        let start = noodles::core::Position::try_from(
            start
                .checked_add(1)
                .wrap_err("start position overflowed when converting to 1-based")?,
        )?;
        let end = usize::try_from(value.end).wrap_err("end position too large for usize")?;
        let end = noodles::core::Position::try_from(end)
            .wrap_err("end position too large for noodles")?;
        Ok(noodles::core::Region::new(value.contig.to_string(), start..=end))
    }
}

/// A complete genomic region that represents a full chromosome or a user-specified region
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum SelectedRegion {
    /// Represents the entire contig, e.g., a full chromosome
    EntireContig(Region),
    UserDefinedSubset {
        /// The region that is a subset of the entire contig
        region: Region,
        /// Last position in contig this region covers
        last_position: u64,
    },
}

impl std::ops::Deref for SelectedRegion {
    type Target = Region;

    fn deref(&self) -> &Self::Target {
        match self {
            SelectedRegion::EntireContig(region) => region,
            SelectedRegion::UserDefinedSubset { region, .. } => region,
        }
    }
}

/// A chunk of a larger genomic region used for processing data in segments
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ChunkRegion {
    pub region: Region,
    /// The last valid position in the full region this chunk belongs to
    pub last_position: u64,
}

impl std::ops::Deref for ChunkRegion {
    type Target = Region;

    fn deref(&self) -> &Self::Target {
        &self.region
    }
}

impl ChunkRegion {
    pub fn pos_to_idx(&self, pos: u32) -> Result<usize> {
        let segment_start_pos = usize::try_from(self.region.start)
            .wrap_err("Segment range does not fit in usize")
            .this_is_a_bug()?;
        usize::try_from(pos)
            .wrap_err("Position does not fit in usize")
            .this_is_a_bug()?
            .checked_sub(segment_start_pos)
            .wrap_err_with(|| format!("Position {pos} is not in segment {}", self.region))
    }
}
