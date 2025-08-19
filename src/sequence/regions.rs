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
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}-{}", self.contig, self.start, self.end)
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
