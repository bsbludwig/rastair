use smol_str::SmolStr;

/// A genomic region with chromosome and coordinates
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Region {
    pub chromosome: SmolStr,
    /// 1-based start position (inclusive)
    pub start: u64,
    /// 1-based end position (inclusive)
    pub end: u64,
}

impl Region {
    /// Returns true if the given position falls within this region's bounds
    pub fn contains(&self, pos: u64) -> bool {
        (self.start..self.end).contains(&pos)
    }
}

/// A complete genomic region that represents a full chromosome or a user-specified region
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct FullRegion(pub Region);

impl std::ops::Deref for FullRegion {
    type Target = Region;

    fn deref(&self) -> &Self::Target {
        &self.0
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
