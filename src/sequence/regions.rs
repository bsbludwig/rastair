use crate::utils::logging::ThisIsABug as _;
use color_eyre::{
    Result,
    eyre::{Context as _, ContextCompat, ensure},
};
use rust_htslib::bam::FetchDefinition;
use seqair_types::SmolStr;
use std::fmt;
use tracing::warn;

/// A genomic region with chromosome and coordinates
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, serde::Serialize, serde::Deserialize)]
pub struct Region {
    pub contig: SmolStr,
    /// 0-based start position (inclusive)
    pub start: u64,
    /// 0-based end position (inclusive)
    pub end: u64,
}

impl Region {
    pub fn range(&self) -> std::ops::RangeInclusive<u64> {
        self.start..=self.end
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

impl<'reg> TryFrom<&'reg Region> for FetchDefinition<'reg> {
    type Error = color_eyre::Report;

    fn try_from(region: &'reg Region) -> Result<Self> {
        let start = i64::try_from(region.start).wrap_err("start is invalid i64")?;
        let end = i64::try_from(region.end.saturating_add(1)).wrap_err("end is invalid i64")?;

        Ok(FetchDefinition::RegionString(region.contig.as_bytes(), start, end))
    }
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

/// Restrict a set of regions to those whose contig exists in the FASTA reference.
///
/// `explicit` marks regions the user named on the CLI: those must all be present,
/// so a missing one is a hard error. Auto-derived regions (the whole BAM header)
/// are instead dropped with a warning — a BAM routinely carries decoy/alt contigs
/// absent from a slimmer FASTA. Bailing when nothing remains surfaces a wholesale
/// reference/naming mismatch (e.g. `chr1` vs `1`) instead of silently writing an
/// empty VCF.
pub(crate) fn retain_fasta_regions(
    regions: Vec<SelectedRegion>,
    fasta_has_contig: impl Fn(&str) -> bool,
    explicit: bool,
) -> Result<Vec<SelectedRegion>> {
    if explicit {
        for region in &regions {
            ensure!(
                fasta_has_contig(&region.contig),
                "Requested contig {:?} is not present in the FASTA reference",
                region.contig.as_str(),
            );
        }
        return Ok(regions);
    }

    let mut skipped: Vec<SmolStr> = Vec::new();
    let kept: Vec<SelectedRegion> = regions
        .into_iter()
        .filter(|region| {
            let present = fasta_has_contig(&region.contig);
            if !present {
                skipped.push(region.contig.clone());
            }
            present
        })
        .collect();

    if !skipped.is_empty() {
        warn!(
            skipped = skipped.iter().map(SmolStr::as_str).collect::<Vec<_>>().join(", "),
            "Skipping {} contig(s) from the BAM header that are absent from the FASTA reference",
            skipped.len(),
        );
    }

    ensure!(
        !kept.is_empty(),
        "None of the BAM contigs are present in the FASTA reference; check that the \
         BAM and FASTA share the same reference and contig naming",
    );

    Ok(kept)
}

/// A chunk of a larger genomic region used for processing data in segments
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, serde::Serialize, serde::Deserialize)]
pub struct ChunkRegion {
    pub region: Region,
    /// The last valid position in the full region this chunk belongs to
    pub last_position: u64,
    /// Number of bases of overlap at the start of this chunk
    pub overlap_start: u64,
    /// Number of bases of overlap at the end of this chunk
    pub overlap_end: u64,
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn test_region_contains() {
        let region = Region { contig: "chr1".into(), start: 100, end: 200 };
        assert!(region.contains(100));
        assert!(region.contains(150));
        assert!(region.contains(200));
        assert!(!region.contains(50));
        assert!(!region.contains(250));
    }

    fn entire(contig: &str) -> SelectedRegion {
        SelectedRegion::EntireContig(Region { contig: contig.into(), start: 1, end: 100 })
    }

    fn fasta_with(contigs: &[&str]) -> impl Fn(&str) -> bool {
        let set: HashSet<String> = contigs.iter().map(|c| c.to_string()).collect();
        move |contig: &str| set.contains(contig)
    }

    #[test]
    fn auto_derived_drops_contigs_missing_from_fasta() -> Result<()> {
        let regions = vec![entire("chr1"), entire("decoy"), entire("chr2")];
        let kept = retain_fasta_regions(regions, fasta_with(&["chr1", "chr2"]), false)?;
        let names: Vec<&str> = kept.iter().map(|r| r.contig.as_str()).collect();
        assert_eq!(names, ["chr1", "chr2"]);
        Ok(())
    }

    #[test]
    fn auto_derived_keeps_all_when_fasta_matches() -> Result<()> {
        let regions = vec![entire("chr1"), entire("chr2")];
        let kept = retain_fasta_regions(regions, fasta_with(&["chr1", "chr2", "chr3"]), false)?;
        assert_eq!(kept.len(), 2);
        Ok(())
    }

    #[test]
    fn auto_derived_bails_when_nothing_matches() {
        let regions = vec![entire("1"), entire("2")];
        // Whole-reference naming mismatch (`chr1` vs `1`): must error, not silently
        // produce an empty result.
        let err = retain_fasta_regions(regions, fasta_with(&["chr1", "chr2"]), false).unwrap_err();
        assert!(err.to_string().contains("None of the BAM contigs"));
    }

    #[test]
    fn explicit_region_missing_from_fasta_is_hard_error() {
        let regions = vec![entire("chr1"), entire("decoy")];
        let err = retain_fasta_regions(regions, fasta_with(&["chr1"]), true).unwrap_err();
        assert!(err.to_string().contains("decoy"));
    }

    #[test]
    fn explicit_regions_all_present_are_kept() -> Result<()> {
        let regions = vec![entire("chr1"), entire("chr2")];
        let kept = retain_fasta_regions(regions, fasta_with(&["chr1", "chr2"]), true)?;
        assert_eq!(kept.len(), 2);
        Ok(())
    }
}
