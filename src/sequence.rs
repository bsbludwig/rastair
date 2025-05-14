// given a region to scan
// - split into segments based on indexes
//   -  can i know the length of the region from fastq?
// - for each segment
//   - read the reference sequence from fastq
//   - for each position
//      - read the pileup
//      - include context of bases before/after
//      - calculate metrics
//
// segments are optional in the beginning

use crate::utils::{
    RegionString,
    file_helpers::{FastaReader, open_fasta},
};
use clap::value_parser;
use clio::ClioPath;
use color_eyre::{
    Result,
    eyre::{Context, ContextCompat, bail},
};
use regions::{ChunkRegion, FullRegion, Region};
use rust_htslib::bam::{self, FetchDefinition, Read as _};
use smol_str::SmolStr;
use tracing::{debug, instrument, trace};

mod chunked;
mod regions;

#[derive(Debug, clap::Args, Clone)]
pub struct SegmentsParams {
    /// A sorted and indexed bam file
    #[arg(value_parser=value_parser!(ClioPath).exists().is_file())]
    pub bam_file: ClioPath,

    /// A sorted and indexed (via samtools faidx) fasta file. Can be bgzip
    /// compressed, but requires both a gzi index and a fai index
    #[arg(short='r', long, value_parser=value_parser!(ClioPath).exists().is_file())]
    pub fasta_file: ClioPath,

    /// Restrict to a specific chromosome or region of a chromosome. Format is
    /// "chr", "chr:start" or "chr:start-end", where start is 1-based and end is
    /// inclusive.
    #[arg(short = 'l', long)]
    pub region: Option<RegionString>,

    /// Number of threads to use for reading the bam file
    #[arg(long, default_value_t = 4)]
    pub threads: usize,

    #[command(flatten)]
    pub segmentation: SegmentationParams,
}

#[derive(Debug, clap::Args, Clone)]
pub struct SegmentationParams {
    /// Maximum length of a segment in bases
    #[arg(long, default_value_t = 1_000_000)]
    pub max_segment_length: u64,

    /// Number of bases to overlap between segments
    #[arg(long, default_value_t = 100)]
    pub segment_overlap: u64,
}

impl SegmentsParams {
    pub fn readers(&self) -> Result<Readers> {
        let fasta = open_fasta(&self.fasta_file)?;

        let mut bam = bam::IndexedReader::from_path(self.bam_file.path())?;
        bam.set_threads(self.threads)?;

        Ok(Readers { fasta, bam, params: self.clone() })
    }
}

pub struct Readers {
    pub fasta: FastaReader,
    pub bam: bam::IndexedReader,
    params: SegmentsParams,
}

impl Readers {
    /// Calculate segments based on configuration parameters
    #[instrument(level = "debug", skip_all)]
    pub fn segments(&self) -> Result<impl Iterator<Item = ChunkRegion> + use<>> {
        let full_regions = if let Some(region) = &self.params.region {
            debug!(?region, "fetching specified region");
            let start = region.start.map(|x| x.get().into()).unwrap_or(1);
            let last_position = {
                // If no end specified, use chromosome length from BAM header
                let header = self.bam.header();
                header
                    .target_len(
                        header.tid(region.chromosome.as_bytes()).wrap_err("failed to fetch tid")?,
                    )
                    .wrap_err("failed to fetch header length")?
            };
            let end = region.end.map(|x| x.get().into()).unwrap_or(last_position);

            // Since the user specified a region, let's go with only that
            vec![FullRegion(Region { chromosome: region.chromosome.clone(), start, end })]
        } else {
            debug!("fetching all regions");
            get_full_regions(self.bam.header()).wrap_err("failed to get full region")?
        };

        if full_regions.is_empty() {
            bail!("No regions found");
        }

        let initial_start = full_regions[0].0.start;
        let chunked = chunked::ChunkedRegions {
            full_regions,
            current_region_idx: 0,
            current_start: initial_start,
            max_length: self.params.segmentation.max_segment_length,
            overlap: self.params.segmentation.segment_overlap,
        };

        Ok(chunked)
    }

    pub fn segment(&mut self, region: &ChunkRegion) -> Result<Segment> {
        let fetch_some_more = 2;
        let last_position_to_fetch =
            region.end.wrapping_add(fetch_some_more).min(region.last_position);

        // Calculate exact capacity needed to avoid reallocations
        let len = usize::try_from(last_position_to_fetch.wrapping_sub(region.start))
            .wrap_err("failed to convert segment length to usize")?;

        trace!(?region, len, "fetching segment");
        let mut seq = Vec::with_capacity(len);
        self.fasta.fetch(&region.chromosome, region.start, last_position_to_fetch)?;
        self.fasta.read(&mut seq)?;
        Ok(Segment { range: region.clone(), sequence: seq })
    }
}

#[derive(Debug)]
pub struct Segment {
    pub range: ChunkRegion,
    pub sequence: Vec<u8>,
}

impl<'seg> TryFrom<&'seg Segment> for FetchDefinition<'seg> {
    type Error = color_eyre::Report;

    fn try_from(segment: &'seg Segment) -> Result<Self> {
        Ok(FetchDefinition::RegionString(
            segment.range.region.chromosome.as_bytes(),
            i64::try_from(segment.range.region.start).wrap_err("start is invalid i64")?,
            i64::try_from(segment.range.region.end).wrap_err("end is invalid i64")?,
        ))
    }
}

#[instrument(level = "debug", skip(header))]
fn get_full_regions(header: &bam::HeaderView) -> Result<Vec<FullRegion>> {
    header
        .target_names()
        .iter()
        .enumerate()
        .filter(|(_, name)| !name.is_empty())
        .map(|(tid, name)| -> Result<FullRegion> {
            let chr = SmolStr::new(
                std::str::from_utf8(name).wrap_err("bam target name not valid UTF-8")?,
            );
            let length = header
                .target_len(u32::try_from(tid).wrap_err("faield to get tid")?)
                .wrap_err("failed to get target length")?;

            Ok(FullRegion(Region {
                chromosome: chr,
                start: 1, // 1-based coordinates
                end: length,
            }))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::path::PathBuf;

    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests").join("data")
    }

    fn get_test_bam() -> ClioPath {
        ClioPath::new(test_data_dir().join("test.bam")).expect("test bam path should be valid")
    }

    fn get_test_fasta() -> ClioPath {
        ClioPath::new(test_data_dir().join("test.fasta.gz"))
            .expect("test fasta path should be valid")
    }

    #[test]
    fn test_segment_reading() -> Result<()> {
        // Create a test region and params
        let region = ChunkRegion {
            region: Region { chromosome: "chr19".into(), start: 6105700, end: 6105800 },
            last_position: 6105900,
        };

        let params = SegmentsParams {
            bam_file: get_test_bam(),
            fasta_file: get_test_fasta(),
            region: None,
            threads: 4,
            segmentation: SegmentationParams { max_segment_length: 1000, segment_overlap: 100 },
        };

        // Initialize readers
        let mut readers = params.readers()?;

        // Test reading a single segment
        let segment = readers.segment(&region)?;

        // Verify segment properties
        assert_eq!(segment.range, region);

        Ok(())
    }

    #[test]
    fn test_calculate_segments() -> Result<()> {
        let params = SegmentsParams {
            bam_file: get_test_bam(),
            fasta_file: get_test_fasta(),
            region: Some("chr19:6105700-6105800".parse().unwrap()),
            threads: 4,
            segmentation: SegmentationParams { max_segment_length: 1000, segment_overlap: 100 },
        };

        let readers = params.readers()?;
        // Collect all segments into a Vec since we need to verify properties across all segments
        let segments = readers.segments()?.collect::<Vec<_>>();

        assert!(!segments.is_empty(), "Should have at least one segment");

        // Check segment properties
        for segment in segments {
            assert_eq!(segment.region.chromosome, "chr19");
            assert!(segment.region.start >= 6105700);
            assert!(segment.region.end <= 6105800);
        }

        Ok(())
    }

    #[test]
    fn test_overlapping_segments() -> Result<()> {
        let params = SegmentsParams {
            bam_file: get_test_bam(),
            fasta_file: get_test_fasta(),
            region: Some("chr19:6105700-6105900".parse().unwrap()),
            threads: 4,
            segmentation: SegmentationParams {
                max_segment_length: 100, // Small max length to force multiple segments
                segment_overlap: 20,     // Known overlap amount
            },
        };

        let mut readers = params.readers()?;
        let segments = readers.segments()?.collect::<Vec<_>>();

        assert!(segments.len() > 1, "Should have multiple segments");

        // Check overlaps between adjacent segments
        for pair in segments.windows(2) {
            let first = readers.segment(&pair[0])?;
            let second = readers.segment(&pair[1])?;

            // Verify overlap amount
            let overlap = first.range.region.end - second.range.region.start;
            assert_eq!(
                overlap, params.segmentation.segment_overlap,
                "Overlap amount should match configured value"
            );

            // Calculate overlap regions accounting for 0-based sequence indexing
            let first_start_idx =
                (first.range.region.end - overlap - first.range.region.start) as usize;
            let first_end_idx = (first.range.region.end - first.range.region.start) as usize;
            let second_start_idx = 0;
            let second_end_idx = overlap as usize;

            // Verify overlapping sequence content matches
            let overlap_first = &first.sequence[first_start_idx..first_end_idx];
            let overlap_second = &second.sequence[second_start_idx..second_end_idx];
            assert_eq!(overlap_first, overlap_second, "Overlapping sequences should match");
        }

        Ok(())
    }
}
