//! Segment processing for genomic analysis.
//!
//! This module provides functionality for:
//! - Dividing genomic regions into manageable segments
//! - Reading reference sequences from FASTA files
//! - Accessing aligned reads from BAM files
//! - Processing segments with configurable overlap between chunks

use crate::utils::{
    Base, RegionString,
    file_helpers::{FastaReader, open_fasta},
};
use clap::value_parser;
use clio::ClioPath;
use color_eyre::{
    Result, Section,
    eyre::{Context, ContextCompat, ensure},
};
pub use regions::{ChunkRegion, FullRegion, Region};
use rust_htslib::bam::{self, FetchDefinition, HeaderView, Read as _};
use smallvec::SmallVec;
use smol_str::SmolStr;
use std::num::NonZeroU32;
use tracing::{debug, instrument, trace};

mod chunked;
mod regions;

#[derive(Debug, clap::Args, Clone)]
pub struct ReaderParams {
    /// Path to sorted and indexed BAM file
    #[arg(value_parser=value_parser!(ClioPath).exists().is_file())]
    pub bam_file: ClioPath,

    /// Path to sorted and indexed (via samtools faidx) FASTA file. Can be bgzip
    /// compressed, but requires both a gzi index and a fai index
    #[arg(short='r', long, value_parser=value_parser!(ClioPath).exists().is_file())]
    pub fasta_file: ClioPath,

    /// Restrict to a specific chromosome or region of a chromosome. Format is
    /// "chr", "chr:start" or "chr:start-end", where start is 1-based and end is
    /// inclusive.
    #[arg(short = 'l', long)]
    pub region: Option<RegionString>,
}

impl ReaderParams {
    pub fn readers(&self) -> Result<Readers> {
        let fasta = open_fasta(&self.fasta_file)?;
        let bam_path = self.bam_file.path();
        let bam = bam::IndexedReader::from_path(bam_path)
            .with_suggestion(|| {
                format!(
                    "Ensure the BAM file is sorted and indexed with \
                    `samtools sort {bam_path:?}` and `samtools index {bam_path:?}`, respectively."
                )
            })
            .note("If you have a .bai file, ensure it is in the same directory as the BAM file.")?;

        Ok(Readers { fasta, bam, params: self.clone() })
    }
}

pub struct Readers {
    pub fasta: FastaReader,
    pub bam: bam::IndexedReader,
    params: ReaderParams,
}

impl Readers {
    /// Calculate segments based on configuration parameters
    #[instrument(level = "debug", skip_all)]
    pub fn segments(
        &self,
        segment_max_length: u64,
        segment_overlap: u64,
    ) -> Result<impl Iterator<Item = ChunkRegion> + use<>> {
        let full_regions = if let Some(region) = &self.params.region {
            vec![
                get_selected_region(region, self.bam.header())
                    .wrap_err("Failed to get selected region from BAM file")?,
            ]
        } else {
            debug!("fetching all regions");
            get_full_regions(self.bam.header())
                .wrap_err("Failed to get all regions from BAM file")?
        };
        ensure!(!full_regions.is_empty(), "No regions found");

        let initial_start = full_regions[0].0.start;
        let chunked = chunked::ChunkedRegions {
            full_regions,
            current_region_idx: 0,
            current_start: initial_start,
            max_length: segment_max_length,
            overlap: segment_overlap,
        };

        Ok(chunked)
    }

    #[instrument(level = "debug", skip_all)]
    pub fn segment(&mut self, region: &ChunkRegion, overfetch: u64) -> Result<Segment> {
        let last_position_to_fetch = region.end.wrapping_add(overfetch).min(region.last_position);

        // Calculate exact capacity needed to avoid reallocations
        let len = usize::try_from(last_position_to_fetch.wrapping_sub(region.start))
            .wrap_err("Failed to convert segment length to usize")?;

        trace!(?region, len, "fetching segment");
        let mut seq = Vec::with_capacity(len);
        self.fasta
            .fetch(&region.contig, region.start, last_position_to_fetch)
            .wrap_err("Failed to fetch region")
            .and_then(|_| self.fasta.read(&mut seq).wrap_err("Failed to read sequence from region"))
            // chain the calls so we can add this nice error:
            .wrap_err_with(|| format!("Failed to get region {} from FASTA file", region.region))?;

        Ok(Segment { range: region.clone(), sequence: seq })
    }
}

#[derive(Debug)]
pub struct Segment {
    pub range: ChunkRegion,
    pub sequence: Vec<u8>,
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

        self.sequence
            .get(start..end)
            .wrap_err_with(|| format!("Failed to read sequence slice {:?}", start..end))
    }
}

impl<'seg> TryFrom<&'seg Segment> for FetchDefinition<'seg> {
    type Error = color_eyre::Report;

    fn try_from(segment: &'seg Segment) -> Result<Self> {
        Ok(FetchDefinition::RegionString(
            segment.range.region.contig.as_bytes(),
            i64::try_from(segment.range.region.start).wrap_err("start is invalid i64")?,
            i64::try_from(segment.range.region.end).wrap_err("end is invalid i64")?,
        ))
    }
}

#[instrument(level = "debug", skip(bam_header))]
fn get_selected_region(region: &RegionString, bam_header: &HeaderView) -> Result<FullRegion> {
    let chromosome = region.chromosome.as_str();
    // If no start position is specified, default to beginning of chromosome
    let start = region.start.map(to_u64).unwrap_or(1);

    let target_id = bam_header
        .tid(region.chromosome.as_bytes())
        .wrap_err_with(|| {
            format!("Failed to fetch target ID for chromosome {} from header", region.chromosome)
        })
        .with_note(|| {
            format!(
                "This usually means the specified chromosome {} is not in the input BAM file",
                region.chromosome
            )
        })?;
    let last_position =
        bam_header.target_len(target_id).wrap_err("Failed to fetch header length")?;
    // If no end specified, use chromosome length from BAM header
    let end = region.end.map(to_u64).unwrap_or(last_position);

    ensure!(
        start <= last_position,
        "Specified start position {end} past the end of chromosome {chromosome}"
    );
    ensure!(
        end <= last_position,
        "Specified end position {end} past the end of chromosome {chromosome}"
    );

    // Since the user specified this region, we're only returning that one
    Ok(FullRegion(Region { contig: region.chromosome.clone(), start, end }))
}

fn to_u64(value: NonZeroU32) -> u64 {
    value.get().into()
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
                std::str::from_utf8(name).wrap_err("BAM target name not valid UTF-8")
                    .note("This is against the BAM specification, please check with the tool that created this file")?,
            );
            let length = header
                .target_len(u32::try_from(tid).wrap_err("Failed to get a target ID that was part of the BAM header")
                    .note("The BAM header might be corrupt")?)
                .wrap_err("Failed to get target length")?;

            Ok(FullRegion(Region {
                contig: chr,
                start: 1, // 1-based coordinates
                end: length,
            }))
        })
        .collect()
}

#[cfg(test)]
#[allow(clippy::cast_possible_truncation)]
mod tests {
    use proptest::{prop_assume, proptest};

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
            region: Region { contig: "chr19".into(), start: 6105700, end: 6105800 },
            last_position: 6105900,
        };

        let params =
            ReaderParams { bam_file: get_test_bam(), fasta_file: get_test_fasta(), region: None };

        // Initialize readers
        let mut readers = params.readers()?;

        // Test reading a single segment
        let segment = readers.segment(&region, 2)?;

        // Verify segment properties
        assert_eq!(segment.range, region);

        Ok(())
    }

    #[test]
    fn test_calculate_segments() -> Result<()> {
        let params = ReaderParams {
            bam_file: get_test_bam(),
            fasta_file: get_test_fasta(),
            region: Some("chr19:6105700-6105800".parse().unwrap()),
        };

        let readers = params.readers()?;
        // Collect all segments into a Vec since we need to verify properties across all segments
        let segments = readers.segments(1000, 100)?.collect::<Vec<_>>();

        assert!(!segments.is_empty(), "Should have at least one segment");

        // Check segment properties
        for segment in segments {
            assert_eq!(segment.region.contig, "chr19");
            assert!(segment.region.start >= 6105700);
            assert!(segment.region.end <= 6105800);
        }

        Ok(())
    }

    #[test]
    fn test_overlapping_segments() -> Result<()> {
        let params = ReaderParams {
            bam_file: get_test_bam(),
            fasta_file: get_test_fasta(),
            region: Some("chr19:6105700-6105900".parse().unwrap()),
        };

        let mut readers = params.readers()?;
        let segment_max_length = 100; // Small max length to force multiple segments
        let segment_overlap = 20;
        let segments = readers.segments(segment_max_length, segment_overlap)?.collect::<Vec<_>>();

        assert!(segments.len() > 1, "Should have multiple segments");

        // Check overlaps between adjacent segments
        for pair in segments.windows(2) {
            let first = readers.segment(&pair[0], 2)?;
            let second = readers.segment(&pair[1], 2)?;

            // Verify overlap amount
            let overlap = first.range.region.end - second.range.region.start;
            assert_eq!(overlap, segment_overlap, "Overlap amount should match configured value");

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

    #[test]
    fn test_segment_to_fetch_definition() -> Result<()> {
        // Create a test segment
        let region = ChunkRegion {
            region: Region { contig: "chr19".into(), start: 6105700, end: 6105800 },
            last_position: 6105900,
        };
        let segment = Segment { range: region.clone(), sequence: vec![65, 66, 67] }; // "ABC"

        // Convert to FetchDefinition
        let fetch_def = FetchDefinition::try_from(&segment)?;

        // Verify conversion worked correctly
        if let FetchDefinition::RegionString(chr, start, end) = fetch_def {
            assert_eq!(chr, b"chr19");
            assert_eq!(start, 6105700);
            assert_eq!(end, 6105800);
        } else {
            panic!("Unexpected FetchDefinition variant");
        }

        Ok(())
    }

    #[test]
    fn test_get_selected_region_variations() -> Result<()> {
        let params =
            ReaderParams { bam_file: get_test_bam(), fasta_file: get_test_fasta(), region: None };

        let readers = params.readers()?;
        let header = readers.bam.header();

        // Test chromosome-only region
        let region_chr_only: RegionString = "chr19".parse().unwrap();
        let full_region = get_selected_region(&region_chr_only, header)?;

        assert_eq!(full_region.0.contig, "chr19");
        assert_eq!(full_region.0.start, 1); // Should default to 1

        // The end should be the chromosome length from the header
        let chr19_tid = header.tid(b"chr19").unwrap();
        let chr19_len = header.target_len(chr19_tid).unwrap();
        assert_eq!(full_region.0.end, chr19_len);

        // Test chromosome with start but no end
        let region_with_start: RegionString = "chr19:100".parse().unwrap();
        let full_region = get_selected_region(&region_with_start, header)?;

        assert_eq!(full_region.0.contig, "chr19");
        assert_eq!(full_region.0.start, 100);
        assert_eq!(full_region.0.end, chr19_len); // Should default to chromosome length

        Ok(())
    }

    #[test]
    fn test_get_selected_region_errors() -> Result<()> {
        let params =
            ReaderParams { bam_file: get_test_bam(), fasta_file: get_test_fasta(), region: None };

        let readers = params.readers()?;
        let header = readers.bam.header();

        // Test non-existent chromosome
        let region_invalid_chr: RegionString = "nonexistent".parse().unwrap();
        let result = get_selected_region(&region_invalid_chr, header);
        assert!(result.is_err());

        // Get valid chromosome length
        let chr19_tid = header.tid(b"chr19").unwrap();
        let chr19_len = header.target_len(chr19_tid).unwrap();

        // Test start position beyond chromosome length
        let invalid_start = chr19_len + 100;
        let region_invalid_start: RegionString = format!("chr19:{invalid_start}").parse().unwrap();
        let result = get_selected_region(&region_invalid_start, header);
        assert!(result.is_err());

        // Test end position beyond chromosome length
        let invalid_end = chr19_len + 100;
        let region_invalid_end: RegionString = format!("chr19:100-{invalid_end}").parse().unwrap();
        let result = get_selected_region(&region_invalid_end, header);
        assert!(result.is_err());

        Ok(())
    }

    #[test]
    fn test_get_full_regions() -> Result<()> {
        let params =
            ReaderParams { bam_file: get_test_bam(), fasta_file: get_test_fasta(), region: None };

        let readers = params.readers()?;
        let header = readers.bam.header();

        let full_regions = get_full_regions(header)?;

        // Should have at least one region
        assert!(!full_regions.is_empty());

        // Verify that chromosome names match what's in the BAM header
        for (i, region) in full_regions.iter().enumerate() {
            let target_name = std::str::from_utf8(header.target_names()[i]).unwrap();
            assert_eq!(region.0.contig, target_name);

            // Start should be 1 (1-based)
            assert_eq!(region.0.start, 1);

            // End should match the chromosome length
            let tid = u32::try_from(i).unwrap();
            let chr_len = header.target_len(tid).unwrap();
            assert_eq!(region.0.end, chr_len);
        }

        Ok(())
    }

    #[test]
    fn test_segment_error_handling() -> Result<()> {
        let params =
            ReaderParams { bam_file: get_test_bam(), fasta_file: get_test_fasta(), region: None };

        let mut readers = params.readers()?;

        // Test with an invalid region (non-existent chromosome)
        let invalid_region = ChunkRegion {
            region: Region { contig: "nonexistent".into(), start: 100, end: 200 },
            last_position: 300,
        };

        let result = readers.segment(&invalid_region, 2);
        assert!(result.is_err());

        Ok(())
    }

    #[test]
    fn test_non_overlapping_segments() -> Result<()> {
        let params = ReaderParams {
            bam_file: get_test_bam(),
            fasta_file: get_test_fasta(),
            region: Some("chr19:6105700-6105900".parse().unwrap()),
        };

        let readers = params.readers()?;

        let segment_max_length = 50; // Small max length to force multiple segments
        let segment_overlap = 0; // No overlap
        let segments = readers.segments(segment_max_length, segment_overlap)?.collect::<Vec<_>>();

        assert!(segments.len() > 1, "Should have multiple segments");

        // Check no overlaps between adjacent segments
        for pair in segments.windows(2) {
            // The first segment should end exactly where the next one starts
            assert_eq!(pair[0].region.end, pair[1].region.start);
        }

        Ok(())
    }

    #[test]
    fn test_sequence_slice() -> Result<()> {
        let segment = Segment {
            range: ChunkRegion {
                region: Region { contig: "chr19".into(), start: 6105700, end: 6105800 },
                last_position: 6105900,
            },
            sequence: b"ATCGG".into(),
        };

        // Test valid slice
        let slice = segment.sequence_slice::<5>(1, 4)?;
        assert_eq!(slice.as_slice(), &[Base::T, Base::C, Base::G]);

        // Test out-of-bounds slice
        let slice = segment.sequence_slice::<5>(0, 10)?;
        assert_eq!(slice.as_slice(), &[Base::A, Base::T, Base::C, Base::G, Base::G]);

        Ok(())
    }

    proptest!(
        #[test]
        fn proptest_sequence_slice(start in 0usize..100, end in 0usize..100, seq in "[ATCG]{0,10}") {
            prop_assume!(start <= end, "Start must be less than or equal to end");
            let segment = Segment {
                range: ChunkRegion {
                    region: Region { contig: "chr19".into(), start: 6105700, end: 6105800 },
                    last_position: 6105900,
                },
                sequence: seq.into_bytes(),
            };

            // Ensure we don't panic on out-of-bounds slices
            segment.sequence_slice::<5>(start, end).unwrap();
        }
    );
}
