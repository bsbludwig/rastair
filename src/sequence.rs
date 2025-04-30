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
// segements are optional in the beginning

use crate::utils::{
    RegionString,
    file_helpers::{FastaReader, open_fasta},
};
use clap::value_parser;
use clio::ClioPath;
use color_eyre::Result;
use rust_htslib::bam::{self, FetchDefinition, Read as _};
use smallarcvec::SmallByteVec;
use smol_str::SmolStr;

#[derive(Debug, clap::Args)]
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

    /// Maximum length of a segment in bases
    #[arg(long, default_value_t = 1_000_000)]
    pub max_segment_length: u64,

    /// Number of bases to overlap between segments
    #[arg(long, default_value_t = 100)]
    pub segment_overlap: u64,
}

impl SegmentsParams {
    pub fn segments(&self) -> Result<Readers> {
        let fasta = open_fasta(&self.fasta_file)?;
        // indexed_reader.fetch("chr19", fetch_range.start, fetch_range.end + 1)?;

        let mut bam = bam::IndexedReader::from_path(self.bam_file.path())?;
        bam.set_threads(8)?;
        if let Some(region) = &self.region {
            bam.fetch(region)?;
        } else {
            bam.fetch(FetchDefinition::All)?;
        }
        // bam.fetch(("chr19", fetch_range.start, fetch_range.end + 1))?;

        Ok(Readers { fasta, bam })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Region {
    pub chromosome: SmolStr,
    pub start: u64,
    pub end: u64,
}

pub struct Readers {
    pub fasta: FastaReader,
    pub bam: bam::IndexedReader,
}

#[derive(Debug, Clone)]
pub struct Segment {
    pub range: Region,
    pub sequence: SmallByteVec,
}

impl Readers {
    /// Calculate segments based on configuration parameters
    pub fn calculate_segments(&mut self, params: &SegmentsParams) -> Result<Vec<Segment>> {
        let full_regions = if let Some(region) = &params.region {
            // If region is specified, start with that
            vec![Region {
                chromosome: region.chromosome.clone(),
                start: region.start.map(|s| u64::from(s.get())).unwrap_or(1),
                end: region.end.map(|e| u64::from(e.get())).unwrap_or_else(|| {
                    // If no end specified, use chromosome length from BAM header
                    let header = self.bam.header();
                    header
                        .target_len(header.tid(region.chromosome.as_bytes()).expect("get tid"))
                        .expect("fetch header length")
                }),
            }]
        } else {
            // If no region specified, create regions for all chromosomes
            get_full_regions(self.bam.header())
        };

        // Chunk up the regions into sub-regions that are at most max_segment_length and have an overlap of segment_overlap
        let chunked_regions = chunk_up(full_regions.as_slice(), params);
        // Create segments from the chunked regions
        chunked_regions.into_iter().map(|region| self.segment(&region)).collect::<Result<Vec<_>>>()
    }

    fn segment(&mut self, region: &Region) -> Result<Segment> {
        let mut seq = Vec::with_capacity(
            usize::try_from(region.end.wrapping_sub(region.start.wrapping_add(1)))
                .expect("failed to convert u64 to usize"),
        );
        self.fasta.fetch(&region.chromosome, region.start, region.end)?;
        self.fasta.read(&mut seq)?;
        Ok(Segment { range: region.clone(), sequence: seq.into() })
    }
}

fn get_full_regions(header: &bam::HeaderView) -> Vec<Region> {
    header
        .target_names()
        .iter()
        .enumerate()
        .filter(|(_, name)| !name.is_empty())
        .map(|(tid, name)| {
            let chr = SmolStr::new(
                std::str::from_utf8(name).expect("bam target name always valid UTF-8"),
            );
            let length =
                header.target_len(u32::try_from(tid).expect("get tid")).expect("get target length");

            Region {
                chromosome: chr,
                start: 1, // 1-based coordinates
                end: length,
            }
        })
        .collect()
}

fn chunk_up(full_regions: &[Region], params: &SegmentsParams) -> Vec<Region> {
    let mut chunked_regions = Vec::<Region>::with_capacity(full_regions.len() * 2); // lower-bound estimate
    for region in full_regions {
        let mut start = region.start;
        while start < region.end {
            let end = (start + params.max_segment_length).min(region.end);
            chunked_regions.push(Region { chromosome: region.chromosome.clone(), start, end });
            start = end.saturating_sub(params.segment_overlap);
        }
    }
    chunked_regions
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
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

    proptest! {
        #[test]
        fn test_chunk_up_properties(
            // Generate chromosome names
            chrom in "[A-Za-z0-9]{1,10}",
            // Generate reasonable region lengths
            start in 1u64..100_000u64,
            length in 1u64..100_000u64,
            // Generate reasonable chunk parameters
            max_length in 10_000u64..100_000u64,
            overlap in 1u64..1000u64
        ) {
            let end = start.saturating_add(length);
            let region = Region {
                chromosome: chrom.into(),
                start,
                end,
            };

            let params = SegmentsParams {
                bam_file: get_test_bam(),
                fasta_file: get_test_fasta(),
                region: None,
                max_segment_length: max_length,
                segment_overlap: overlap,
            };

            let chunks = chunk_up(&[region.clone()], &params);

            // Properties that should hold for all chunk configurations:
            prop_assert!(!chunks.is_empty(), "Chunks should not be empty");

            // First chunk should start at original start
            prop_assert_eq!(chunks[0].start, region.start, "First chunk should start at region start");

            // Last chunk should end at original end
            prop_assert_eq!(chunks.last().unwrap().end, region.end, "Last chunk should end at region end");

            // All chunks should have same chromosome
            prop_assert!(chunks.iter().all(|c| c.chromosome == region.chromosome), "All chunks should have same chromosome");

            // No chunk should exceed max length
            prop_assert!(chunks.iter().all(|c| c.end - c.start <= max_length), "No chunk should exceed max length");

            // All chunks except last should be max_length
            prop_assert!(chunks[..chunks.len()-1].iter().all(|c| c.end - c.start == max_length),
                        "All chunks except last should be max_length");

            // Adjacent chunks should overlap by segment_overlap (if not the last chunk)
            prop_assert!(chunks.windows(2).all(|w| {
                w[1].start == w[0].end.saturating_sub(overlap)
            }), "Adjacent chunks should have correct overlap");

            // Chunks should cover entire region
            for (i, chunk) in chunks.windows(2).enumerate() {
                prop_assert!(chunk[0].end >= chunk[1].start,
                    "Gap between chunks {} and {}", i, i+1);
            }
        }
    }

    #[test]
    fn test_segment_reading() -> Result<()> {
        // Create a test region and params
        let region = Region { chromosome: "chr19".into(), start: 6105700, end: 6105800 };

        let params = SegmentsParams {
            bam_file: get_test_bam(),
            fasta_file: get_test_fasta(),
            region: None,
            max_segment_length: 1000,
            segment_overlap: 100,
        };

        // Initialize readers
        let mut readers = params.segments()?;

        // Test reading a single segment
        let segment = readers.segment(&region)?;

        // Verify segment properties
        assert_eq!(segment.range, region);
        assert!(!segment.sequence.is_empty(), "Sequence should not be empty");
        let expected_len =
            usize::try_from(region.end - region.start).expect("region length should fit in usize");
        assert_eq!(segment.sequence.len(), expected_len);

        Ok(())
    }

    #[test]
    fn test_get_full_regions() -> Result<()> {
        // Initialize BAM reader to get header
        let bam = bam::IndexedReader::from_path(get_test_bam().path())?;
        let regions = get_full_regions(bam.header());

        // Verify regions
        assert!(!regions.is_empty(), "Should have at least one region");

        // Check that all regions are valid
        for region in &regions {
            assert!(!region.chromosome.is_empty(), "Chromosome name should not be empty");
            assert!(region.start > 0, "Start should be 1-based");
            assert!(region.end >= region.start, "End should be >= start");
        }

        Ok(())
    }

    #[test]
    fn test_calculate_segments() -> Result<()> {
        let params = SegmentsParams {
            bam_file: get_test_bam(),
            fasta_file: get_test_fasta(),
            region: Some("chr19:6105700-6105800".parse().unwrap()),
            max_segment_length: 1000,
            segment_overlap: 100,
        };

        let mut readers = params.segments()?;
        let segments = readers.calculate_segments(&params)?;

        assert!(!segments.is_empty(), "Should have at least one segment");

        // Check segment properties
        for segment in segments {
            assert_eq!(segment.range.chromosome, "chr19");
            assert!(segment.range.start >= 6105700);
            assert!(segment.range.end <= 6105800);
            let expected_len = usize::try_from(segment.range.end - segment.range.start)
                .expect("segment length should fit in usize");
            assert_eq!(segment.sequence.len(), expected_len);
        }

        Ok(())
    }
}
