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
use color_eyre::{Result, eyre::bail};
use rust_htslib::bam::{self, FetchDefinition, Read as _};
use smol_str::SmolStr;
use tracing::{debug, instrument, trace};

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

    /// Number of threads to use for reading the bam file
    #[arg(long, default_value_t = 4)]
    pub threads: usize,
}

impl SegmentsParams {
    pub fn readers(&self) -> Result<Readers> {
        let fasta = open_fasta(&self.fasta_file)?;

        let mut bam = bam::IndexedReader::from_path(self.bam_file.path())?;
        bam.set_threads(self.threads)?;

        Ok(Readers { fasta, bam })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct FullRegion {
    /// The chromosome name, e.g. "chr19"
    pub chromosome: SmolStr,
    /// The start position of the region, inclusive
    pub start: u64,
    /// The end position of the region, exclusive
    pub end: u64,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ChunkRegion {
    /// The chromosome name, e.g. "chr19"
    pub chromosome: SmolStr,
    /// The start position of the region, inclusive
    pub start: u64,
    /// The end position of the region, exclusive
    pub end: u64,
    /// The last base of the region, exclusive
    pub last_position: u64,
}

impl ChunkRegion {
    pub fn contains(&self, pos: u64) -> bool {
        (self.start..self.end).contains(&pos)
    }
}

pub struct Readers {
    pub fasta: FastaReader,
    pub bam: bam::IndexedReader,
}

#[derive(Debug)]
pub struct Segment {
    pub range: ChunkRegion,
    pub sequence: Vec<u8>,
}

impl<'seg> From<&'seg Segment> for FetchDefinition<'seg> {
    fn from(segment: &'seg Segment) -> Self {
        FetchDefinition::RegionString(
            segment.range.chromosome.as_bytes(),
            i64::try_from(segment.range.start).expect("start is valid i64"),
            i64::try_from(segment.range.end).expect("end is valid i64"),
        )
    }
}

impl Readers {
    /// Calculate segments based on configuration parameters
    #[instrument(level = "debug", skip_all)]
    pub fn segments(
        &self,
        params: &SegmentsParams,
    ) -> Result<impl Iterator<Item = ChunkRegion> + use<>> {
        let full_regions = if let Some(region) = &params.region {
            debug!(?region, "fetching specified region");
            let start = region.start.map(|x| x.get().into()).unwrap_or(1);
            let last_position = {
                // If no end specified, use chromosome length from BAM header
                let header = self.bam.header();
                header
                    .target_len(header.tid(region.chromosome.as_bytes()).expect("get tid"))
                    .expect("fetch header length")
            };
            let end = region.end.map(|x| x.get().into()).unwrap_or(last_position);

            // Since the user specified a region, let's go with only that
            vec![FullRegion { chromosome: region.chromosome.clone(), start, end }]
        } else {
            debug!("fetching all regions");
            get_full_regions(self.bam.header())
        };

        if full_regions.is_empty() {
            bail!("No regions found");
        }

        let initial_start = full_regions[0].start;
        let chunked = ChunkedRegions {
            full_regions,
            current_region_idx: 0,
            current_start: initial_start,
            max_length: params.max_segment_length,
            overlap: params.segment_overlap,
        };

        Ok(chunked)
    }

    pub fn segment(&mut self, region: &ChunkRegion) -> Result<Segment> {
        let fetch_some_more = 2;
        let last_position_to_fetch =
            region.end.wrapping_add(fetch_some_more).min(region.last_position);

        // Calculate exact capacity needed to avoid reallocations
        let len = usize::try_from(last_position_to_fetch.wrapping_sub(region.start))
            .expect("failed to convert segment length to usize");

        trace!(?region, len, "fetching segment");
        let mut seq = Vec::with_capacity(len);
        self.fasta.fetch(&region.chromosome, region.start, last_position_to_fetch)?;
        self.fasta.read(&mut seq)?;
        Ok(Segment { range: region.clone(), sequence: seq })
    }
}

struct ChunkedRegions {
    full_regions: Vec<FullRegion>,
    current_region_idx: usize,
    current_start: u64,
    max_length: u64,
    overlap: u64,
}

impl Iterator for ChunkedRegions {
    type Item = ChunkRegion;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.current_region_idx >= self.full_regions.len() {
                return None;
            }

            let full_region = &self.full_regions[self.current_region_idx];
            if self.current_start >= full_region.end {
                // Move to next region when we've finished the current one
                self.current_region_idx += 1;
                if self.current_region_idx < self.full_regions.len() {
                    self.current_start = self.full_regions[self.current_region_idx].start;
                }
                continue;
            }

            let end = self.current_start.saturating_add(self.max_length).min(full_region.end);
            let chunk = ChunkRegion {
                chromosome: full_region.chromosome.clone(),
                start: self.current_start,
                end,
                last_position: full_region.end,
            };

            self.current_start = end;
            if self.current_start < full_region.end {
                self.current_start = self.current_start.saturating_sub(self.overlap);
            }

            return Some(chunk);
        }
    }
}

#[instrument(level = "debug", skip(header))]
fn get_full_regions(header: &bam::HeaderView) -> Vec<FullRegion> {
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

            FullRegion {
                chromosome: chr,
                start: 1, // 1-based coordinates
                end: length,
            }
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

    // proptest! {
    //     #[test]
    //     fn test_chunk_up_properties(
    //         // Generate chromosome names
    //         chrom in "[A-Za-z0-9]{1,10}",
    //         // Generate reasonable region lengths
    //         start in 1u64..100_000u64,
    //         length in 1u64..100_000u64,
    //         // Generate reasonable chunk parameters
    //         max_length in 10_000u64..100_000u64,
    //         overlap in 1u64..1000u64
    //     ) {
    //         let end = start.saturating_add(length);
    //         let region = Region {
    //             chromosome: chrom.into(),
    //             start,
    //             end,
    //         };

    //         let params = SegmentsParams {
    //             bam_file: get_test_bam(),
    //             fasta_file: get_test_fasta(),
    //             region: None,
    //             max_segment_length: max_length,
    //             segment_overlap: overlap,
    //         };

    //         // Generate chunks using the iterator
    //         let chunks: Vec<Region> = chunk_up_iter(vec![region.clone()], &params).collect();

    //         // Properties that should hold for all chunk configurations:
    //         prop_assert!(!chunks.is_empty(), "Chunks should not be empty");

    //         // First chunk should start at original start
    //         prop_assert_eq!(chunks[0].start, region.start, "First chunk should start at region start");

    //         // Last chunk should end at original end
    //         prop_assert_eq!(chunks.last().unwrap().end, region.end, "Last chunk should end at region end");

    //         // All chunks should have same chromosome
    //         prop_assert!(chunks.iter().all(|c| c.chromosome == region.chromosome), "All chunks should have same chromosome");

    //         // No chunk should exceed max length
    //         prop_assert!(chunks.iter().all(|c| c.end - c.start <= max_length), "No chunk should exceed max length");

    //         // All chunks except last should be max_length unless they would exceed region end
    //         prop_assert!(chunks[..chunks.len()-1].iter().all(|c| {
    //             c.end - c.start == max_length || c.end == region.end
    //         }), "All chunks except last should be max_length or reach region end");

    //         // Adjacent chunks should overlap by segment_overlap (if not the last chunk)
    //         prop_assert!(chunks.windows(2).all(|w| {
    //             w[1].start == w[0].end.saturating_sub(overlap)
    //         }), "Adjacent chunks should have correct overlap");

    //         // Chunks should cover entire region
    //         for (i, chunk) in chunks.windows(2).enumerate() {
    //             prop_assert!(chunk[0].end >= chunk[1].start,
    //                 "Gap between chunks {} and {}", i, i+1);
    //         }
    //     }
    // }

    #[test]
    fn test_segment_reading() -> Result<()> {
        // Create a test region and params
        let region = ChunkRegion {
            chromosome: "chr19".into(),
            start: 6105700,
            end: 6105800,
            last_position: 6105900,
        };

        let params = SegmentsParams {
            bam_file: get_test_bam(),
            fasta_file: get_test_fasta(),
            region: None,
            max_segment_length: 1000,
            segment_overlap: 100,
            threads: 4,
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
            max_segment_length: 1000,
            segment_overlap: 100,
            threads: 4,
        };

        let readers = params.readers()?;
        // Collect all segments into a Vec since we need to verify properties across all segments
        let segments = readers.segments(&params)?.collect::<Vec<_>>();

        assert!(!segments.is_empty(), "Should have at least one segment");

        // Check segment properties
        for segment in segments {
            assert_eq!(segment.chromosome, "chr19");
            assert!(segment.start >= 6105700);
            assert!(segment.end <= 6105800);
        }

        Ok(())
    }

    #[test]
    fn test_overlapping_segments() -> Result<()> {
        let params = SegmentsParams {
            bam_file: get_test_bam(),
            fasta_file: get_test_fasta(),
            region: Some("chr19:6105700-6105900".parse().unwrap()),
            max_segment_length: 100, // Small max length to force multiple segments
            segment_overlap: 20,     // Known overlap amount
            threads: 4,
        };

        let readers = params.readers()?;
        let segments = readers.segments(&params)?.collect::<Vec<_>>();

        assert!(segments.len() > 1, "Should have multiple segments");

        // Check overlaps between adjacent segments
        // for pair in segments.windows(2) {
        //     let first = readers.segment(&pair[0])?;
        //     let second = readers.segment(&pair[1])?;

        //     // Verify overlap amount
        //     let overlap = first.range.end - second.range.start;
        //     assert_eq!(overlap, params.segment_overlap);

        //     // Verify overlapping sequence content matches
        //     let overlap_first = &first.sequence[first.sequence.len() - overlap as usize..];
        //     let overlap_second = &second.sequence[..overlap as usize];
        //     assert_eq!(overlap_first, overlap_second, "Overlapping sequences should match");
        // }

        Ok(())
    }
}
