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
use color_eyre::{
    Result,
    eyre::{Context as _, ContextCompat as _},
};
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
            let header = self.bam.header();
            header
                .target_names()
                .iter()
                .enumerate()
                .filter(|(_, name)| !name.is_empty())
                .map(|(tid, name)| {
                    let chr = SmolStr::new(
                        std::str::from_utf8(name).wrap_err("bam target name is not valid UTF-8")?,
                    );
                    let length = header
                        .target_len(u32::try_from(tid).expect("get tid"))
                        .wrap_err("Failed to get target length for target we just read")?;

                    Ok(Region {
                        chromosome: chr,
                        start: 1, // 1-based coordinates
                        end: length,
                    })
                })
                .collect::<Result<Vec<Region>>>()?
        };

        // Chunk up the regions into sub-regions that are at most max_segment_length and have an overlap of segment_overlap
        let mut chunked_regions = Vec::<Region>::with_capacity(full_regions.len() * 2); // lower-bound estimate
        for region in full_regions {
            let mut start = region.start;
            while start < region.end {
                let end = (start + params.max_segment_length).min(region.end);
                chunked_regions.push(Region { chromosome: region.chromosome.clone(), start, end });
                start = end.saturating_sub(params.segment_overlap);
            }
        }
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
