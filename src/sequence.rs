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
    #[arg(value_name="BAM_FILE", value_parser=value_parser!(ClioPath).exists().is_file())]
    bam_file: ClioPath,

    /// A sorted and indexed (via samtools faidx) fasta file. Can be bgzip
    /// compressed, but requires both a gzi index and a fai index
    #[arg(short='r', long, value_name="FASTA_FILE", required=true, value_parser=value_parser!(ClioPath).exists().is_file())]
    fasta_file: ClioPath,

    /// Restrict to a specific chromosome or region of a chromosome. Format is
    /// "chr", "chr:start" or "chr:start-end", where start is 1-based and end is
    /// inclusive.
    #[arg(short = 'l', long)]
    region: Option<RegionString>,
}

impl SegmentsParams {
    pub fn segments(&self) -> Result<Segments> {
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

        Ok(Segments { fasta, bam })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Region {
    pub chromosome: SmolStr,
    pub start: u32,
    pub end: u32,
}

pub struct Segments {
    pub fasta: FastaReader,
    pub bam: bam::IndexedReader,
}

#[derive(Debug, Clone)]
pub struct Segment {
    pos: Region,
    sequence: SmallByteVec,
}

pub struct SegmentIterator {}

impl Iterator for SegmentIterator {
    type Item = Segment;

    fn next(&mut self) -> Option<Self::Item> {
        todo!()
    }
}
