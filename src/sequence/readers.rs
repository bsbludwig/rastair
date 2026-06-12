use crate::{
    sequence::{
        ChunkRegion, Region, SelectedRegion, chunked::ChunkedRegions, segementation::Segment,
    },
    utils::{
        CliRegionInput, RegionString, cli,
        file_helpers::{FastaReader, open_fasta},
    },
};
use clap::value_parser;
use clio::ClioPath;
use color_eyre::{
    Result, Section,
    eyre::{Context, ContextCompat, ensure},
};
use rust_htslib::bam::{self, HeaderView, Read as _};
use seqair_types::SmolStr;
use tracing::{debug, instrument, trace};

#[derive(Debug, clap::Args, Clone)]
pub struct ReaderParams {
    /// Path to sorted and indexed BAM or CRAM file
    #[arg(value_parser=value_parser!(ClioPath).exists().is_file(), value_hint=clap::ValueHint::FilePath)]
    #[arg(help_heading = cli::sections::INPUT)]
    pub bam_file: ClioPath,

    /// Path to sorted and indexed (via samtools faidx) FASTA file. Can be bgzip
    /// compressed, but requires both a gzi index and a fai index
    #[arg(short='r', long, value_parser=value_parser!(ClioPath).exists().is_file(), value_hint=clap::ValueHint::FilePath)]
    #[arg(help_heading = cli::sections::INPUT)]
    pub fasta_file: ClioPath,

    /// Restrict processing to specific genomic regions.
    ///
    /// Accepts either space-separated region strings or a single BED file:
    /// - Region strings: `chr`, `chr:start`, `chr:start-end` (1-based inclusive)
    /// - Multiple regions separated by whitespace: `"chr1 chr2:100-200"`
    /// - BED files with `@` prefix: `@targets.bed`
    #[arg(short = 'l', long = "region", value_parser = clap::value_parser!(CliRegionInput))]
    #[arg(help_heading = cli::sections::INPUT)]
    pub regions: Option<CliRegionInput>,
}

impl ReaderParams {
    pub fn readers(&self) -> Result<Readers> {
        let fasta = open_fasta(&self.fasta_file)?;
        let bam_path = self.bam_file.path();
        let mut bam = bam::IndexedReader::from_path(bam_path)
            .with_suggestion(|| {
                format!(
                    "Ensure the BAM/CRAM file is sorted and indexed with \
                    `samtools sort {bam_path:?}` and `samtools index {bam_path:?}`, respectively."
                )
            })
            .note("If you have a .bai/.crai file, ensure it is in the same directory as the BAM/CRAM file.")?;
        bam.set_reference(self.fasta_file.path())
            .wrap_err("Failed to set FASTA reference for BAM reader")
            .note("Rastair itself already opened the FASTA file successfully, this error is from the BAM reader implementation")?;

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
        let full_regions = if let Some(input) = &self.params.regions {
            input
                .regions()
                .iter()
                .map(|r| {
                    get_selected_region(r, self.bam.header())
                        .wrap_err("Failed to get selected region from BAM file")
                })
                .collect::<Result<Vec<_>>>()?
        } else {
            debug!("fetching all regions");
            get_full_regions(self.bam.header())
                .wrap_err("Failed to get all regions from BAM file")?
        };
        ensure!(!full_regions.is_empty(), "No regions found");

        let initial_start = full_regions[0].start;
        let chunked = ChunkedRegions {
            full_regions,
            current_region_idx: 0,
            current_start: initial_start,
            max_length: segment_max_length,
            overlap: segment_overlap,
        };

        Ok(chunked)
    }

    /// Fetch a segment from the FASTA file, with optional overfetching
    #[instrument(level = "debug", skip_all)]
    pub fn segment(&mut self, region: &ChunkRegion, overfetch: u64) -> Result<Segment> {
        let last_position_to_fetch = region.end.wrapping_add(overfetch).min(region.last_position);

        // Calculate exact capacity needed to avoid reallocations
        let len = usize::try_from(last_position_to_fetch.wrapping_sub(region.start))
            .wrap_err("Failed to convert segment length to usize")?;

        trace!(?region, len, "fetching segment");
        let seq = self
            .fasta
            .fetch_seq(&region.contig, region.start, last_position_to_fetch)
            .wrap_err_with(|| format!("Failed to get region {} from FASTA file", region.region))?;

        Ok(Segment {
            range: region.clone(),
            sequence: seq,
            overlap_start: region.overlap_start,
            overlap_end: region.overlap_end,
        })
    }
}

#[instrument(level = "debug", skip(bam_header))]
fn get_selected_region(region: &RegionString, bam_header: &HeaderView) -> Result<SelectedRegion> {
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
    Ok(SelectedRegion::UserDefinedSubset {
        region: Region { contig: region.chromosome.clone(), start, end },
        last_position,
    })
}

fn to_u64(value: seqair_types::Pos1) -> u64 {
    value.as_u64()
}

#[instrument(level = "debug", skip(header))]
fn get_full_regions(header: &bam::HeaderView) -> Result<Vec<SelectedRegion>> {
    header
        .target_names()
        .iter()
        .enumerate()
        .filter(|(_, name)| !name.is_empty())
        .map(|(tid, name)| -> Result<SelectedRegion> {
            let chr = SmolStr::new(
                std::str::from_utf8(name).wrap_err("BAM target name not valid UTF-8")
                    .note("This is against the BAM specification, please check with the tool that created this file")?,
            );
            let length = header
                .target_len(u32::try_from(tid).wrap_err("Failed to get a target ID that was part of the BAM header")
                    .note("The BAM header might be corrupt")?)
                .wrap_err("Failed to get target length")?;

            Ok(SelectedRegion::EntireContig(Region {
                contig: chr,
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
    fn test_get_selected_region_variations() -> Result<()> {
        let params =
            ReaderParams { bam_file: get_test_bam(), fasta_file: get_test_fasta(), regions: None };

        let readers = params.readers()?;
        let header = readers.bam.header();

        // Test chromosome-only region
        let region_chr_only: RegionString = "chr19".parse().unwrap();
        let full_region = get_selected_region(&region_chr_only, header)?;

        assert_eq!(full_region.contig, "chr19");
        assert_eq!(full_region.start, 1); // Should default to 1

        // The end should be the chromosome length from the header
        let chr19_tid = header.tid(b"chr19").unwrap();
        let chr19_len = header.target_len(chr19_tid).unwrap();
        assert_eq!(full_region.end, chr19_len);

        // Test chromosome with start but no end
        let region_with_start: RegionString = "chr19:100".parse().unwrap();
        let full_region = get_selected_region(&region_with_start, header)?;

        assert_eq!(full_region.contig, "chr19");
        assert_eq!(full_region.start, 100);
        assert_eq!(full_region.end, chr19_len); // Should default to chromosome length

        Ok(())
    }

    #[test]
    fn test_get_selected_region_errors() -> Result<()> {
        let params =
            ReaderParams { bam_file: get_test_bam(), fasta_file: get_test_fasta(), regions: None };

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
            ReaderParams { bam_file: get_test_bam(), fasta_file: get_test_fasta(), regions: None };

        let readers = params.readers()?;
        let header = readers.bam.header();

        let full_regions = get_full_regions(header)?;

        // Should have at least one region
        assert!(!full_regions.is_empty());

        // Verify that chromosome names match what's in the BAM header
        for (i, region) in full_regions.iter().enumerate() {
            let target_name = std::str::from_utf8(header.target_names()[i]).unwrap();
            assert_eq!(region.contig, target_name);

            // Start should be 1 (1-based)
            assert_eq!(region.start, 1);

            // End should match the chromosome length
            let tid = u32::try_from(i).unwrap();
            let chr_len = header.target_len(tid).unwrap();
            assert_eq!(region.end, chr_len);
        }

        Ok(())
    }
}
