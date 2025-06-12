//! Helpers for writing concise tests.

use crate::{
    call::variants::VariantCandidatePileup,
    sequence::{SegmentationParams, SegmentsParams},
    utils::RegionString,
};
use clio::ClioPath;
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use std::num::NonZeroU32;

/// Fetch a pileup for a specific variant position in the test BAM file.
pub fn variant_pileup(chr: &str, start: u32) -> Result<VariantCandidatePileup> {
    let region = RegionString {
        chromosome: chr.into(),
        start: Some(NonZeroU32::new(start.saturating_sub(100).max(1)).unwrap()),
        end: None,
    };
    let p = SegmentsParams {
        bam_file: ClioPath::new("tests/data/test.bam").unwrap(),
        fasta_file: ClioPath::new("tests/data/test.fasta.gz").unwrap(),
        region: Some(region.clone()),
        threads: 1,
        segmentation: SegmentationParams { segment_max_length: 1000, segment_overlap: 0 },
    };
    let mut readers = p.readers().wrap_err("failed to fetch segments")?;

    let chunk = readers.segments()?.next().wrap_err("failed to fetch segment")?;

    let pileups = chunk.process(&mut readers).wrap_err("failed to process region")?;
    let pileup = pileups
        .into_iter()
        .find(|p| p.pos == start)
        .ok_or_else(|| color_eyre::eyre::eyre!("No pileup found at position {}", start))?;

    Ok(pileup)
}
