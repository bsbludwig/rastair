//! Helpers for writing concise tests.

use crate::{
    call::variants::VariantCandidatePileup,
    sequence::{SegmentationParams, SegmentsParams},
    utils::RegionString,
};
use clio::ClioPath;
use color_eyre::{
    Section,
    eyre::{ContextCompat as _, Result, WrapErr},
};
use std::num::NonZeroU32;

/// Fetch a pileup for a specific variant position (0-based!) in the test BAM file.
///
/// When comparing this to IGV, please keep in mind that IGV and VCF files use
/// 1-based positions, so the `pos` parameter is off by one compared to what you
/// see there.
pub fn variant_pileup(chr: &str, pos: u32) -> Result<VariantCandidatePileup> {
    let region = RegionString {
        chromosome: chr.into(),
        // Make sure to fetch some context around the position for metrics
        start: Some(NonZeroU32::new(pos.saturating_sub(60).max(1)).unwrap()),
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
        .find(|p| p.pos == pos)
        .ok_or_else(|| color_eyre::eyre::eyre!("No variant at {chr}:{pos}"))
        .note("Variant pileups are only built when at least one base differs from the reference")?;

    Ok(pileup)
}
