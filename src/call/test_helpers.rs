//! Helpers for writing concise tests.

use crate::{
    call::{
        process::{IncludeAllCpGs, PileupMappingParams},
        variants::VariantCandidatePileup,
    },
    sequence::{Readers, SegmentationParams, SegmentsParams},
    utils::RegionString,
};
use clio::ClioPath;
use color_eyre::{
    Section,
    eyre::{ContextCompat as _, Result, WrapErr},
};
use std::num::NonZeroU32;

pub(crate) fn test_readers(chr: &str, pos: u32) -> Result<Readers> {
    let region = RegionString {
        chromosome: chr.into(),
        // Make sure to fetch some context around the position for metrics
        start: Some(NonZeroU32::new(pos.saturating_sub(60).max(1)).unwrap()),
        end: None,
    };
    let p = SegmentsParams {
        bam_file: ClioPath::new("tests/data/test.bam").unwrap(),
        fasta_file: ClioPath::new("tests/data/test.fasta.gz").unwrap(),
        region: Some(region),
        segmentation: SegmentationParams { segment_max_length: 1000, segment_overlap: 0 },
    };
    p.readers().wrap_err("failed to fetch segments")
}

/// Fetch a pileup for a specific variant position (0-based!) in the test BAM file.
///
/// When comparing this to IGV, please keep in mind that IGV and VCF files use
/// 1-based positions, so the `pos` parameter is off by one compared to what you
/// see there.
pub(crate) fn variant_pileup(chr: &str, pos: u32) -> Result<VariantCandidatePileup> {
    let mut readers = test_readers(chr, pos)?;

    let chunk = readers.segments()?.next().wrap_err("failed to fetch segment")?;

    let pileups = chunk
        .process(
            &mut readers,
            &PileupMappingParams {
                include_cpgs: IncludeAllCpGs::Yes,
                keep_overlapping_reads: false,
                read_masking: Default::default(),
                read_flags: Default::default(),
            },
        )
        .wrap_err("failed to process region")?;
    let pileup = pileups
        .into_iter()
        .find(|p| p.pos == pos)
        .ok_or_else(|| color_eyre::eyre::eyre!("No variant at {chr}:{pos}"))
        .note("Variant pileups are only built when at least one base differs from the reference")?;

    Ok(pileup)
}
