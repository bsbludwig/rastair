//! Helpers for writing concise tests.

use crate::{
    call::{
        process::{IncludeAllCpGs, PileupMappingParams, pileup_mapper},
        variants::{SeenBases, VariantCandidatePileup},
    },
    sequence::{Readers, SegmentationParams, SegmentsParams},
    utils::RegionString,
};
use clio::ClioPath;
use color_eyre::{
    Section,
    eyre::{ContextCompat as _, Result, WrapErr},
};
use rust_htslib::bam::{FetchDefinition, Read as _, pileup::Pileup};
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
        threads: 1,
        segmentation: SegmentationParams { segment_max_length: 1000, segment_overlap: 0 },
    };
    p.readers().wrap_err("failed to fetch segments")
}

pub(crate) fn pileup(chr: &str, pos: u32) -> Result<Pileup> {
    let mut readers = test_readers(chr, pos)?;
    let chunk = readers.segments()?.next().wrap_err("failed to fetch segment")?;
    let segment = readers.segment(&chunk)?;
    FetchDefinition::try_from(&segment)
        .and_then(|r| readers.bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))?;

    readers
        .bam
        .pileup()
        .filter_map(|p| p.ok())
        .find(|p| p.pos() == pos)
        // .map(|pile| (SeenBases(pile.alignments().filter_map(pileup_mapper).collect()), pile))
        .wrap_err_with(|| format!("No pileup found at {chr}:{pos} in BAM file"))
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
