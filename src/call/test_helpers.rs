//! Helpers for writing concise tests.
#![allow(unused, clippy::unwrap_in_result, reason = "test code")]
#![cfg_attr(coverage_nightly, coverage(off))]

#[cfg(not(feature = "experimental-seqair"))]
use crate::call::pileup::Pileup;
use crate::{
    call::{
        process::{PileupMappingParams, get_pileups},
        variant_calling::{ReadFlags, VariantCallingParams},
    },
    metrics::PileupMetrics,
    sequence::{PileupReaders, ReaderParams, Readers, Segment},
    utils::{CliRegionInput, RegionString},
};
use clio::ClioPath;
use color_eyre::{
    Section,
    eyre::{ContextCompat as _, Result, WrapErr, eyre},
};
use seqair_types::Pos1;
use std::rc::Rc;

impl ReaderParams {
    pub fn test_data() -> Self {
        Self {
            bam_file: ClioPath::new("tests/data/test.bam").unwrap(),
            fasta_file: ClioPath::new("tests/data/test.fasta.gz").unwrap(),
            regions: None,
        }
    }

    pub fn around(&mut self, chr: &str, pos: u32) -> &mut Self {
        let region = RegionString {
            chromosome: chr.into(),
            start: Some(Pos1::new(pos.saturating_sub(60).max(1)).unwrap()),
            end: None,
        };
        self.regions = Some(CliRegionInput::from_region(region));
        self
    }

    pub fn test_with(bam: &str, fasta: &str) -> Self {
        Self {
            bam_file: ClioPath::new(bam).unwrap(),
            fasta_file: ClioPath::new(fasta).unwrap(),
            regions: None,
        }
    }

    #[cfg(not(feature = "experimental-seqair"))]
    pub fn pileup(&self, chr: &str, pos: u32) -> Result<(Rc<Segment>, Pileup)> {
        let region = RegionString {
            chromosome: chr.into(),
            start: Some(Pos1::new(pos.saturating_sub(60).max(1)).unwrap()),
            end: None,
        };
        let params = Self { regions: Some(CliRegionInput::from_region(region)), ..self.clone() };
        let mut readers = params.pileup_readers().wrap_err("failed to fetch segments")?;
        let chunk = readers.segments(1000, 0)?.next().wrap_err("failed to fetch segment")?;

        let pileup_mapping_params = PileupMappingParams::default();
        let (segment, pileups) = get_pileups(&mut readers, &chunk, &pileup_mapping_params)
            .wrap_err("failed to process region")?;
        let pileup = pileups
            .into_iter()
            .find(|p| p.pos == pos)
            .ok_or_else(|| eyre!("No variant at {chr}:{pos}"))
            .note(
                "Variant pileups are only built when at least one base differs from the reference",
            )?;

        Ok((segment, pileup))
    }

    #[cfg(feature = "experimental-seqair")]
    pub fn pileup(&self, chr: &str, pos: u32) -> Result<(Rc<Segment>, PileupMetrics)> {
        let region = RegionString {
            chromosome: chr.into(),
            start: Some(Pos1::new(pos.saturating_sub(60).max(1)).unwrap()),
            end: None,
        };
        let params = Self { regions: Some(CliRegionInput::from_region(region)), ..self.clone() };
        let mut readers = params.pileup_readers().wrap_err("failed to fetch segments")?;
        let chunk = readers.segments(1000, 0)?.next().wrap_err("failed to fetch segment")?;

        let pileup_mapping_params = PileupMappingParams::default();
        let (segment, pileups) = get_pileups(&mut readers, &chunk, &pileup_mapping_params)
            .wrap_err("failed to process region")?;
        let pileup = pileups
            .into_iter()
            .find(|p| p.pos == pos)
            .ok_or_else(|| eyre!("No variant at {chr}:{pos}"))
            .note(
                "Variant pileups are only built when at least one base differs from the reference",
            )?;

        Ok((segment, pileup))
    }
}

pub(crate) fn test_readers(chr: &str, pos: u32) -> Result<Readers> {
    let region = RegionString {
        chromosome: chr.into(),
        // Make sure to fetch some context around the position for metrics
        start: Some(Pos1::new(pos.saturating_sub(60).max(1)).unwrap()),
        end: None,
    };
    let p = ReaderParams {
        bam_file: ClioPath::new("tests/data/test.bam").unwrap(),
        fasta_file: ClioPath::new("tests/data/test.fasta.gz").unwrap(),
        regions: Some(CliRegionInput::from_region(region)),
    };
    p.readers().wrap_err("failed to fetch segments")
}

/// Fetch a pileup for a specific variant position (0-based!) in the test BAM file.
///
/// When comparing this to IGV, please keep in mind that IGV and VCF files use
/// 1-based positions, so the `pos` parameter is off by one compared to what you
/// see there.
#[cfg(not(feature = "experimental-seqair"))]
pub(crate) fn variant_pileup(chr: &str, pos: u32) -> Result<(Rc<Segment>, Pileup)> {
    ReaderParams::test_data().pileup(chr, pos)
}

#[cfg(feature = "experimental-seqair")]
pub(crate) fn variant_pileup(chr: &str, pos: u32) -> Result<(Rc<Segment>, PileupMetrics)> {
    ReaderParams::test_data().pileup(chr, pos)
}
