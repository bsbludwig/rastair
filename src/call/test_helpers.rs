//! Helpers for writing concise tests.
#![allow(unused, clippy::unwrap_in_result)]

use crate::{
    call::{
        pileup::Pileup,
        process::{IncludeAllCpGs, PileupMappingParams},
        variant_calling::{ReadFlags, VariantCallingParams},
    },
    sequence::{ReaderParams, Readers},
    utils::RegionString,
};
use clio::ClioPath;
use color_eyre::{
    Section,
    eyre::{ContextCompat as _, Result, WrapErr, eyre},
};
use std::{num::NonZeroU32, rc::Rc};

impl ReaderParams {
    pub fn test_data() -> Self {
        Self {
            bam_file: ClioPath::new("tests/data/test.bam").unwrap(),
            fasta_file: ClioPath::new("tests/data/test.fasta.gz").unwrap(),
            region: None,
        }
    }

    pub fn test_with(bam: &str, fasta: &str) -> Self {
        Self {
            bam_file: ClioPath::new(bam).unwrap(),
            fasta_file: ClioPath::new(fasta).unwrap(),
            region: None,
        }
    }

    pub fn pileup(&self, chr: &str, pos: u32) -> Result<Pileup> {
        let region = RegionString {
            chromosome: chr.into(),
            start: Some(NonZeroU32::new(pos.saturating_sub(60).max(1)).unwrap()),
            end: None,
        };
        let params = Self { region: Some(region), ..self.clone() };
        let mut readers = params.readers().wrap_err("failed to fetch segments")?;
        let chunk = readers.segments(1000, 0)?.next().wrap_err("failed to fetch segment")?;

        let pileup_mapping_params = PileupMappingParams {
            include_cpgs: IncludeAllCpGs::Yes,
            variant_calling: VariantCallingParams::default(),
        };
        let (segment, pileups) = chunk
            .process(&mut readers, &pileup_mapping_params)
            .wrap_err("failed to process region")?;
        let pileup = pileups
            .into_iter()
            .find(|p| p.pos == pos)
            .ok_or_else(|| eyre!("No variant at {chr}:{pos}"))
            .note(
                "Variant pileups are only built when at least one base differs from the reference",
            )?;

        Ok(pileup)
    }
}

pub(crate) fn test_readers(chr: &str, pos: u32) -> Result<Readers> {
    let region = RegionString {
        chromosome: chr.into(),
        // Make sure to fetch some context around the position for metrics
        start: Some(NonZeroU32::new(pos.saturating_sub(60).max(1)).unwrap()),
        end: None,
    };
    let p = ReaderParams {
        bam_file: ClioPath::new("tests/data/test.bam").unwrap(),
        fasta_file: ClioPath::new("tests/data/test.fasta.gz").unwrap(),
        region: Some(region),
    };
    p.readers().wrap_err("failed to fetch segments")
}

/// Fetch a pileup for a specific variant position (0-based!) in the test BAM file.
///
/// When comparing this to IGV, please keep in mind that IGV and VCF files use
/// 1-based positions, so the `pos` parameter is off by one compared to what you
/// see there.
pub(crate) fn variant_pileup(chr: &str, pos: u32) -> Result<Pileup> {
    ReaderParams::test_data().pileup(chr, pos)
}
