use crate::{
    call::{pileup::Pileup, variant_calling::VariantCallingParams},
    sequence::{ChunkRegion, Readers, Segment},
};
use color_eyre::eyre::{Result, WrapErr};
use rust_htslib::bam::{FetchDefinition, Read as _};
use std::{ops::Deref, rc::Rc};
use tracing::{Level, instrument, trace, warn};

#[derive(Debug, Clone)]
pub struct PileupMappingParams {
    pub include_cpgs: IncludeAllCpGs,
    pub variant_calling: VariantCallingParams,
}

impl Deref for PileupMappingParams {
    type Target = VariantCallingParams;

    fn deref(&self) -> &Self::Target {
        &self.variant_calling
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IncludeAllCpGs {
    Yes,
    No,
}

impl Deref for IncludeAllCpGs {
    type Target = bool;

    fn deref(&self) -> &Self::Target {
        match self {
            IncludeAllCpGs::Yes => &true,
            IncludeAllCpGs::No => &false,
        }
    }
}

impl ChunkRegion {
    /// Process the chunk region to collect pileups
    ///
    /// # Returns
    /// - The segment corresponding to the chunk region
    /// - An iterator over pileups in the region
    #[instrument(level = "info", skip_all)]
    pub fn process(
        &self,
        readers: &mut Readers,
        params: &PileupMappingParams,
    ) -> Result<(Rc<Segment>, impl Iterator<Item = Pileup>)> {
        let segment = readers.segment(self, 2).wrap_err("failed to fetch segment")?;
        trace!(len = segment.sequence.len(), "Processing region");

        // Fetch the pileups for the segment
        FetchDefinition::try_from(&segment.region)
            .wrap_err("Could not convert region string")
            .and_then(|r| readers.bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))
            .wrap_err_with(|| format!("Could not fetch region `{}` from BAM file", self.region))?;

        let segment = Rc::new(segment);
        let segment_clone = segment.clone();

        // Go over each column in the pileup and collect variant candidates
        let piles = readers
            .bam
            .pileup()
            .filter_map(|p| {
                if tracing::enabled!(Level::TRACE) {
                    match p {
                        Ok(p) => Some(p),
                        Err(e) => {
                            trace!(%e, "Failed to read pileup, skipping");
                            None
                        }
                    }
                } else {
                    p.ok()
                }
            })
            .filter(|p| {
                // Filter out pileups that are not in the region of interest
                self.contains(u64::from(p.pos()))
            })
            .map(move |pile| {
                Pileup::from_hts(&pile, segment.clone(), params).wrap_err_with(|| {
                    format!("Failed to get candidate from pileup at position {}", pile.pos())
                })
            })
            .filter_map(|res| match res {
                Ok(x) => Some(x),
                Err(error) => {
                    warn!(error = format!("{error:#}"), "Failed to get pileup, skipping");
                    None
                }
            });
        Ok((segment_clone, piles))
    }
}
