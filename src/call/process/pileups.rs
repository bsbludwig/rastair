use crate::{
    call::{
        pileup::{
            Pileup,
            from_hts::{ReadMismatchCache, ReadOrientationCache},
            overlapping_reads::NameCollector,
        },
        require_tags::TagRequirement,
        variant_calling::VariantCallingParams,
    },
    sequence::{ChunkRegion, Readers, Segment},
    utils::logging::ThisIsABug,
};
use color_eyre::eyre::{Result, WrapErr};
use rust_htslib::bam::{FetchDefinition, Read as _};
use std::{ops::Deref, rc::Rc};
use tracing::{Level, debug, instrument, trace, warn};

#[derive(better_default::Default)]
pub struct PileupMappingParams {
    pub variant_calling: VariantCallingParams,
    pub require_tags: TagRequirement,
    /// Ignore indels within this distance from read ends.
    #[default(20)]
    pub indel_end_of_read_cutoff: usize,
    /// Homopolymer repeat length cutoff for indel depth filtering.
    #[default(3)]
    pub indel_repeat_limit: usize,
    /// Maximum number of non-TAPS mismatches allowed on a read supporting an indel.
    #[default(5)]
    pub indel_max_mismatches: u32,
}

impl Deref for PileupMappingParams {
    type Target = VariantCallingParams;

    fn deref(&self) -> &Self::Target {
        &self.variant_calling
    }
}

#[instrument(level = "info", skip_all)]
pub fn get_pileups(
    readers: &mut Readers,
    region: &ChunkRegion,
    params: &PileupMappingParams,
) -> Result<(Rc<Segment>, impl Iterator<Item = Pileup>)> {
    let segment = readers.segment(region, 2).wrap_err("Failed to fetch segment")?;
    debug!(len = segment.sequence.len(), "Processing region");

    // Fetch the pileups for the segment
    FetchDefinition::try_from(&segment.region)
        .wrap_err("Could not convert region string")
        .this_is_a_bug()
        .and_then(|r| readers.bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))
        .wrap_err_with(|| format!("Could not fetch region `{}` from BAM file", region.region))?;

    // Install a pileup-level filter so that reads can be rejected once at read
    // time instead of once per pileup column they span.
    {
        let read_flags = params.read_flags.clone();
        let unpaired = params.unpaired;
        let tag_filter = params.require_tags.clone();
        readers.bam.set_pileup_filter(move |record| {
            read_flags.filter_flags(record.flags(), unpaired) && tag_filter.allows(&record)
        });
    }

    let segment = Rc::new(segment);
    let segment_clone = segment.clone();

    // Go over each column in the pileup from htslib and build our own pileup
    let mut pileup = readers.bam.pileup();
    pileup.set_max_depth(params.max_coverage);
    let mut collector = NameCollector::new(params);
    let mut orientation_cache = ReadOrientationCache::default();
    let mut mismatch_cache = ReadMismatchCache::default();
    let piles = pileup
        .filter_map(|p| match p {
            Ok(p) => Some(p),
            Err(e) => {
                if tracing::enabled!(Level::TRACE) {
                    trace!(%e, "Failed to read pileup, skipping");
                }
                None
            }
        })
        .filter(|p| {
            // We might get pileups from htslib that are not in the region of
            // interest (but actually before it). Since our segments only cover
            // the specified region, we can just skip these (we won't have the
            // reference sequence for this anyway). They'll be part of the
            // next/previous segment anyway.
            region.contains(u64::from(p.pos()))
        })
        .map(move |pile| {
            Pileup::from_hts(
                &pile,
                segment.clone(),
                params,
                &mut collector,
                &mut orientation_cache,
                &mut mismatch_cache,
            )
            .wrap_err_with(|| {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::ReaderParams;

    #[test]
    fn test_reading_bounds() -> Result<()> {
        // check that we can read exactly the right positions
        let params = ReaderParams {
            region: Some("chr19:6105700-6105800".parse().unwrap()),
            ..ReaderParams::test_data()
        };
        let mut readers = params.readers()?;
        let segments: Vec<_> = readers.segments(10_000, 100)?.collect();
        readers.segment(&segments[0], 0)?;

        let pileup_mapping_params = PileupMappingParams::default();
        let (_segment, pileups) = get_pileups(&mut readers, &segments[0], &pileup_mapping_params)?;
        let pileups: Vec<_> = pileups.collect();

        assert!(!pileups.is_empty());
        assert_eq!(pileups.first().unwrap().pos, 6_105_700);
        assert_eq!(pileups.last().unwrap().pos, 6_105_800);

        Ok(())
    }
}
