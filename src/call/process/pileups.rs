use crate::{
    call::{require_tags::TagRequirement, variant_calling::VariantCallingParams},
    sequence::{ChunkRegion, Segment},
};
use color_eyre::eyre::Result;
use std::{ops::Deref, rc::Rc};

#[cfg(not(feature = "experimental-seqair"))]
use crate::{
    call::pileup::{Pileup, from_hts::PileupScratch},
    sequence::Readers,
    utils::logging::ThisIsABug,
};
#[cfg(not(feature = "experimental-seqair"))]
use color_eyre::eyre::WrapErr;
#[cfg(not(feature = "experimental-seqair"))]
use rust_htslib::bam::{FetchDefinition, Read as _};
#[cfg(not(feature = "experimental-seqair"))]
use tracing::{Level, debug, instrument, trace, warn};

#[cfg(feature = "experimental-seqair")]
use crate::{
    call::pileup::{Pileup, overlapping_reads::NameCollector},
    sequence::{PileupReaders, ReferenceWindow},
};
#[cfg(feature = "experimental-seqair")]
use color_eyre::eyre::{ContextCompat as _, WrapErr as _};
#[cfg(feature = "experimental-seqair")]
use seqair::reader::SegmentOptions;
#[cfg(feature = "experimental-seqair")]
use seqair_types::{Base, Pos0};
#[cfg(feature = "experimental-seqair")]
use std::{num::NonZeroU32, sync::Arc};
#[cfg(feature = "experimental-seqair")]
use tracing::{debug, instrument, warn};

#[derive(better_default::Default)]
pub struct PileupMappingParams {
    pub variant_calling: VariantCallingParams,
    pub require_tags: TagRequirement,
    /// Ignore indels within this distance from read ends.
    #[default(0)]
    pub indel_end_of_read_cutoff: usize,
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

#[cfg(not(feature = "experimental-seqair"))]
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
    let mut scratch = PileupScratch::new(params);
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
            Pileup::from_hts(&pile, segment.clone(), params, &mut scratch).wrap_err_with(|| {
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

#[cfg(feature = "experimental-seqair")]
#[instrument(level = "info", skip_all)]
pub fn get_pileups(
    readers: &mut PileupReaders,
    region: &ChunkRegion,
    params: &PileupMappingParams,
) -> Result<(Rc<Segment>, impl Iterator<Item = Pileup>)> {
    // Build the rastair Segment (FASTA fetch only) and extract ref window.
    let segment = readers.segment(region, 2).wrap_err("Failed to fetch segment")?;
    debug!(len = segment.sequence.len(), "Processing region (seqair)");

    let ref_window = {
        let bases: Arc<[Base]> =
            segment.sequence.iter().map(|&b| Base::from(b)).collect::<Vec<_>>().into();
        let start = Pos0::try_from(region.start).wrap_err("region start out of Pos0 range")?;
        ReferenceWindow { bases, start }
    };

    // Install the reference window so compute() can do mismatch inference.
    readers.inner_mut().customize_mut().reference = Some(ref_window);
    readers.inner_mut().customize_mut().read_flags = params.read_flags.clone();
    readers.inner_mut().customize_mut().unpaired = params.unpaired;
    readers.inner_mut().customize_mut().tag_requirement = params.require_tags.clone();
    readers.inner_mut().customize_mut().repeat_limit = params.indel_repeat_limit;
    readers.inner_mut().customize_mut().guess_orientation = params.guess_read_orientation;

    // Build the seqair Segment covering [region.start .. region.end + overfetch].
    // seqair's segment `end` is an *inclusive* 0-based position and must be
    // `<= contig_last_pos` (= contig length − 1). `region.last_position` is the
    // contig length, so clamp to `last_position - 1`; htslib tolerates ends one
    // past the contig, seqair rejects them.
    let overfetch: u64 = 2;
    let contig_last_pos = region.last_position.saturating_sub(1);
    let last = region.end.saturating_add(overfetch).min(contig_last_pos);
    let start = Pos0::try_from(region.start).wrap_err("region start out of Pos0 range")?;
    let end = Pos0::try_from(last).wrap_err("region end out of Pos0 range")?;
    let len_u32 =
        u32::try_from(last.saturating_sub(region.start).saturating_add(1)).unwrap_or(u32::MAX);
    const ONE: NonZeroU32 = NonZeroU32::MIN;
    let max_len = NonZeroU32::new(len_u32.max(1)).unwrap_or(ONE);
    let opts = SegmentOptions::new(max_len);
    let mut segs = readers
        .inner_mut()
        .segments((region.contig.as_str(), start, end), opts)
        .wrap_err("Failed to plan seqair segments")?;
    let seqair_seg = segs.next().wrap_err("No seqair segment for region")?;

    // Fetch BAM records + FASTA into PileupEngine (compute() runs here).
    let mut guard =
        readers.inner_mut().pileup(&seqair_seg).wrap_err("Failed to start seqair pileup")?;
    guard.set_max_depth(params.max_coverage);

    let segment = Rc::new(segment);
    let mut collector = NameCollector::new(params);
    let mut pileups: Vec<Pileup> = Vec::new();

    while let Some(col) = guard.pileups() {
        let pos = col.pos().as_u64();
        if !region.contains(pos) {
            continue;
        }
        match Pileup::from_seqair(&col, segment.clone(), params, &mut collector) {
            Ok(p) => pileups.push(p),
            Err(error) => {
                warn!(error = format!("{error:#}"), pos, "Failed to get pileup, skipping");
            }
        }
    }

    Ok((segment, pileups.into_iter()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::ReaderParams;
    use color_eyre::eyre::Result;

    #[test]
    fn test_reading_bounds() -> Result<()> {
        // check that we can read exactly the right positions
        let params = ReaderParams {
            regions: Some("chr19:6105700-6105800".parse().unwrap()),
            ..ReaderParams::test_data()
        };
        let mut readers = params.pileup_readers()?;
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
