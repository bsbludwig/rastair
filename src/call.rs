use crate::{
    sequence::{ChunkRegion, Readers, Segment, SegmentsParams},
    utils::TryAsBase as _,
};
use color_eyre::eyre::{Context, ContextCompat, Result};
use rastair2_vcf::Vcf;
use rust_htslib::bam::{
    FetchDefinition, Read as _,
    pileup::{Alignment, Pileup},
};
use smallvec::SmallVec;
use tracing::{debug, info, instrument, trace, warn};

use filtering::threshold::ThresholdConfig;
use methylation_event_writer::MethylationEventWriter;
use variants::{PositionInRead, SeenBase, SeenBases};

mod methylation;
mod methylation_event_writer;
pub mod scores;
pub mod variants;
mod filtering {
    pub mod threshold;
}
pub mod vcf;
pub mod vcf_writer;

// Re-exports for VCF writer
pub use scores::VariantCandidatePileupMetrics;
pub use variants::VariantCandidatePileup;

#[derive(Debug, clap::Args)]
pub struct CallParams {
    #[command(flatten)]
    segments: SegmentsParams,

    #[command(flatten)]
    thresholds: filtering::threshold::ThresholdConfig,

    #[command(flatten)]
    vcf: vcf_writer::Params,
}

#[instrument(level = "debug", skip(params))]
pub fn call(params: &CallParams) -> Result<()> {
    let mut readers = params.segments.readers().wrap_err("failed to fetch segments")?;

    let mut regions_seen = 0;
    let regions: Vec<ChunkRegion> =
        readers.segments().wrap_err("Could not fetch segments from BAM file")?.collect();
    if regions.is_empty() {
        warn!("No segments found in BAM file, nothing to do");
        return Ok(());
    }
    debug!("Going to process {} segments", regions.len());

    let mut vcf_writer = params.vcf.vcf_writer(&regions).wrap_err("failed to create VCF writer")?;

    regions.into_iter().try_for_each(|region| {
        regions_seen += 1;
        process_region(&region, &mut readers, &params.thresholds, &mut vcf_writer)
            .wrap_err_with(|| format!("failed to process region {}", region.region))
    })?;

    info!("Processed {regions_seen} segments");
    info!("Wrote output to {}", params.vcf.vcf_output.display());

    return Ok(());
}

#[instrument(level = "info", skip_all, fields(region=%region.region))]
fn process_region(
    region: &ChunkRegion,
    readers: &mut Readers,
    thresholds: &ThresholdConfig,
    output: &mut Vcf<vcf::Record>,
) -> Result<()> {
    let segment = readers.segment(region).wrap_err("failed to fetch segment")?;
    trace!(len = segment.sequence.len(), "Processing region");

    FetchDefinition::try_from(&segment)
        .wrap_err("Could not convert region string")
        .and_then(|r| readers.bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))
        .wrap_err_with(|| format!("Could not fetch region `{}` from BAM file", region.region))?;
    let piles = readers
        .bam
        .pileup()
        .filter_map(|p| p.ok())
        .filter(|p| {
            // Filter out pileups that are not in the region of interest
            region.contains(u64::from(p.pos()))
        })
        .flat_map(|pile| {
            collect_candidate(&pile, &segment, &segment.range)
                .wrap_err_with(|| {
                    format!("Failed to get candidate from pileup at position {}", pile.pos())
                })
                .transpose()
        })
        .map(|pile| -> Result<(VariantCandidatePileup, VariantCandidatePileupMetrics)> {
            let pile = pile?;
            let metrics = pile.metrics().wrap_err("Failed to calculate metrics")?;
            Ok((pile, metrics))
        })
        .filter_map(|res| match res {
            Ok(x) => Some(x),
            Err(error) => {
                warn!(%error, "Failed to get pileup, skipping");
                None
            }
        })
        .collect::<Vec<_>>();

    if piles.is_empty() {
        trace!("No candidate piles found in region, skipping");
        return Ok(());
    } else {
        let count = readable::num::Unsigned::from(piles.len());
        let bytes = readable::byte::Byte::from(
            piles.len()
                * std::mem::size_of::<(VariantCandidatePileup, VariantCandidatePileupMetrics)>(),
        );
        debug!(%count, %bytes, "Collected candidates");
    }

    // At this point we have both the pileup and the metrics for each pileup.
    // We can now gather metrics for each of the reads that were used to create the pileup.

    // TODO: Emit variant candidates here if user asked for it

    // Now we can analyze possible methylation events
    let mut iterator = piles.into_iter().peekable();
    loop {
        let Some((pile, metrics)) = iterator.next() else {
            // last item in iter
            break;
        };
        // This returns a reference, so it's easist to do a `loop` here instead
        // of writing a fancy adaptor
        let _next = iterator.peek();

        // TODO: Use `next` here
        if pile.likely_methylation_event(&metrics, thresholds) {
            MethylationEventWriter(&pile, &metrics)
                .write(output)
                .wrap_err("Failed to write methylation event to VCF")?;
        }
    }

    Ok(())
}

#[instrument(level = "trace", skip_all)]
fn collect_candidate(
    pile: &Pileup,
    segment: &Segment,
    region: &ChunkRegion,
) -> Result<Option<VariantCandidatePileup>> {
    let segment_start_pos =
        usize::try_from(segment.range.start).expect("segment range fits in usize");
    let idx = usize::try_from(pile.pos())
        .expect("position fits in usize")
        .checked_sub(segment_start_pos)
        .wrap_err_with(|| {
            format!(
                "pile position {} is not in segment range {}..{}",
                pile.pos(),
                segment.range.start,
                segment.range.end
            )
        })?;
    let bases = SeenBases(pile.alignments().filter_map(pileup_mapper).collect());
    let reference_base =
        segment.sequence.get(idx).wrap_err("failed to get reference base")?.as_base()?;
    let next_base = segment.sequence.get(idx + 1).and_then(|x| x.as_base().ok());
    if bases.is_variant_candidate() {
        // info!(?bases, pos = pile.pos(), ?reference_base, ?next_base, "found pile of interest");
        Ok(Some(VariantCandidatePileup {
            chrom: region.chromosome.clone(),
            pos: pile.pos(),
            bases,
            reference_base,
            next_base,
        }))
        // info!(?pileup, metrics=?pileup.metrics(), "variant candidate");
    } else if bases.matches(reference_base) {
        // Matches reference base
        // boring.
        // trace!(?bases, pos = pile.pos(), ?reference_base, ?next_base, "pile matches reference");
        Ok(None)
    } else {
        warn!(
            ?bases,
            pos = pile.pos(),
            ?reference_base,
            ?next_base,
            "pile does not match reference but is also not interesting"
        );
        Ok(None)
    }
}

fn pileup_mapper(a: Alignment<'_>) -> Option<SeenBase> {
    let pos = a.qpos()?;
    let record = a.record();
    if !record.is_proper_pair() {
        // fixme: maybe be more lenient here
        return None;
    }
    if record.is_quality_check_failed() {
        return None;
    }
    // fixme: understand this better:
    // if record.cigar().iter().any(|c| matches!(c, Cigar::SoftClip(_))) {
    //     return None;
    // }

    Some(SeenBase {
        qname: SmallVec::from(record.qname()),
        base: record.seq()[pos].as_base().ok()?, // fixme: handle error or at least check usual error modes
        qual: *record.qual().get(pos)?,
        mapq: record.mapq(),
        reverse: record.is_reverse(),
        position: PositionInRead {
            pos: u32::try_from(pos).expect("position fits in u32"),
            read_length: u32::try_from(record.seq().len()).expect("read length fits in u32"),
        },
    })
}
