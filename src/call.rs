use crate::{
    call::{variants::VariantCandidatePileup, vcf::Filters},
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
use tracing::{Level, debug, info, instrument, trace, warn};

use filtering::threshold::ThresholdConfig;
use variants::{PositionInRead, SeenBase, SeenBases};

mod variant_calling;
pub mod variants;
mod filtering {
    pub mod threshold;
}
mod metrics;
pub mod vcf;
pub mod vcf_writer;

#[derive(Debug, clap::Args)]
pub struct CallParams {
    #[command(flatten)]
    segments: SegmentsParams,

    #[command(flatten)]
    thresholds: filtering::threshold::ThresholdConfig,

    #[command(flatten)]
    vcf: vcf_writer::Params,
}

/// Read BAM + FASTA and call variants and methylation events
#[instrument(level = "debug", skip(params))]
pub fn call(params: &CallParams) -> Result<()> {
    // Initialize readers for BAM and FASTA files
    let mut readers = params.segments.readers().wrap_err("failed to fetch segments")?;

    // Get segments that are small enough to process in RAM
    let mut regions_seen = 0;
    let regions: Vec<ChunkRegion> =
        readers.segments().wrap_err("Could not fetch segments from BAM file")?.collect();
    if regions.is_empty() {
        warn!("No segments found in BAM file, nothing to do");
        return Ok(());
    }
    debug!("Going to process {} segments", regions.len());

    // Create a VCF writer for the output
    let mut vcf_writer = params.vcf.vcf_writer(&regions).wrap_err("failed to create VCF writer")?;

    // Process each region and write results to the VCF
    // TODO: For multithreaded processing, collect data in order and write in main thread
    regions.into_iter().try_for_each(|region| {
        regions_seen += 1;
        process_region(&region, &mut readers, &params.thresholds)
            .and_then(|piles| {
                piles.into_iter().try_for_each(|pile| {
                    write_pileup(&pile, &mut vcf_writer).wrap_err_with(|| {
                        format!("Failed to write pileup for {}:{}", pile.chrom, pile.pos)
                    })
                })
            })
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
    _thresholds: &ThresholdConfig,
) -> Result<Vec<VariantCandidatePileup>> {
    let segment = readers.segment(region).wrap_err("failed to fetch segment")?;
    trace!(len = segment.sequence.len(), "Processing region");

    // Fetch the pileups for the segment
    FetchDefinition::try_from(&segment)
        .wrap_err("Could not convert region string")
        .and_then(|r| readers.bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))
        .wrap_err_with(|| format!("Could not fetch region `{}` from BAM file", region.region))?;

    // Go over each column in the pileup and collect variant candidates
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
        .filter_map(|res| match res {
            Ok(x) => Some(x),
            Err(error) => {
                warn!(%error, "Failed to get pileup, skipping");
                None
            }
        })
        .collect::<Vec<_>>();

    if tracing::enabled!(Level::DEBUG) {
        if piles.is_empty() {
            trace!("No candidate piles found in region, skipping");
            return Ok(piles);
        } else {
            let count = readable::num::Unsigned::from(piles.len());
            let bytes = readable::byte::Byte::from(
                piles.len() * std::mem::size_of::<VariantCandidatePileup>(),
            );
            debug!(%count, %bytes, "Collected candidates");
        }
    }

    // TODO: Add at least threshold filtering for methylation calling :)

    Ok(piles)
}

/// Is this pileup a candidate for a variant?
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

    if bases.matches(reference_base) {
        // Matches reference base, boring.
        // trace!(?bases, pos = pile.pos(), ?reference_base, ?next_base, "pile matches reference");
        Ok(None)
    } else {
        Ok(Some(VariantCandidatePileup {
            chrom: region.chromosome.clone(),
            pos: pile.pos(),
            bases,
            reference_base,
            sequence_before: segment
                .sequence_slice::<2>(idx.saturating_sub(2), idx)
                .wrap_err("failed to get sequence before pileup")?,
            sequence_after: segment
                .sequence_slice::<2>(idx + 1, idx + 3)
                .wrap_err("failed to get sequence after pileup")?,
        }))
    }
}

/// Collect info from a pileup alignment
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

/// Write a pileup to the VCF output
#[instrument(level = "trace", skip_all)]
fn write_pileup(pile: &VariantCandidatePileup, output: &mut Vcf<vcf::Record>) -> Result<()> {
    // TODO: Maybe pull out the metrics
    let metrics = pile.metrics().wrap_err("Failed to calculate metrics")?;
    let calling_metrics = pile.calling_metrics().wrap_err("Failed to calculate calling metrics")?;

    let record = vcf::Record {
        fixed_fields: pile.fixed_fields(),
        // FIXME: Add filters based on thresholds
        filters: Filters::new().add(rastair2_vcf::standard_fields::PASS),
        info: metrics,
        samples: smallvec::smallvec![calling_metrics],
    };
    output.add(record).wrap_err("Failed to add record")
}
