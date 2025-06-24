use crate::{
    call::{methylation::params::MethylationCallingParams, variant_calling::VariantCallingParams},
    sequence::{ChunkRegion, SegmentsParams},
    vcf::{self},
    vcf_writer,
};
use color_eyre::eyre::{Result, WrapErr};
use rastair2_vcf::Vcf;
use tracing::{debug, info, instrument, warn};

pub mod methylation;
pub mod metrics;
pub mod process;
pub mod variant_calling;
pub mod variants;

#[cfg(test)]
pub mod test_helpers;

#[derive(Debug, clap::Args)]
pub struct CallParams {
    #[command(flatten)]
    segments: SegmentsParams,

    #[command(flatten)]
    variant_calling: VariantCallingParams,

    #[command(flatten)]
    methylation: MethylationCallingParams,

    #[command(flatten)]
    vcf: vcf_writer::Params,

    /// Number of threads to use for processing the BAM file, 0 means auto-detect
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
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
    } else if regions[0].region.len() < 2 {
        warn!(region=%regions[0].region, "Given range is one base long, this will not yield any results for context-specific methylation calling.");
    }

    debug!("Going to process {} segments", regions.len());

    // Create a VCF writer for the output
    let mut vcf_writer = params.vcf.vcf_writer(&regions).wrap_err("failed to create VCF writer")?;

    let pileup_mapping_params = process::PileupMappingParams {
        include_cpgs: params.methylation.should_include_all_cpgs(),
        keep_overlapping_reads: params.variant_calling.keep_overlapping_reads,
    };

    // Process each region and write results to the VCF
    // TODO: For multithreaded processing, have readers per thread, collect data in order, and write in main thread
    regions.into_iter().try_for_each(|region| {
        regions_seen += 1;
        region
            .process(&mut readers, &pileup_mapping_params)
            .and_then(|piles| {
                piles.into_iter().try_for_each(|pile| {
                    pile.variant_metrics(&params.variant_calling)
                        .wrap_err("Failed to collect metrics")
                        .and_then(|record| {
                            params.methylation.call(record).wrap_err("Failed to call methylation")
                        })
                        .and_then(|record| {
                            write_pileup(&record, &mut vcf_writer)
                                .wrap_err("failed to write to VCF")
                        })
                        .wrap_err_with(|| {
                            format!("Failed to process pileup {}:{}", pile.chrom(), pile.pos)
                        })
                })
            })
            .wrap_err_with(|| format!("failed to process region {}", region.region))
    })?;

    drop(vcf_writer); // Ensure all data is flushed to the output file
    info!(file = %params.vcf.vcf_output, "Wrote output");

    Ok(())
}

/// Write a pileup to the VCF output
#[instrument(level = "trace", skip_all)]
fn write_pileup(record: &vcf::Record, output: &mut Vcf<vcf::Record>) -> Result<()> {
    output.add(record).wrap_err("Failed to add record")
}
