use crate::{
    call::methylation::params::MethylationCallingParams,
    sequence::{ChunkRegion, SegmentsParams},
    vcf::{self},
    vcf_writer,
};
use color_eyre::eyre::{Result, WrapErr};
use rastair2_vcf::Vcf;
use tracing::{debug, info, instrument, warn};

mod methylation;
mod metrics;
mod process;
mod variant_calling;
pub mod variants;

#[cfg(test)]
pub mod test_helpers;

#[derive(Debug, clap::Args)]
pub struct CallParams {
    #[command(flatten)]
    segments: SegmentsParams,

    #[command(flatten)]
    methylation: MethylationCallingParams,

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
    // TODO: For multithreaded processing, have readers per thread, collect data in order, and write in main thread
    regions.into_iter().try_for_each(|region| {
        regions_seen += 1;
        region
            .process(&mut readers)
            .and_then(|piles| {
                piles.into_iter().try_for_each(|pile| {
                    pile.variant_metrics()
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
    info!("Wrote output to {}", params.vcf.vcf_output.display());

    return Ok(());
}

/// Write a pileup to the VCF output
#[instrument(level = "trace", skip_all)]
fn write_pileup(record: &vcf::Record, output: &mut Vcf<vcf::Record>) -> Result<()> {
    output.add(record).wrap_err("Failed to add record")
}
