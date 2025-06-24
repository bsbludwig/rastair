use crate::{
    call::{methylation::params::MethylationCallingParams, variant_calling::VariantCallingParams},
    sequence::{ChunkRegion, Readers, SegmentsParams},
    vcf, vcf_writer,
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr, eyre};
use rayon::prelude::*;
use std::{ops::Mul as _, thread};
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
    ///
    /// Note that VCF writing might use additional threads internally for compression.
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
}

/// Read BAM + FASTA and call variants and methylation events
#[instrument(level = "debug", skip(params))]
pub fn call(params: &CallParams) -> Result<()> {
    // Initialize readers for BAM and FASTA files
    let readers = params.segments.readers().wrap_err("failed to fetch segments")?;

    // Get segments that are small enough to process in RAM
    let regions: Vec<ChunkRegion> =
        readers.segments().wrap_err("Could not fetch segments from BAM file")?.collect();
    if regions.is_empty() {
        warn!("No segments found in BAM file, nothing to do");
        return Ok(());
    } else if regions[0].region.len() < 2 {
        warn!(region=%regions[0].region, "Given range is one base long, this will not yield any results for context-specific methylation calling.");
    }

    debug!("Going to process {} segments", regions.len());

    // Process each region and write results to the VCF
    //
    // There are two aspects to this: Collecting the variant candidates and
    // writing them to the VCF. To use all CPU available, we use rayon to
    // process the regions in parallel. From there, we send ready-made VCF
    // records to a special writer thread that only deals with writing the VCF
    // file. (In reality, writing also uses multiple threads internally for
    // parallel compression, but we don't have to care about that here.)

    // The connection between the processing threads and the VCF writer this
    // ordered channel. It buffers `Vec<vcf::Record>`s, alongside the index from
    // the parallel iterator.

    let (vcf_sender, vcf_receiver) = {
        // At least 5x buffer for VCF records to account for reordering and processing time
        let buffer_size = params.threads.max(2).mul(5);
        ordered_channel::bounded(buffer_size)
    };

    // Create a VCF writer for the output
    let writer_thread = thread::Builder::new()
        .name("vcf_writer".to_string())
        .spawn({
            let vcf_output = params.vcf.vcf_output.clone();
            let mut vcf_writer =
                params.vcf.vcf_writer(&regions).wrap_err("failed to create VCF writer")?;
            move || -> Result<()> {
                for records in vcf_receiver {
                    for record in records {
                        vcf_writer.add(&record).wrap_err("failed to write record to VCF")?;
                    }
                }

                drop(vcf_writer); // Ensure all data is flushed to the output file
                info!(file = %vcf_output, "Wrote output");

                Ok(())
            }
        })
        .wrap_err("failed to spawn VCF writer thread")?;

    // Run this in a custom rayon thread pool to control the number of threads
    // and be able to tweak parameters when profiling
    rayon::ThreadPoolBuilder::new()
        .thread_name(|idx| format!("worker-{idx}"))
        .num_threads(params.threads)
        .build()
        .wrap_err("failed to create thread pool for rayon")?
        .install(move || {
            regions.into_iter().enumerate().par_bridge().try_for_each_with(
                (vcf_sender, params),
                |(vcf_sender, params), (index, region)| {
                    process_region(index, region, vcf_sender, params)
                },
            )
        })?;

    writer_thread
        .join()
        .map_err(|e| eyre!("{e:?}"))
        .wrap_err("failed to join VCF writer thread")?
        .wrap_err("writer thread error")?;

    Ok(())
}

thread_local! {
    static READERS: std::cell::RefCell<Option<Readers>> = const { std::cell::RefCell::new(None) };
}

#[instrument(level = "debug", skip_all, fields(region=%region.region))]
fn process_region(
    index: usize,
    region: ChunkRegion,
    vcf_sender: &mut ordered_channel::Sender<Vec<vcf::Record>>,
    params: &CallParams,
) -> Result<()> {
    let res = READERS.with(|tl| -> Result<Vec<vcf::Record>> {
        let mut tl_borrow = tl.borrow_mut();

        // Initialize thread-local resources if not already done
        if tl_borrow.is_none() {
            *tl_borrow = Some(params.segments.readers().wrap_err("fetch readers")?);
        }
        let readers = tl_borrow.as_mut().wrap_err("failed to access thread-local resources")?;

        let pileup_mapping_params = process::PileupMappingParams {
            include_cpgs: params.methylation.should_include_all_cpgs(),
            keep_overlapping_reads: params.variant_calling.keep_overlapping_reads,
        };

        let piles =
            region.process(readers, &pileup_mapping_params).wrap_err("Failed to process region")?;

        let records = piles
            .into_iter()
            .map(|pile| {
                pile.variant_metrics(&params.variant_calling)
                    .wrap_err("Failed to collect metrics")
                    .and_then(|record| {
                        params.methylation.call(record).wrap_err("Failed to call methylation")
                    })
                    .wrap_err_with(|| {
                        format!("Failed to process pileup {}:{}", pile.chrom(), pile.pos)
                    })
            })
            .collect::<Result<Vec<_>>>()
            .wrap_err("Failed to collect metrics")?;

        Ok(records)
    });

    let records = match res {
        Ok(records) => records,
        Err(e) => {
            warn!(e = format!("{e:#}"), "Failed to process region");
            Vec::new() // we still want to send an empty vector to the VCF writer to increment the index
        }
    };

    vcf_sender.send(index, records).wrap_err("Failed to send records to VCF writer")?;

    Ok(())
}
