use crate::{
    call::{methylation::params::MethylationCallingParams, variant_calling::VariantCallingParams},
    io::vcf_writer,
    sequence::{ChunkRegion, Readers, SegmentsParams},
    vcf,
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr, eyre};
use rayon::prelude::*;
use smol_str::SmolStr;
use std::{
    ops::Mul as _,
    thread::{self, available_parallelism},
};
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

    /// Number of threads to use for processing the BAM file. Will use all
    /// available threads when not specified.
    ///
    /// Note that VCF writing might use additional threads internally for compression.
    /// This can be overwritten with `--vcf-threads`.
    #[arg(short='@', long = "threads", default_value_t = available_parallelism().map(|n|n.get()).unwrap_or(2).max(1))]
    pub threads: usize,
}

/// Read BAM + FASTA and call variants and methylation events
#[instrument(level = "debug", skip(params))]
pub fn call(params: &CallParams) -> Result<()> {
    // Initialize readers for BAM and FASTA files
    let readers = params.segments.readers().wrap_err("Failed to fetch segments")?;

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
    // This is done in parallel to speed up the processing, so here are a few
    // comments on this works. There are two aspects to this: Collecting the
    // variant candidates and writing them to the VCF.
    //
    // To use all CPU available, we use rayon to process the regions in
    // parallel. From there, we send ready-made VCF records to a special writer
    // thread that only deals with writing the VCF file.
    let writer_threads = params.vcf.vcf_threads;
    let worker_threads = params.threads.saturating_sub(writer_threads.get()).max(1);

    // The connection between the processing threads and the VCF writer this
    // ordered channel. It buffers `Vec<vcf::Record>`s, alongside the index from
    // the parallel iterator.
    let (vcf_sender, vcf_receiver) = {
        // At least 10x buffer for VCF records to account for reordering and processing time
        let buffer_size = worker_threads.mul(10);
        ordered_channel::bounded(buffer_size)
    };

    // We're going to go over the regions in parallel, and add the index of the
    // region here for the ordered channel.
    let regions_iter = regions.iter().enumerate();

    // Create a VCF writer for the output
    let writer_thread = thread::Builder::new()
        .name("writer".to_string())
        .spawn({
            let vcf_output = params.vcf.vcf_output.clone();
            let metadata = [
                format!("rastair2Version={}", env!("CARGO_PKG_VERSION")),
                format!(
                    "rastair2Command={}",
                    std::env::args().skip(1).collect::<Vec<_>>().join(" ")
                ),
                format!("reference={}", params.segments.fasta_file),
            ];
            let mut writer =
                params.vcf.writer(&regions, &metadata).wrap_err("Failed to create VCF writer")?;
            move || -> Result<()> {
                // The segments we get have some overlap between them, so we
                // need to ensure that we don't write the same record multiple
                // times.
                let mut last_seen_chrom: Option<SmolStr> = None;
                let mut last_seen_pos: Option<u32> = None;

                // Since we only have the region index to ensure order, each
                // processing thread will send a vector of VCF records when it's
                // done with a region.
                for records in vcf_receiver {
                    for record in &records {
                        let record: &vcf::Record = record;

                        // Skip records that are already seen
                        if last_seen_chrom.as_ref() == Some(&record.fixed_fields.chrom)
                            && last_seen_pos >= Some(record.fixed_fields.pos)
                        {
                            continue;
                        }
                        // Seen a new record, update the last seen
                        last_seen_chrom = Some(record.fixed_fields.chrom.clone());
                        last_seen_pos = Some(record.fixed_fields.pos);

                        writer.add(record).wrap_err("Failed to write record to VCF")?;
                    }
                }

                // Ensure all data is flushed to the output file
                drop(writer);

                info!(file = %vcf_output, "Wrote output");
                Ok(())
            }
        })
        .wrap_err("Failed to spawn VCF writer thread")?;

    // Run this in a custom rayon thread pool to control the number of threads
    // and be able to tweak parameters when profiling
    rayon::ThreadPoolBuilder::new()
        .thread_name(|idx| format!("worker-{idx}"))
        .num_threads(worker_threads)
        .build()
        .wrap_err("Failed to create thread pool for rayon")?
        .install(move || {
            regions_iter.par_bridge().try_for_each_with(
                (vcf_sender, params),
                |(vcf_sender, params), (index, region)| {
                    process_region(index, region, vcf_sender, params)
                },
            )
        })?;

    writer_thread
        .join()
        .map_err(|e| eyre!("{e:?}"))
        .wrap_err("Failed to join VCF writer thread")?
        .wrap_err("writer thread error")?;

    Ok(())
}

thread_local! {
    /// Readers for the BAM and FASTA files, initialized per thread to avoid
    /// re-opening files or having a lock
    static READERS: std::cell::RefCell<Option<Readers>> = const { std::cell::RefCell::new(None) };
}

#[instrument(level = "debug", skip_all, fields(region=%region.region))]
fn process_region(
    index: usize,
    region: &ChunkRegion,
    vcf_sender: &mut ordered_channel::Sender<Vec<vcf::Record>>,
    params: &CallParams,
) -> Result<()> {
    // Use thread-local readers to avoid re-opening files in each thread
    let res = READERS.with(|local_readers| -> Result<Vec<vcf::Record>> {
        let mut local_readers = local_readers.borrow_mut();
        let readers = {
            // Initialize thread-local readers first time the thread accesses them
            if local_readers.is_none() {
                let readers = params
                    .segments
                    .readers()
                    .wrap_err("Failed to open readers in worker thread")?;
                *local_readers = Some(readers);
            }
            local_readers.as_mut().wrap_err("Failed to access thread-local resources")?
        };

        let pileup_mapping_params = process::PileupMappingParams {
            include_cpgs: params.methylation.should_include_all_cpgs(),
            keep_overlapping_reads: params.variant_calling.keep_overlapping_reads,
        };

        let piles =
            region.process(readers, &pileup_mapping_params).wrap_err("Failed to process region")?;

        let records = piles
            .into_iter()
            .filter(|pile| {
                // Filter out piles that are not CpG if requested
                !params.variant_calling.cpgs_only || pile.is_cpg
            })
            .map(|pile| {
                pile.variant_metrics(&params.variant_calling)
                    .wrap_err("Failed to collect metrics")
                    .and_then(|record| {
                        params.methylation.call(record).wrap_err("Failed to call methylation")
                    })
                    // chain the steps above to add outer error context
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
            // We still send an empty vector to the channel to increment the index
            Vec::new()
        }
    };

    vcf_sender.send(index, records).wrap_err("Failed to send records to VCF writer")?;

    Ok(())
}
