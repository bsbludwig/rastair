use crate::{
    bed::BedWriter,
    call::{methylation::params::MethylationCallingParams, variant_calling::VariantCallingParams},
    io::vcf_writer,
    sequence::{ChunkRegion, Readers, SegmentsParams},
    vcf,
};
use clio::ClioPath;
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr, ensure, eyre};
use rayon::prelude::*;
use smol_str::SmolStr;
use std::{
    ops::Mul as _,
    thread::{self, available_parallelism},
};
use tracing::{debug, info, instrument, warn};

pub mod denovo_cpg;
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
    denovo_cpg: denovo_cpg::DenovoParams,

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

    /// Output BED file with the called methylation events.
    #[arg(long = "bed")]
    pub bed_output: Option<ClioPath>,
}

/// Read BAM + FASTA and call variants and methylation events
#[instrument(level = "debug", skip(params))]
pub fn call(params: &CallParams) -> Result<()> {
    ensure!(
        Some(&params.vcf.vcf_output) != params.bed_output.as_ref(),
        "Can't write both VCF and BED output to the same file. Please specify different output files."
    );

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

            let bed_output = params.bed_output.clone();
            let mut bed_writer = bed_output
                .as_ref()
                .map(|bed_output| {
                    BedWriter::new(bed_output).wrap_err("Failed to create BED writer")
                })
                .transpose()?;
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
                        if last_seen_chrom.as_ref() == Some(&record.main.chrom)
                            && last_seen_pos >= Some(record.main.pos)
                        {
                            continue;
                        }
                        // Seen a new record, update the last seen
                        last_seen_chrom = Some(record.main.chrom.clone());
                        last_seen_pos = Some(record.main.pos);

                        writer.add(record).wrap_err("Failed to write record to VCF")?;

                        if let Some(bed_writer) = bed_writer.as_mut()
                            && (*record.info.in_cp_g || *record.info.de_novo_cp_g_candidate)
                        {
                            // Write the record to the BED file if requested
                            bed_writer
                                .write_record(
                                    &record
                                        .try_into()
                                        .wrap_err("Failed to convert record to BED format")?,
                                )
                                .wrap_err("Failed to write record to BED")?;
                        }
                    }
                }

                // Ensure all data is flushed to the output file
                drop(writer);

                info!(file = %vcf_output, "Wrote VCF output");
                if let Some(bed_output) = bed_output {
                    info!(file = %bed_output, "Wrote BED output");
                }
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
            read_masking: params.variant_calling.read_masking.clone(),
            read_flags: params.variant_calling.read_flags.clone(),
        };

        let piles = region.process(readers, &pileup_mapping_params)?;

        let mut records = piles
            .iter()
            .map(|pile| pile.variant_metrics(&params.variant_calling))
            .collect::<Result<Vec<_>>>()
            .wrap_err("Failed to collect metrics")?;

        // Call methylation events if requested
        let record_len = records.len();
        for i in 0..record_len {
            let (before, current, after) = surrounding_records(&mut records, i);

            params.denovo_cpg.filter(current).wrap_err("Failed to add filters for de-novo CpGs")?;

            params
                .methylation
                .call(current, before, after) // Might also add filters
                .wrap_err("Failed to call methylation")?;
        }

        if params.variant_calling.cpgs_only {
            // Filter out piles that are not CpG if requested
            records.retain(|record| *record.info.in_cp_g || *record.info.de_novo_cp_g_candidate);
        }

        Ok(records)
    });

    let records = match res {
        Ok(records) => records,
        Err(e) => {
            warn!(error = format!("{e:#}"), "Failed to process region");
            // We still send an empty vector to the channel to increment the index
            Vec::new()
        }
    };

    vcf_sender.send(index, records).wrap_err("Failed to send records to VCF writer")?;

    Ok(())
}

/// Get the surrounding records for a given index in the records slice.
fn surrounding_records(
    records: &mut [vcf::Record],
    index: usize,
) -> (Option<&vcf::Record>, &mut vcf::Record, Option<&vcf::Record>) {
    // To appease the borrow checker and get a mutable reference to the current record,
    // we split the records into three parts.
    let (left, right) = records.split_at_mut(index);
    let (current_slice, next_slice) = right.split_at_mut(1);
    let current = &mut current_slice[0];

    let before = left.last();
    let after = next_slice.first();
    // we might not have the direct neighbors
    let before = before.filter(|r| {
        r.main.chrom == current.main.chrom && Some(r.main.pos) == current.main.pos.checked_sub(1)
    });
    let after = after.filter(|r| {
        r.main.chrom == current.main.chrom && Some(r.main.pos) == current.main.pos.checked_add(1)
    });

    (before, current, after)
}
