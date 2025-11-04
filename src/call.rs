use crate::{
    bed::{
        rastair1::{BedParams, BedRecordsConvertParams, Rastair1BedFormat},
        writer::BedWriter,
    },
    call::{
        methylation::params::MethylationCallingParams, ml::MachineLearning,
        variant_calling::VariantCallingParams,
    },
    io::vcf_writer,
    metrics2::{MetricsForAlt, PileupMetrics},
    sequence::{ChunkRegion, ReaderParams, Readers},
    utils::{Surrounding, cli, logging::ThisIsABug as _, surrounding_pileups},
    vcf::{self, lowDp},
};
use clio::ClioPath;
use color_eyre::{
    Section,
    eyre::{ContextCompat as _, Result, WrapErr, ensure, eyre},
};
use ordered_channel::Receiver;
use rastair_vcf::VcfFilter;
use rayon::prelude::*;
use smol_str::SmolStr;
use std::{
    ops::Mul as _,
    thread::{self, available_parallelism},
};
use tracing::{debug, info, instrument, trace, warn};

pub mod denovo_cpg;
pub mod methylation;
pub mod metrics;
pub mod ml;
pub mod process;
mod record_filters;
pub mod variant_calling;
pub mod variants;

#[cfg(test)]
pub mod test_helpers;

#[derive(Debug, clap::Args, serde::Serialize)]
pub struct CallParams {
    // --- Input parameters ---
    #[command(flatten)]
    #[serde(skip)]
    pub segments: ReaderParams,
    #[command(flatten)]
    #[serde(skip)]
    pub segmentation: SegmentationParams,

    // --- Calling parameters ---
    #[command(flatten)]
    pub variant_calling: VariantCallingParams,
    #[command(flatten)]
    pub denovo_cpg: denovo_cpg::DenovoParams,
    #[command(flatten)]
    pub methylation: MethylationCallingParams,
    #[command(flatten)]
    pub ml: ml::MachineLearningParams,

    // --- Output parameters ---
    #[command(flatten)]
    pub record_filters: record_filters::RecordFilters,

    #[command(flatten)]
    #[serde(skip)]
    pub vcf: vcf_writer::VcfParams,

    #[command(flatten)]
    #[serde(skip)]
    pub bed: BedParams,

    // --- Other runtime parameters ---
    /// Number of threads to use for processing the BAM file. Will use all
    /// available threads when not specified.
    ///
    /// Note that VCF writing might use additional threads internally for compression.
    /// This can be overwritten with `--vcf-threads`.
    #[arg(short='@', long = "threads", default_value_t = available_parallelism().map(|n|n.get()).unwrap_or(2).max(1))]
    #[arg(help_heading = cli::sections::PROCESSING)]
    #[serde(skip)]
    pub total_threads: usize,
}

impl CallParams {
    fn figure_out_outputs(&mut self) -> Result<()> {
        let user_chose_output = self.vcf.vcf.is_some() || self.bed.bed.is_some();

        if user_chose_output {
            ensure!(
                self.vcf.vcf.as_ref() != self.bed.bed.as_ref(),
                "Can't write both VCF and BED output to the same file. Please specify different output files."
            );

            // If the user called rastair with something like `-o test.bed` (or
            // `-o test.bed.gz`), this is technically wrong: `-o` is short for
            // `--vcf` not for `--bed`.
            //
            // But we're gonna be nice about it and not error out but set the
            // `bed` field with that value instead (if no other `--bed` value is
            // given).
            if self.bed.bed.is_none()
                && let Some(vcf_filename) = self.vcf.vcf.as_ref()
                && let Some(filename) = vcf_filename.file_name()
                && let Some(filename) = filename.to_str()
                && (filename.ends_with(".bed") || filename.ends_with(".bed.gz"))
            {
                warn!(file=%vcf_filename, "VCF output file name ends with `.bed`/`.bed.gz`, did you mean to use `--bed` instead of `-o`/`--vcf`? Assuming you meant `--bed` and switching the output accordingly.");
                debug!(bed=?self.bed.bed, vcf=?self.vcf.vcf, "Switching output from VCF to BED");
                self.bed.bed = self.vcf.vcf.take();
            }
        } else if self.record_filters.cpgs_only {
            // Default to BED output if only CpGs are requested
            self.bed.bed = Some(ClioPath::std());
        } else {
            // Default to VCF output if no output is specified
            self.vcf.vcf = Some(ClioPath::std());
        }

        Ok(())
    }
}

#[derive(Debug, clap::Args, Clone)]
pub struct SegmentationParams {
    /// Maximum length of a segment in bases
    ///
    /// Used for splitting work between threads. Tweak this to adjust memory
    /// usage.
    #[arg(long, default_value_t = 100_000)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub segment_max_length: u64,

    /// Number of bases to overlap between segments
    ///
    /// Helpful to avoid missing variants at the edges of segments.
    #[arg(long, default_value_t = 200)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub segment_overlap: u64,
}

/// Read BAM + FASTA and call variants and methylation events
#[instrument(level = "debug", skip(params))]
pub fn call(mut params: CallParams) -> Result<()> {
    params.figure_out_outputs().wrap_err("Unclear output choice")?;
    let params = &params; // make params immutable for threads

    // Initialize readers for BAM and FASTA files
    let readers = params.segments.readers().wrap_err("Failed to read BAM/FASTA files")?;

    // Get segments that are small enough to process in RAM
    let regions: Vec<ChunkRegion> = readers
        .segments(params.segmentation.segment_max_length, params.segmentation.segment_overlap)
        .wrap_err("Could not fetch segments from BAM file")?
        .collect();
    if regions.is_empty() {
        warn!("No segments found in BAM file, nothing to do");
        return Ok(());
    } else if regions[0].region.len() < 2 {
        warn!(region=%regions[0].region, "Given range is one base long, this will not yield any results for context-specific methylation calling.");
    }

    // Init ML model if requested
    let ml = params.ml.init().wrap_err("Failed to initialize machine learning model")?;

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
    let worker_threads = params.total_threads.saturating_sub(writer_threads.get()).max(1);
    debug!(
        "Gonna use {} threads: {} for processing, {} for writing VCF",
        params.total_threads,
        worker_threads,
        writer_threads.get()
    );

    // The connection between the processing threads and the VCF writer is this
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
    let writer_thread =
        build_writer(params, &regions, vcf_receiver).wrap_err("VCF writer error")?;

    // Run this in a custom rayon thread pool to control the number of threads
    // and be able to tweak parameters when profiling
    rayon::ThreadPoolBuilder::new()
        .thread_name(|idx| format!("worker-{idx}"))
        .num_threads(worker_threads)
        .start_handler(|idx| trace!(idx, "Starting worker thread"))
        .exit_handler(|idx| trace!(idx, "Closing worker thread"))
        .build()
        .wrap_err("Failed to create thread pool for rayon")?
        .install(move || {
            regions_iter.par_bridge().try_for_each_with(
                (vcf_sender, params),
                |(vcf_sender, params), (index, region)| {
                    process_region_wrapper(index, region, vcf_sender, params, &ml)
                },
            )
        })
        .wrap_err("Failed to process regions in parallel")
        .note("Rastair might have still written an (incomplete) output file")?;

    writer_thread
        .join()
        .map_err(|_e| eyre!("Writer thread crashed"))
        .this_is_a_bug()? // this error is a panic in the thread
        .wrap_err("Error in writer thread")?; // this error is from actual result returned by the thread

    Ok(())
}

/// Build the VCF writer thread
fn build_writer(
    params: &CallParams,
    regions: &[ChunkRegion],
    vcf_receiver: Receiver<Vec<PileupMetrics>>,
) -> Result<thread::JoinHandle<Result<()>>> {
    let vcf_output = params.vcf.vcf.clone();
    let vcf_filter = params.record_filters.clone();
    let metadata = [
        format!("rastairVersion={}", env!("CARGO_PKG_VERSION")),
        format!("rastairCommand={}", std::env::args().skip(1).collect::<Vec<_>>().join(" ")),
        format!(
            "rastairConfig={}",
            serde_json::to_string(params)
                .wrap_err("Failed to serialize config to JSON")
                .this_is_a_bug()?
        ),
        format!("reference={}", params.segments.fasta_file),
    ];
    let mut vcf_writer =
        params.vcf.writer(regions, &metadata).wrap_err("Failed to create VCF writer")?;

    let bed = params.bed.clone();
    let mut bed_writer = bed.writer().wrap_err("Failed to create BED writer")?;
    let bed_params =
        BedRecordsConvertParams { ml_threshold: params.ml.ml, filters: bed.filters.clone() };

    // Spawn the actual VCF writer thread. Everything in here is driven by the
    // incoming records from the processing threads.
    //
    // The result returned from this thread is evaluated when the handle is joined.
    thread::Builder::new()
        .name("writer".to_string())
        .spawn(move || -> Result<()> {
            /// The segments we get have some overlap between them, so we need
            /// to ensure that we don't write the same record multiple times.
            #[derive(Default)]
            struct LastSeen {
                contig: Option<SmolStr>,
                pos: Option<u32>,
            }

            impl LastSeen {
                /// If this is new, returns true and updates the last seen record
                fn is_new(&mut self, contig: SmolStr, pos: u32) -> bool {
                    if self.contig.as_ref() == Some(&contig) && self.pos >= Some(pos) {
                        false
                    } else {
                        self.contig = Some(contig);
                        self.pos = Some(pos);
                        true
                    }
                }
            }

            let mut last_seen = LastSeen::default();

            let mut write = |record: &PileupMetrics| -> Result<()> {
                let vcf_record =
                    record.to_vcf_record().wrap_err("Failed to convert metrics to VCF record")?;

                if let Some(vcf_writer) = vcf_writer.as_mut()
                    && vcf_filter.matches(&vcf_record)
                {
                    vcf_writer.add(&vcf_record).wrap_err("Failed to write VCF record")?;
                }

                if let Some(bed_writer) = bed_writer.as_mut()
                    && (*vcf_record.info.in_cp_g || *vcf_record.info.de_novo_cp_g_candidate)
                    && let Some(bed_record) =
                        Rastair1BedFormat::from_record(&vcf_record, &bed_params)
                            .wrap_err("Failed to convert VCF record to BED format")
                            .this_is_a_bug()?
                {
                    bed_writer
                        .write_record(&bed_record)
                        .wrap_err("Failed to write record to BED")?;
                }

                Ok(())
            };

            // Since we only have the region index to ensure order, each
            // processing thread will send a vector of VCF records when it's
            // done with a region.
            for records in vcf_receiver {
                'current_batch: for record in &records {
                    if !last_seen.is_new(record.contig(), record.pos()) {
                        continue 'current_batch;
                    }
                    if let Err(e) = write(record) {
                        warn!(error = format!("{e:#}"), "Failed to write record, skipping");
                    }
                }
            }

            if let Some(vcf_output) = vcf_output.as_ref() {
                drop(vcf_writer);
                info!(file = %vcf_output, "Wrote VCF output");
            }
            if let Some(bed_output) = bed.bed.as_ref()
                && let Some(bed_writer) = bed_writer
            {
                bed_writer.close().wrap_err("Failed to close BED writer")?;
                info!(file = %bed_output, "Wrote BED output");
            }
            Ok(())
        })
        .wrap_err("Failed to spawn VCF writer thread")
}

/// Wrapper function for processing a region in a thread-safe manner.
///
/// Calls [`process_region`] with thread-local readers and ships the result to
/// the VCF writer.
#[instrument(level = "debug", skip_all, fields(region=%region.region))]
fn process_region_wrapper(
    index: usize,
    region: &ChunkRegion,
    vcf_sender: &mut ordered_channel::Sender<Vec<PileupMetrics>>,
    params: &CallParams,
    ml: &MachineLearning,
) -> Result<()> {
    thread_local! {
        /// Readers for the BAM and FASTA files, initialized per thread to avoid
        /// re-opening files or having a lock
        static READERS: std::cell::RefCell<Option<Readers>> = const { std::cell::RefCell::new(None) };
    }

    // Use thread-local readers to avoid re-opening files in each thread
    let res = READERS.with(|local_readers| -> Result<Vec<PileupMetrics>> {
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

        process_region(readers, region, params, ml)
    });

    let records = match res {
        Ok(records) => records,
        Err(e) => {
            warn!(error = format!("{e:#}"), "Failed to process region");
            // We still send an empty vector to the channel to increment the index
            Vec::new()
        }
    };

    if let Err(err) =
        vcf_sender.send(index, records).wrap_err("Failed to send records to VCF writer")
    {
        trace!(
            error = format!("{err:#}"),
            "Failed to send records to VCF writer, probably because the channel is closed"
        );
    }

    Ok(())
}

/// Analyse pileups in a region
fn process_region(
    readers: &mut Readers,
    region: &ChunkRegion,
    params: &CallParams,
    ml: &MachineLearning,
) -> Result<Vec<PileupMetrics>> {
    let pileup_mapping_params = process::PileupMappingParams {
        include_cpgs: params.methylation.should_include_all_cpgs(),
        variant_calling: params.variant_calling.clone(),
    };

    let pileups = region.process(readers, &pileup_mapping_params)?;
    let mut pileups: Vec<PileupMetrics> = pileups
        .into_iter()
        .map(PileupMetrics::try_from)
        .filter_map(|x: Result<PileupMetrics>| match x {
            Err(e) => {
                warn!(error = format!("{e:#}"), "failed to calculate metric, skipping");
                None
            }
            Ok(x) => Some(x),
        })
        .collect();

    // let mut records = piles
    //     .iter()
    //     .map(|pile| pile.variant_metrics(&params.variant_calling))
    //     .collect::<Result<Vec<_>>>()
    //     .wrap_err("Failed to collect metrics")?;

    // Call methylation events if requested
    // let record_len = metrics.len();
    // for i in 0..record_len {
    //     let (before, current, after) = surrounding_pileups(&mut metrics, i);

    //     if *current.info.read_depth < params.variant_calling.v_min_depth {
    //         current.filters.add_all(vcf::lowDp);
    //     }
    //     params.denovo_cpg.filter(current).wrap_err("Failed to add filters for de-novo CpGs")?;
    //     params.denovo_cpg.add_if_adjecent(current, before, after);

    //     params
    //         .methylation
    //         .call(current, before, after) // Might also add filters
    //         .wrap_err("Failed to call methylation")?;
    // }

    let pileups_len = pileups.len();
    for i in 0..pileups_len {
        let surrounding = surrounding_pileups(&mut pileups, i);

        // params.denovo_cpg.add_if_adjecent(&mut surrounding);

        let Surrounding { current, .. } = surrounding;
        // params.denovo_cpg.filter(current).wrap_err("Failed to add filters for de-novo CpGs")?;

        if current.pos_metrics.read_depth < params.variant_calling.v_min_depth {
            current.pos_filters.push(lowDp.filter());
        }
    }

    // Filter out piles that are not CpG if requested. We're doing this here to
    // not waste processing time on records we will discard anyway.
    if params.record_filters.cpgs_only {
        pileups.retain(|p| p.pileup.is_cpg || *p.pos_metrics.de_novo_cpg_candidate);
    }

    if !ml.disabled {
        /// Filter out very unlikely alts before running slow ML
        fn pre_ml_filter(c: &MetricsForAlt) -> bool {
            c.metrics.pos_metrics.read_depth > 1 && *c.metrics.pos_metrics.mapq > 5.
        }

        let pileups_len = pileups.len();
        for i in 0..pileups_len {
            let Surrounding { before, current, after } = surrounding_pileups(&mut pileups, i);

            for alt_base in current.alts() {
                let alt = current
                    .alt_metrics(alt_base)
                    .wrap_err("Failed to get alt metrics")
                    .this_is_a_bug()?;
                if pre_ml_filter(&alt)
                    && let Some(pred) = ml.predict2(&alt, before, after)
                    && pred.pass()
                {
                    let filters = current
                        .alt_filters_mut(alt_base)
                        .wrap_err("Failed to get mutable alt metrics")
                        .this_is_a_bug()?;
                    filters.ml.replace(pred.prediction);
                }
            }
        }
    }

    // Okay, here is what we have:
    // - `metrics`: The pileup metrics for each position
    //   - `pos_metrics`: The position-level metrics
    //   - `ref_metrics`: The reference allele metrics
    //   - `alt_metrics`: The alt allele metrics
    //      - `ml`: The ML prediction that this is a true variant
    // Now, we need to convert these into VCF records.

    Ok(vec![])
}
