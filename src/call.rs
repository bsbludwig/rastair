//! This is the most complex part of Rastair
//!
//! Its data flow is as follows:
//!
//! - Input: Paths for BAM and FASTA, config parameters, output path and format
//! - Split the genome into segments and process them in parallel:
//!   1. Load reads from BAM overlapping the segment
//!   2. Build pileups for each position in the segment
//!   3. Calculate metrics for each pileup
//!   4. Pre-filter positions (e.g. only keep CpG sites)
//!   5. Call variants based on the metrics
//! - Output: Write variants to output file in specified format (VCF/BCF)
//!   1. One output thread collects records for each segment
//!   2. Filter be given criteria
//!   3. Convert to the output format
//!   4. Write to file in order

use crate::{
    bed::rastair1::BedParams,
    call::{
        methylation::params::MethylationCallingParams,
        pileup::{Pileup, SimpleRead},
        process::{calculate_pileup_metrics, get_pileups},
        variant_calling::VariantCallingParams,
    },
    io::vcf_writer,
    metrics::{self, MethylationEvidenceStrandInfo, PileupMetrics, ml::types::MachineLearning},
    sequence::{ChunkRegion, ReaderParams, Readers, Segment, SegmentationParams},
    utils::{PileupMetricsIteratorExt, cli, logging::ThisIsABug as _},
};
use clio::ClioPath;
use color_eyre::{
    Section,
    eyre::{ContextCompat as _, Result, WrapErr, ensure, eyre},
};
use rayon::prelude::*;
use std::{ops::Mul as _, rc::Rc, thread::available_parallelism};
use tracing::{Level, debug, instrument, trace, warn};

pub mod denovo_cpg;
pub mod methylation;
pub mod ml;
pub mod pileup;
mod record_filters;
pub mod variant_calling;
mod writer;

pub use record_filters::RecordFilters;
pub use writer::writer_thread;

// Jump in here if you want to know how the processing of regions works
pub mod process;

#[cfg(test)]
pub mod test_helpers;
#[cfg(test)]
pub mod tests;

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
    #[arg(short='@', long = "threads", env = "RASTAIR_THREADS", default_value_t = available_parallelism().map(|n|n.get()).unwrap_or(2).max(1))]
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

        if self.bed.bed.is_some() && self.vcf.vcf.is_none() {
            debug!("Only BED output requested, filtering for CpG/de-novo CpG sites only");
            self.record_filters.cpgs_only = true;
        }

        Ok(())
    }
}

/// Read BAM + FASTA and call variants and methylation events
#[instrument(level = "debug", skip(params))]
pub fn call(mut params: CallParams) -> Result<()> {
    params.figure_out_outputs().wrap_err("Unclear output choice")?;
    params.segmentation.sanitize();
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
        writer_thread(params, &regions, vcf_receiver).wrap_err("VCF writer error")?;

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
                    // This is where the actual processing happens!
                    process_region_wrapper(index, region, vcf_sender, params, &ml)
                },
            )
        })
        .wrap_err("Failed to process regions in parallel")
        .note("Rastair might have still written an (incomplete) output file")?;

    writer_thread
        .join()
        .map_err(|_| eyre!("Writer thread crashed"))
        .this_is_a_bug()? // this error is a panic in the thread
        .wrap_err("Error in writer thread")?; // this error is from actual result returned by the thread

    Ok(())
}

/// Wrapper function for processing a region in a thread-safe manner.
///
/// Calls [`process_region`] with thread-local readers and ships the result to
/// the VCF writer.
#[instrument(level = "info", skip_all, fields(region=%region.region))]
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
    let records = READERS.with(|local_readers| -> Result<Vec<PileupMetrics>> {
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
            local_readers
                .as_mut()
                .wrap_err("Failed to access thread-local resources")
                .this_is_a_bug()?
        };

        // This is the actual processing of the region

        // NOTE: There are some filters applied here to ignore certain reads.
        let pileup_mapping_params =
            process::PileupMappingParams { variant_calling: params.variant_calling.clone() };
        let (segment, pileups) = get_pileups(readers, region, &pileup_mapping_params)?;

        let res = process_region(segment, pileups, params, ml);

        // Handle processing errors gracefully to not crash the whole processing
        match res {
            Ok(records) => Ok(records),
            Err(e) => {
                warn!(error = format!("{e:#}"), "Failed to process region");
                // We still send an empty vector to the channel to increment the index
                Ok(Vec::new())
            }
        }
    })?;

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
    segment: Rc<Segment>,
    pileups: impl Iterator<Item = Pileup>,
    params: &CallParams,
    ml: &MachineLearning,
) -> Result<Vec<PileupMetrics>> {
    // Calculate metrics for each pileup.
    let threshold_filters = process::ThresholdFilterParams {
        variant_calling: params.variant_calling.clone(),
        methylation: params.methylation.thresholds.clone(),
        denovo_cpg: params.denovo_cpg.clone(),
    };

    macro_rules! log_failed_and_skip {
        ($msg:expr) => {
            |x: Result<PileupMetrics>| match x {
                Err(e) => {
                    warn!(error = format!("{e:#}"), $msg);
                    None
                }
                Ok(x) => Some(x),
            }
        };
    }

    let pileups: Vec<PileupMetrics> = calculate_pileup_metrics(pileups, &segment)
        .filter_map(log_failed_and_skip!("failed to calculate metric, skipping"))
        .map_surrounding(process::set_denovo_adj)
        .filter_map(log_failed_and_skip!("failed to set denovo adjacency, skipping"))
        .map(|mut current| {
            current.pos_metrics.extended.methylation_strand_info =
                MethylationEvidenceStrandInfo::from_pileup(&current);
            current
        })
        .filter(|p| params.record_filters.pre_filter(p))
        .map_surrounding(|b, c, a| {
            // More filters: Add ML metrics if requested
            process::add_ml_metrics(b, c, a, ml)
        })
        .filter_map(log_failed_and_skip!("failed to calculate ML score, skipping"))
        .map(|mut pileup| {
            // Add 'simple' filters based on the collected metrics
            process::apply_threshold_filters(&mut pileup, &threshold_filters)
                .wrap_err("Failed to apply threshold filters")?;
            Ok(pileup)
        })
        .filter_map(log_failed_and_skip!("failed to add threshold filters, skipping"))
        .map_surrounding(|b, c, a| {
            // For CpG sites and de-novo CpG sites, if one position is pass, mark
            // corresponding as pass as well
            process::propagate_denovo_pass_flags(b, c, a, params.ml.threshold())
        })
        .filter_map(log_failed_and_skip!("failed to propagate CpG pass flags, skipping"))
        .map(|mut pileup| {
            // Finally, set the actual variant calls based on all metrics and filters
            process::set_alt_calls(&mut pileup, params.ml.threshold())?;
            process::add_position_tags(&mut pileup);
            Ok(pileup)
        })
        .filter_map(log_failed_and_skip!("failed to set alt calls, skipping"))
        .map(|mut pileup| {
            // Set "extended" metrics that depend on the segment and params. This is
            // done in a separate step since it uses the pileup as well as the ML score.
            pileup.pos_metrics.extended.genotype =
                pileup.estimate_genotype(params.ml.threshold(), params.variant_calling.error_model);
            pileup.pos_metrics.extended.methylated =
                metrics::methylation::call(&pileup)?.unwrap_or_default();
            Ok(pileup)
        })
        .filter_map(log_failed_and_skip!("failed to calculate extended metrics, skipping"))
        .filter(|p| only_core_positions(&segment, p))
        .collect();

    // At this point, we have collected all metrics for the pileups in this
    // region. The recipient is responsible for further filtering based on
    // filters and writing them to the VCF or BED file.

    if tracing::enabled!(Level::DEBUG) {
        if pileups.is_empty() {
            debug!("No relevant pileups found in region");
        } else {
            let count_piles = readable::num::Unsigned::from(pileups.len());
            let pile_size = pileups.len() * std::mem::size_of::<PileupMetrics>();
            let read_size = pileups.iter().map(|p| p.pileup.reads.len()).sum::<usize>()
                * std::mem::size_of::<SimpleRead>();
            let bytes = readable::byte::Byte::from(pile_size + read_size);
            debug!(%count_piles, %bytes, "Collected pileup metrics");
        }
    }

    Ok(pileups)
}

fn only_core_positions(segment: &Segment, p: &PileupMetrics) -> bool {
    let pos = u64::from(p.pos());
    let core_start = segment.region.start + segment.overlap_start;
    let core_end = segment.region.end.saturating_sub(segment.overlap_end);

    pos >= core_start && pos < core_end
}
