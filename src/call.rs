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
//!   1. The calling thread receives the records of each segment, in order
//!   2. Filter be given criteria
//!   3. Convert to the output format
//!   4. Write to file in order

#[cfg(any(not(feature = "experimental-seqair"), test))]
use crate::call::pileup::Pileup;
#[cfg(any(not(feature = "experimental-seqair"), test))]
use crate::call::process::calculate_pileup_metrics;
use crate::runtime::{fault_injection, progress, segments, threads};
use crate::{
    bed::rastair1::BedParams,
    call::{
        methylation::params::MethylationCallingParams, pileup::SimpleRead, process::get_pileups,
        require_tags::RequireTagsParams, variant_calling::VariantCallingParams,
    },
    io::vcf_writer,
    metrics::{self, MethylationEvidenceStrandInfo, PileupMetrics, ml::types::MachineLearning},
    sequence::{
        ChunkRegion, PileupReaders, ReaderParams, ReaderSource, Segment, SegmentationParams,
    },
    utils::{cli, map_surrounding},
};
use clio::ClioPath;
use color_eyre::{
    Section,
    eyre::{Result, WrapErr, ensure},
};
use std::{rc::Rc, thread::available_parallelism};
use tracing::{Level, debug, instrument, trace, warn};

pub mod denovo_cpg;
pub mod methylation;
pub mod ml;
pub mod pileup;
mod record_filters;
pub(crate) mod require_tags;
pub mod variant_calling;
mod writer;

pub use record_filters::{PreFilterInputs, RecordFilters};
use writer::SegmentWriter;

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
    #[command(flatten)]
    #[serde(flatten)]
    pub require_tags: RequireTagsParams,

    // --- Calling parameters ---
    #[command(flatten)]
    pub variant_calling: VariantCallingParams,
    #[command(flatten)]
    pub indel: variant_calling::IndelParams,
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

    #[cfg(not(feature = "experimental-seqair"))]
    if params.methylation.rescue_soft_clip_cpg {
        warn!(
            "--rescue-soft-clip-cpg has no effect on the default (htslib) backend; \
             build with the `experimental-seqair` feature to use it"
        );
    }

    let params = &params; // make params immutable for threads

    // Initialize readers for BAM and FASTA files
    let readers = params.segments.pileup_readers().wrap_err("Failed to read BAM/FASTA files")?;

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

    debug!("Going to process {} segments", regions.len());

    let readers = ReaderSource::from(readers);
    let readers = &readers;

    progress::register_signal_handler();
    let panics = threads::PanicCheck::start();

    // Process each region and write results to the VCF
    //
    // This is done in parallel to speed up the processing, so here are a few
    // comments on this works. There are two aspects to this: Collecting the
    // variant candidates and writing them to the VCF.
    //
    // To use all CPU available, we use rayon to process the regions in
    // parallel. The ready-made records come back to this thread in segment
    // order, which writes them.
    let writer_threads = params.vcf.vcf_threads;
    let worker_threads = params.total_threads.saturating_sub(writer_threads.get()).max(1);

    // Needs the worker count to size the inference queue, so it cannot be built
    // before now.
    let ml =
        params.ml.init(worker_threads).wrap_err("Failed to initialize machine learning model")?;

    debug!(
        "Gonna use {} threads: {} for processing, {} for writing VCF",
        params.total_threads,
        worker_threads,
        writer_threads.get()
    );

    let mut writer = SegmentWriter::new(params, &regions).wrap_err("VCF writer error")?;
    let failed_segments = &segments::FailedSegments::default();

    // Run this in a custom rayon thread pool to control the number of threads
    // and be able to tweak parameters when profiling
    let pool = rayon::ThreadPoolBuilder::new()
        .thread_name(|idx| format!("worker-{idx}"))
        .num_threads(worker_threads)
        .start_handler(|idx| trace!(idx, "Starting worker thread"))
        .exit_handler(|idx| trace!(idx, "Closing worker thread"))
        .build()
        .wrap_err("Failed to create thread pool for rayon")?;
    let processed = segments::process_in_order(
        &pool,
        &regions,
        || readers.fork().wrap_err("Failed to open readers in worker thread"),
        |readers, region| Ok(process_region_wrapper(region, readers, params, &ml, failed_segments)),
        |records| writer.write(records),
    )
    .wrap_err("Failed to process regions in parallel")
    .note("Output files are incomplete and were left under their `.partial` names");

    // Close the outputs even if processing failed, so they hold a valid prefix
    let (closed, outputs) = writer.close();
    segments::first_error(processed, closed)?;
    panics.finish()?;
    failed_segments.check(regions.len())?;
    for output in outputs {
        output.commit()?;
    }

    Ok(())
}

/// Processes one segment with the worker's readers. A segment that fails to
/// process is counted and yields no records; the run goes on.
#[instrument(level = "info", skip_all, fields(region=%region.region))]
fn process_region_wrapper(
    region: &ChunkRegion,
    readers: &mut PileupReaders,
    params: &CallParams,
    ml: &MachineLearning,
    failed_segments: &segments::FailedSegments,
) -> Vec<PileupMetrics> {
    let mut records = process_segment(readers, region, params, ml).unwrap_or_else(|error| {
        failed_segments.record(&region.region, error);
        Vec::new()
    });
    // The pipeline's in-place collects keep the capacity sized for every
    // covered base, about five times what survives the filters, and this vec
    // may wait behind a slow segment before the writer frees it.
    records.shrink_to_fit();
    records
}

/// The actual processing of a segment
fn process_segment(
    readers: &mut PileupReaders,
    region: &ChunkRegion,
    params: &CallParams,
    ml: &MachineLearning,
) -> Result<Vec<PileupMetrics>> {
    fault_injection::fault_point(fault_injection::FaultPoint::Worker)?;

    // NOTE: There are some filters applied here to ignore certain reads.
    let pileup_mapping_params = process::PileupMappingParams {
        variant_calling: params.variant_calling.clone(),
        require_tags: params.require_tags.filter(),
        call_indels: params.indel.enabled(),
        indel_max_mismatches: params.indel.indel_max_mismatches,
        indel_end_of_read_cutoff: params.indel.indel_end_of_read_cutoff,
        segment_max_bytes: params.segmentation.segment_max_bytes,
        rescue_soft_clip_cpg: params.methylation.rescue_soft_clip_cpg,
        early_reject: Some(params.record_filters.clone()),
        ..Default::default()
    };

    #[cfg(not(feature = "experimental-seqair"))]
    {
        let (segment, pileups) = get_pileups(readers, region, &pileup_mapping_params)?;
        process_region(segment, pileups, params, ml)
    }
    #[cfg(feature = "experimental-seqair")]
    {
        let (segment, metrics) = get_pileups(readers, region, &pileup_mapping_params)?;
        process_pre_built_metrics(segment, metrics, params, ml)
    }
}

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

#[cfg(any(not(feature = "experimental-seqair"), test))]
/// Analyse pileups in a region
fn process_region(
    segment: Rc<Segment>,
    pileups_iter: impl Iterator<Item = Pileup>,
    params: &CallParams,
    ml: &MachineLearning,
) -> Result<Vec<PileupMetrics>> {
    // One per covered position, as on the seqair path — see the note in
    // `process::pileups::get_pileups`. Doubling into a vec of 592-byte entries
    // both copies and overshoots.
    let mut pileups: Vec<PileupMetrics> =
        Vec::with_capacity(usize::try_from(segment.range.len()).unwrap_or(0));
    pileups.extend(
        calculate_pileup_metrics(pileups_iter, &segment)
            .filter_map(log_failed_and_skip!("failed to calculate metric, skipping")),
    );
    map_surrounding(
        &mut pileups,
        process::set_denovo_adj,
        "failed to set denovo adjacency, skipping",
    );
    set_strand_info_and_prefilter(&mut pileups, params);

    process_collected_pileups(segment, pileups, params, ml)
}

#[cfg(feature = "experimental-seqair")]
fn process_pre_built_metrics(
    segment: Rc<Segment>,
    pileups: impl Iterator<Item = PileupMetrics>,
    params: &CallParams,
    ml: &MachineLearning,
) -> Result<Vec<PileupMetrics>> {
    let mut pileups: Vec<PileupMetrics> = pileups.collect();
    map_surrounding(
        &mut pileups,
        process::set_denovo_adj,
        "failed to set denovo adjacency, skipping",
    );
    set_strand_info_and_prefilter(&mut pileups, params);
    process_collected_pileups(segment, pileups, params, ml)
}

/// The step both backends run between de-novo adjacency and the shared
/// pipeline: strand info depends on the adjacency flags just set, and the
/// pre-filter on the alts, so the order matters.
fn set_strand_info_and_prefilter(pileups: &mut Vec<PileupMetrics>, params: &CallParams) {
    for pileup in pileups.iter_mut() {
        pileup.pos_metrics.extended.methylation_strand_info =
            MethylationEvidenceStrandInfo::from_pileup(pileup);
    }
    pileups.retain(|p| params.record_filters.pre_filter(p));
}

fn process_collected_pileups(
    segment: Rc<Segment>,
    mut pileups: Vec<PileupMetrics>,
    params: &CallParams,
    ml: &MachineLearning,
) -> Result<Vec<PileupMetrics>> {
    let threshold_filters = process::ThresholdFilterParams {
        variant_calling: params.variant_calling.clone(),
        methylation: params.methylation.thresholds.clone(),
        denovo_cpg: params.denovo_cpg.clone(),
    };

    if params.indel.enabled() {
        for p in &mut pileups {
            if let Some(ref mut d) = p.indel_data {
                let tract = u32::from(d.homopolymer_run.max(d.dinucleotide_run));
                d.calls = variant_calling::indel_calling::call_indels(
                    &d.counts,
                    &params.indel,
                    params.indel.use_ml(ml.enabled()),
                    tract,
                    params.indel.rescues_hom_ref(),
                );
            }
        }
    }

    // Pass 2: ML prediction on the inference thread when there is one, and on
    // this thread otherwise or if the GPU failed. Both score the region as one
    // batch per model; see `score_on_cpu` for why that matters on the CPU.
    let score_indels = params.indel.needs_ml_scores(ml.enabled());
    if !process::score_on_gpu(&mut pileups, ml, score_indels) {
        process::score_on_cpu(&mut pileups, ml, score_indels)?;
    }

    if params.indel.rescues_hom_ref() {
        for p in &mut pileups {
            if let Some(ref mut d) = p.indel_data {
                variant_calling::indel_calling::rescue_hom_ref(&mut d.calls, params.ml.threshold());
            }
        }
    }

    let mut pileups: Vec<PileupMetrics> = pileups
        .into_iter()
        .map(|mut pileup| {
            process::apply_threshold_filters(&mut pileup, &threshold_filters)
                .wrap_err("Failed to apply threshold filters")?;
            Ok(pileup)
        })
        .filter_map(log_failed_and_skip!("failed to add threshold filters, skipping"))
        .collect();
    // For CpG sites and de-novo CpG sites, if one position is pass, mark
    // corresponding as pass as well
    map_surrounding(
        &mut pileups,
        |b, c, a| process::propagate_denovo_pass_flags(b, c, a, params.ml.threshold()),
        "failed to propagate CpG pass flags, skipping",
    );

    let pileups: Vec<PileupMetrics> = pileups
        .into_iter()
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
            pileup.pos_metrics.extended.methylation_strand_info =
                MethylationEvidenceStrandInfo::from_pileup_with_methylation(&pileup);

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
            let read_size = pileups.iter().map(|p| p.pos_metrics.depth as usize).sum::<usize>()
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
