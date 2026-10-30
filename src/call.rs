//! `rastair call`.
//!
//! The genome is cut into overlapping segments that worker threads process in
//! parallel; the calling thread writes their records in segment order. Per
//! segment (`process_segment`):
//!
//! 1. `segment_columns` builds each column from seqair, drops the ones with
//!    no alt, indel or CpG, finishes the rest into [`PileupMetrics`], calls
//!    the indels, and marks the CpGs called indels make.
//! 2. The forests score every alt; the SNV engine scores the ones they kept,
//!    and the engine forests score those again.
//! 3. `process_collected_pileups`: alt calls with the engine's genotype, CpGs
//!    of called alts, `0/0` genotypes, methylation, phasing.
//!
//! The writer emits the records the record filters keep as VCF/BCF and BED.

#[cfg(test)]
use crate::call::process::calculate_pileup_metrics;
#[cfg(test)]
use crate::call::tests::pileup::Pileup;
use crate::runtime::{fault_injection, progress, segments, threads};
use crate::{
    bed::rastair1::BedParams,
    call::{
        methylation::MethylationCallingParams, require_tags::RequireTagsParams,
        variant_calling::VariantCallingParams,
    },
    io::vcf_writer,
    metrics::{
        self, Alt, PileupMetrics,
        ml::types::{Cascade, MachineLearning},
    },
    sequence::{ChunkRegion, ReaderParams, Segment, SegmentationParams, SeqairReaders},
    utils::{cli, map_surrounding},
};
use clio::ClioPath;
use color_eyre::{
    Section,
    eyre::{Result, WrapErr, ensure},
};
use seqair_types::Base;
use std::{rc::Rc, thread::available_parallelism};
use tracing::{Level, debug, instrument, trace, warn};

pub mod methylation;
pub mod ml;
pub mod pair_hmm;
pub mod phasing;
pub mod pileup;
mod record_filters;
pub(crate) mod require_tags;
pub mod variant_calling;
mod writer;

#[cfg(test)]
pub use record_filters::pre_filter;
pub use record_filters::{ColumnEvidence, RecordFilters};
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
    pub methylation: MethylationCallingParams,
    #[command(flatten)]
    pub ml: ml::MachineLearningParams,
    #[command(flatten)]
    #[serde(flatten)]
    pub phasing: phasing::PhasingParams,

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
    // NOTE: There are some filters applied here to ignore certain reads.
    pub(crate) fn pileup_mapping(&self) -> process::PileupMappingParams {
        process::PileupMappingParams {
            variant_calling: self.variant_calling.clone(),
            require_tags: self.require_tags.filter(),
            engine: process::Engine::PairHmm,
            indel_floor: self.indel.min_indel_ao,
            segment_max_bytes: self.segmentation.segment_max_bytes,
            rescue_soft_clip_cpg: self.methylation.rescue_soft_clip_cpg,
            columns: process::ColumnKeep::Candidates,
            uncovered_cpgs: self.record_filters.uncovered_cpgs(),
            gpu: self.ml.gpu,
        }
    }

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
    let mut readers =
        params.segments.pileup_readers().wrap_err("Failed to read BAM/FASTA files")?;
    readers.run_reads_per_base = readers.mapped_reads_per_base();
    if readers.run_reads_per_base.is_none() {
        warn!(
            "The alignment index has no read counts (a CRAM, or no .bai/.csi), so the phasing \
             depth guard compares each segment with itself and misses repeats longer than one \
             segment"
        );
    }

    // Get segments that are small enough to process in RAM
    let regions: Vec<ChunkRegion> = readers
        .segments(params.segmentation.segment_max_length, params.segmentation.segment_overlap)
        .wrap_err("Could not fetch segments from BAM file")?
        .collect();
    let Some(first) = regions.first() else {
        warn!("No segments found in BAM file, nothing to do");
        return Ok(());
    };
    if first.region.len() < 2 {
        warn!(region=%first.region, "Given range is one base long, this will not yield any results for context-specific methylation calling.");
    }

    debug!("Going to process {} segments", regions.len());

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
    readers: &mut SeqairReaders,
    params: &CallParams,
    ml: &Cascade,
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

/// Score the alts `only` takes with the forests, on the inference thread if
/// there is one.
fn score_forest(
    pileups: &mut [PileupMetrics],
    ml: &MachineLearning,
    only: impl Fn(&PileupMetrics, Base) -> bool + Copy,
) -> Result<()> {
    // One batch per model either way; see `score_on_cpu` for why that
    // matters on the CPU.
    match process::score_on_gpu(pileups, ml, only) {
        process::GpuScoring::Scored => {}
        process::GpuScoring::LeftToCaller => process::score_on_cpu(pileups, ml, only)?,
    }
    Ok(())
}

/// The actual processing of a segment
fn process_segment(
    readers: &mut SeqairReaders,
    region: &ChunkRegion,
    params: &CallParams,
    ml: &Cascade,
) -> Result<Vec<PileupMetrics>> {
    fault_injection::fault_point(fault_injection::FaultPoint::Worker)?;

    let mapping = params.pileup_mapping();
    let SegmentColumns { segment, mut columns, held } = segment_columns(readers, region, &mapping)?;
    let depth_scale = readers.depth_scale(held.depth);
    score_forest(&mut columns, &ml.forests, |_, _| true)?;
    // The cascade: the engine scores the alts the forest kept, which gives
    // them their genotype, and the engine forests score exactly those again
    // with the engine's features. The others keep their first score.
    held.score_snv_clusters(&segment, &mut columns, &mapping);
    if let Some(engine) = &ml.engine {
        score_forest(&mut columns, engine, |p, alt| p.engine_scored(alt))?;
    }
    process_collected_pileups(segment, columns, params, depth_scale)
}

/// Every column of a segment, overlaps included, carrying everything the
/// forests read: `call` filters and scores these, `ml train` labels them, so a
/// quantity added here reaches both.
pub(crate) struct SegmentColumns {
    pub segment: Rc<Segment>,
    pub columns: Vec<PileupMetrics>,
    /// The reads the SNV engine scores the columns' SNV clusters on, under
    /// [`process::Engine::PairHmm`].
    pub held: process::HeldReads,
}

pub(crate) fn segment_columns(
    readers: &mut SeqairReaders,
    region: &ChunkRegion,
    mapping: &process::PileupMappingParams,
) -> Result<SegmentColumns> {
    let (segment, mut columns, held) =
        process::get_pileups_held(readers, region, mapping).wrap_err("Failed to build pileups")?;
    // The indels are called, and the forests read the CpGs they make.
    map_surrounding(&mut columns, process::mark_called_cpgs);
    Ok(SegmentColumns { segment, columns, held })
}

/// The value of a calling step that succeeded on a column; a failed one is
/// logged with `context` and its column dropped.
fn log_failed_and_skip<T>(context: &str) -> impl Fn(Result<T>) -> Option<T> + '_ {
    move |result| result.inspect_err(|error| warn!(error = format!("{error:#}"), "{context}")).ok()
}

/// Analyse hand-built pileups in a region as `call` does the reader's.
#[cfg(test)]
fn process_region(
    segment: Rc<Segment>,
    pileups_iter: impl Iterator<Item = Pileup>,
    params: &CallParams,
    ml: &Cascade,
) -> Result<Vec<PileupMetrics>> {
    let pileups: Vec<Pileup> = pileups_iter.collect();
    let reads = tests::utils::read_lines(&pileups);
    let mut pileups = metrics_of(pileups.into_iter(), &segment);
    pre_filter(&mut pileups);
    score_forest(&mut pileups, &ml.forests, |_, _| true)?;
    tests::utils::score_snv_lines(&mut pileups, &reads, &segment);
    if let Some(engine) = &ml.engine {
        // The stand-in scores alts by read count; production's cascade takes
        // only those the first stage scored at the floor.
        let admitted = |p: &PileupMetrics, alt: Base| {
            p.engine_scored(alt)
                && p.alts.iter().any(|a| {
                    a.base == alt
                        && a.verdict.ml().is_some_and(|ml| {
                            *ml >= crate::call::pileup::from_seqair::CASCADE_MIN_SCORE
                        })
                })
        };
        score_forest(&mut pileups, engine, admitted)?;
    }
    process_collected_pileups(segment, pileups, params, phasing::DepthScale::SegmentOnly)
}

#[cfg(test)]
fn metrics_of(pileups: impl Iterator<Item = Pileup>, segment: &Segment) -> Vec<PileupMetrics> {
    calculate_pileup_metrics(pileups, segment)
        .filter_map(log_failed_and_skip("failed to calculate metric, skipping"))
        .collect()
}

fn process_collected_pileups(
    segment: Rc<Segment>,
    pileups: Vec<PileupMetrics>,
    params: &CallParams,
    depth_scale: phasing::DepthScale,
) -> Result<Vec<PileupMetrics>> {
    // Before the overlap is trimmed: a site in it links reads reaching into
    // the core.
    let phase = |pileups: &mut [PileupMetrics]| {
        phasing::phase_segment(
            pileups,
            |pos| segment.is_core(pos),
            depth_scale,
            params.phasing.phase_min_switch_quality,
        );
    };
    let mut pileups = call_columns(
        pileups,
        params.ml.ml,
        params.methylation.m_max_coverage,
        params.variant_calling.error_model,
        phase,
    );
    pileups.retain(|p| only_core_positions(&segment, p));

    // At this point, we have collected all metrics for the pileups in this
    // region. The recipient is responsible for further filtering based on
    // filters and writing them to the VCF or BED file.

    if tracing::enabled!(Level::DEBUG) {
        if pileups.is_empty() {
            debug!("No relevant pileups found in region");
        } else {
            let count_piles = readable::num::Unsigned::from(pileups.len());
            let bytes =
                readable::byte::Byte::from(pileups.len() * std::mem::size_of::<PileupMetrics>());
            debug!(%count_piles, %bytes, "Collected pileup metrics");
        }
    }

    Ok(pileups)
}

/// The calls on scored columns: each alt's with the engine's genotype, then
/// `phase`, then the CpGs called alts make (which a pair's phase decides
/// where no cluster holds both), then the `0/0` genotype of a column with no
/// called alt (a de-novo partner's reads weigh as a CpG's), and methylation.
fn call_columns(
    mut pileups: Vec<PileupMetrics>,
    ml_threshold: seqair_types::Probability,
    m_max_coverage: std::num::NonZeroU32,
    error_model: variant_calling::ErrorModel,
    phase: impl FnOnce(&mut [PileupMetrics]),
) -> Vec<PileupMetrics> {
    for pileup in &mut pileups {
        process::set_alt_calls(pileup, ml_threshold, m_max_coverage);
    }
    phase(&mut pileups);
    map_surrounding(&mut pileups, process::mark_called_cpgs);

    pileups
        .into_iter()
        .map(|mut pileup| {
            if !pileup.alts.iter().any(Alt::called) {
                pileup.pos_metrics.genotype = pileup.hom_ref_genotype(error_model);
            }
            pileup.pos_metrics.methylated = metrics::methylation::call(&pileup)?;
            Ok(pileup)
        })
        .filter_map(log_failed_and_skip("failed to call the column, skipping"))
        .collect()
}

fn only_core_positions(segment: &Segment, p: &PileupMetrics) -> bool {
    segment.is_core(p.pos)
}
