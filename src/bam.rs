//! `rastair bam`: rewrite a BAM with per-read methylation tags.
//!
//! Two modes: `legacy` writes XR/XG/XM and leaves SEQ alone; `standard` writes
//! MM/ML and rewrites SEQ (T→C on OT, A→G on OB). Tags are built in
//! [`base_modification`].
//!
//! **XM and MM/ML encode per-read methylation, not position-level calls.** A
//! CpG with beta 0.3 (`methylated: Some(false)` in [`RastairCall::Cpg`]) still
//! has methylated reads that must show `Z` in XM and appear in MM/ML; the call
//! only decides whether the position is a CpG site, and a CpG without a beta
//! (`methylated: None`) is one too. Per-read methylation comes from the
//! observed base alone: T at an OT C, A at an OB G
//! ([`CpgSide::read_methylated`], which `per-read` shares).
//!
//! **MM/ML and XM count different reads.** MM/ML only encode modifications at
//! C bases of the stored SEQ, and of two mates overlapping a CpG only one has
//! a C there (the other has the G), while XM annotates both. So MM/ML readers
//! such as modkit see about half the reads XM counting does: fractions agree,
//! counts differ ~2:1. Compare the two modes by fraction, at a minimum
//! coverage.
//!
//! **The input is read front to back**, so every record comes out, in input
//! order, unplaced ones included, and no index is needed. A reader thread
//! (`Source`) cuts it into `Batch`es: `BATCH_RECORDS` records of BAM or
//! SAM, or one CRAM container, which the worker that takes it decodes, so
//! CRAM decoding runs in parallel. A worker splits its batch into runs of
//! records within `LOOKUP_MAX_SPAN` bases of one contig (unplaced records join
//! any run), looks up the calls (and, for de-novo CpGs in legacy mode, the
//! reference) over each run's span once, and rewrites the records in place;
//! the writer writes the batches in order ([`segments::process_in_order`]).
//! Unsorted input is rewritten correctly but slowly: read by name, nearly
//! every record is a run of its own, with its own calls query (and FASTA
//! fetch in legacy mode). Sort by coordinate first for speed.
//!
//! `--region` reads the records overlapping it through the index (BAI, CSI or
//! CRAI), each once, also when it overlaps two regions.
//!
//! A record whose aux block is malformed fails the run, naming the read.
//!
//! A record of unknown strand (paired, but neither first nor second) carries
//! no methylation; an unmapped or unplaced one has no aligned bases, so it gets
//! the tags of a read without methylation evidence.

use crate::{
    bed::reader::{RastairBedReader, RastairCall},
    sequence::{ReaderParams, Region},
    utils::logging::ThisIsABug as _,
    utils::{
        CliRegionInput, cli,
        file_helpers::{FastaReader, open_fasta},
    },
    vcf::CpgSide,
};
use clap::{Parser, value_parser};
use clio::ClioPath;
use color_eyre::{
    Section as _,
    eyre::{Context, ContextCompat as _, Report, Result, ensure},
};
use rustc_hash::FxHashMap;
use seqair::{
    bam::{
        BamHeader, RawBamRecord, aligned_pairs::MatchPosition, header::PgRecord,
        writer::BamWriterBuilder,
    },
    cram::stream::{CramContainer, CramContainerDecoder, CramStreamReader},
    io::CompressionLevel,
    reader::{AlignmentStream, CramInput, IndexedReader, Tid},
};
use seqair_types::{Base, IupacBase, Pos0, QPos, SmallVec, Strand};
use std::{
    fmt,
    io::{self, Write as _},
    num::NonZeroU8,
    sync::{
        Arc,
        atomic::{AtomicU64, Ordering},
        mpsc::SyncSender,
    },
    thread::available_parallelism,
};

mod base_modification;
use crate::runtime::{
    fault_injection, partial_output::PartialOutput, progress, progress::ProgressTracker, segments,
    threads,
};
use base_modification::remove_modifications;
pub use base_modification::{
    MethylatedPositions, MethylationContext, XmAnnotation, XrTags, determine_context,
};
use tracing::instrument;

/// Subcommands for `rastair bam`
#[derive(Debug, clap::Subcommand)]
pub enum BamSubcommand {
    /// Write modBAM with MM/ML tags as specified by the SAM 4.5 spec
    /// This will rewrite SEQ to un-modify bases that have methylation evidence.
    Standard(BamRewriteArgs),
    /// Write BAM with "legacy" XR/XG/XM tags, compatible with tools like DRAGEN
    /// and Bismark.
    Legacy(BamRewriteArgs),
}

/// Controls whether SEQ is rewritten and which tags are produced
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BamMode {
    Standard,
    Legacy,
}

#[derive(Debug, Parser)]
pub struct BamRewriteArgs {
    #[command(flatten)]
    segments: ReaderParams,
    /// Number of threads: `N - 1` rewrite records and `N` more compress the
    /// output, as `samtools view -@ N` does; three more read and inflate the
    /// input.
    #[arg(short='@', long = "threads", default_value_t = available_parallelism().map(|n|n.get()).unwrap_or(2).max(1))]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub threads: usize,

    /// Rastair's calls to determine methylation
    #[arg(value_parser=value_parser!(ClioPath).exists().is_file(), value_hint=clap::ValueHint::FilePath)]
    #[arg(help_heading = cli::sections::INPUT)]
    calls_file: ClioPath,

    /// Output file
    #[arg(short = 'o', long, default_value = "-")]
    #[arg(help_heading = cli::sections::OUTPUT)]
    output: ClioPath,

    /// BGZF compression level of the output, from 0 (uncompressed) to 12
    /// (smallest, slowest). Higher levels cost the compression threads (`-@`)
    /// more time; samtools writes 6.
    #[arg(long, default_value_t = 1, value_parser = value_parser!(u8).range(0..=12))]
    #[arg(help_heading = cli::sections::OUTPUT)]
    compression_level: u8,
}

/// Records per [`Batch`] of BAM or SAM: about what a 10 kb segment held at
/// 30x. A CRAM batch is one container.
const BATCH_RECORDS: usize = 2_000;
/// The longest reference span one calls and reference lookup covers.
const LOOKUP_MAX_SPAN: u64 = 1_000_000;

/// [`BATCH_RECORDS`]; debug builds take `RASTAIR_BAM_BATCH_RECORDS`, so tests
/// cross batch boundaries with small inputs.
fn batch_records() -> usize {
    static SIZE: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *SIZE.get_or_init(|| test_size("RASTAIR_BAM_BATCH_RECORDS").unwrap_or(BATCH_RECORDS))
}

/// [`LOOKUP_MAX_SPAN`]; debug builds take `RASTAIR_BAM_LOOKUP_SPAN`, so tests
/// cross run boundaries on short contigs.
fn lookup_max_span() -> u64 {
    static SPAN: std::sync::OnceLock<u64> = std::sync::OnceLock::new();
    *SPAN.get_or_init(|| test_size("RASTAIR_BAM_LOOKUP_SPAN").unwrap_or(LOOKUP_MAX_SPAN))
}

fn test_size<T: std::str::FromStr>(var: &str) -> Option<T> {
    if cfg!(debug_assertions) { std::env::var(var).ok()?.parse().ok() } else { None }
}
/// Threads inflating the input's BGZF blocks.
const INFLATE_THREADS: NonZeroU8 = NonZeroU8::MIN.saturating_add(1);

#[tracing::instrument(level = "info", skip_all, fields(
    ?mode,
))]
pub fn rewrite(params: &BamRewriteArgs, mode: BamMode) -> Result<()> {
    progress::register_signal_handler();
    let panics = threads::PanicCheck::start();

    let Input { source, header: input_header, decoder } = Input::open(&params.segments)?;
    let mut header = BamHeader::from_template(&input_header);
    header.add_pg(rastair_pg()).wrap_err("Failed to add rastair's @PG line to the header")?;

    let worker_threads = params.threads.saturating_sub(1).max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .thread_name(|idx| format!("bam-worker-{idx}"))
        .num_threads(worker_threads)
        .build()
        .wrap_err("Failed to create thread pool for BAM rewrite")?;

    let output = PartialOutput::new("BAM output", &params.output);
    let output_file = output.write_path()?;
    let sink: Box<dyn io::Write> = if output_file.is_std() {
        Box::new(io::stdout().lock())
    } else {
        let file = std::fs::File::create(output_file.path())
            .wrap_err_with(|| format!("Failed to create `{}`", output_file.path().display()))?;
        Box::new(io::BufWriter::new(file))
    };
    let mut writer = BamWriterBuilder::to_writer(sink, &header)
        .compression_level(
            CompressionLevel::try_from(i32::from(params.compression_level))
                .wrap_err("invalid `--compression-level`")?,
        )
        // Level-1 BGZF takes about two thirds of the rewrite's CPU, so a fixed
        // few threads capped the run; idle compressors only wait.
        .compression_threads(params.threads)
        .build()
        .wrap_err("failed to create writer")?;

    let (total, unit) = source.progress_total();
    let mut progress = ProgressTracker::counting(total, unit);
    let window = source.window(worker_threads);
    let processed = std::thread::scope(|scope| {
        let (send, batches) = std::sync::mpsc::sync_channel(worker_threads);
        let reader = std::thread::Builder::new()
            .name("bam-reader".to_owned())
            .spawn_scoped(scope, move || source.send_batches(&send))
            .wrap_err("Failed to start the BAM reader thread")?;
        let processed = segments::process_in_order(
            &pool,
            window,
            batches.into_iter().map(Pulled),
            || Worker::open(params, mode, decoder.as_ref()),
            |worker, Pulled(batch)| {
                fault_injection::fault_point(fault_injection::FaultPoint::Worker)?;
                let batch = batch?;
                let done = batch.progress;
                worker.rewrite(&input_header, batch, mode).map(|records| (records, done))
            },
            |(records, done)| {
                for record in &records {
                    writer.write_raw(record).wrap_err("failed to write record to new BAM file")?;
                }
                progress.done_up_to(done);
                fault_injection::fault_point(fault_injection::FaultPoint::Writer)
            },
        );
        // A stopped run drops the receiver, which ends the reader
        let read = match reader.join() {
            Ok(()) => Ok(()),
            Err(_) => {
                Err(color_eyre::eyre::eyre!("The BAM reader panicked (see the crash report above)"))
            }
        };
        segments::first_error(processed, read)
    })
    .wrap_err("Failed to process the BAM file in parallel");

    // Close the output even if processing failed, so it holds a valid prefix
    let closed = writer
        .finish()
        .map_err(Report::new)
        .and_then(|(mut sink, _)| sink.flush().map_err(Report::new))
        .wrap_err("Failed to write the BAM file, the disk may be full");
    segments::first_error(processed, closed)?;
    panics.finish()?;
    output.commit()?;

    Ok(())
}

/// Counts the bytes read from the input, for progress.
struct CountingReader {
    inner: std::fs::File,
    read: Arc<AtomicU64>,
}

impl io::Read for CountingReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.read.fetch_add(n as u64, Ordering::Relaxed);
        Ok(n)
    }
}

/// The input, opened: where its records come from, and what the workers
/// need to finish them.
struct Input {
    source: Source,
    header: BamHeader,
    /// For CRAM read front to back: the template each worker forks its
    /// decoder from
    decoder: Option<CramContainerDecoder>,
}

/// Where the records come from, in output order.
enum Source {
    /// BAM or SAM, front to back
    Stream { stream: AlignmentStream, input_bytes: u64, read: Arc<AtomicU64> },
    /// CRAM, front to back, as containers the workers decode
    Cram { reader: CramStreamReader<CramInput>, input_bytes: u64, read: Arc<AtomicU64> },
    /// `--region` by region through the index
    Regions { reader: IndexedReader, regions: Regions },
}

impl Input {
    fn open(params: &ReaderParams) -> Result<Self> {
        let path = params.bam_file.path();
        let open_error = || format!("Failed to read `{}`", path.display());
        if params.regions.is_some() {
            let reader = IndexedReader::open_with_reference(path, params.fasta_file.path())
                .wrap_err_with(open_error)
                .suggestion("`--region` reads through the index: `samtools index` creates it")?;
            let regions = Regions::new(params.regions.as_ref(), reader.header())?;
            let header = BamHeader::from_template(reader.header());
            return Ok(Self { source: Source::Regions { reader, regions }, header, decoder: None });
        }

        let file = std::fs::File::open(path).wrap_err_with(open_error)?;
        let input_bytes = file.metadata().wrap_err_with(open_error)?.len();
        let read = Arc::new(AtomicU64::new(0));
        let input = CountingReader { inner: file, read: Arc::clone(&read) };
        let inflate_pool = rayon::ThreadPoolBuilder::new()
            .thread_name(|idx| format!("bam-inflate-{idx}"))
            .num_threads(std::num::NonZeroUsize::from(INFLATE_THREADS).get())
            .build()
            .wrap_err("Failed to create thread pool for reading the BAM file")?;
        let stream = AlignmentStream::builder(input)
            .threads_in(INFLATE_THREADS, inflate_pool)
            .reference(params.fasta_file.path())
            // The names of CRAM records the file stores none for, as htslib
            // names them after the path it opened
            .name_prefix(
                path.file_name()
                    .map_or_else(|| "-".into(), |name| name.to_string_lossy().into_owned()),
            )
            .build()
            .wrap_err_with(open_error)?;
        let header = BamHeader::from_template(stream.header());
        // CRAM containers are decoded on the workers, BAM and SAM records come
        // decoded from the stream
        let (source, decoder) = match stream.into_cram() {
            Ok((reader, decoder)) => (Source::Cram { reader, input_bytes, read }, Some(decoder)),
            Err(stream) => (Source::Stream { stream, input_bytes, read }, None),
        };
        Ok(Self { source, header, decoder })
    }
}

impl Source {
    /// The progress the reader reports: input bytes, or bases of `--region`.
    fn progress_total(&self) -> (usize, &'static str) {
        let kilo = |n: u64| usize::try_from(n / 1000).unwrap_or(usize::MAX);
        match self {
            Self::Stream { input_bytes, .. } | Self::Cram { input_bytes, .. } => {
                (kilo(*input_bytes), "kB of input")
            }
            Self::Regions { regions, .. } => {
                (kilo(regions.0.iter().map(|span| span.end - span.start).sum()), "kb of regions")
            }
        }
    }

    /// How many batches may run ahead of the one being written: for CRAM one
    /// per worker, as a decoded container is often 10,000 records, five BAM
    /// batches at about 1 kB each in memory; four per worker otherwise.
    fn window(&self, workers: usize) -> Option<usize> {
        matches!(self, Self::Cram { .. }).then_some(workers)
    }

    /// Read the input to the end and send it as batches, or the error that
    /// stopped it; until the receiver hangs up.
    fn send_batches(self, send: &SyncSender<Result<Batch>>) {
        let mut number = 0_u64;
        let mut emit = |chunk: Chunk, progress: u64| {
            number += 1;
            let kilo = usize::try_from(progress / 1000).unwrap_or(usize::MAX);
            send.send(Ok(Batch { number, chunk, progress: kilo })).is_ok()
        };
        let result = match self {
            Self::Stream { mut stream, read, .. } => (|| {
                let mut records = Vec::with_capacity(batch_records());
                loop {
                    let mut record = RawBamRecord::default();
                    let more = match stream.read_record(&mut record) {
                        Ok(more) => more,
                        Err(error) => {
                            // The records before it are written, as a prefix
                            if !records.is_empty() {
                                emit(Chunk::Records(records), read.load(Ordering::Relaxed));
                            }
                            return Err(error).wrap_err("Failed to read the input BAM file");
                        }
                    };
                    if more {
                        records.push(record);
                    }
                    if (records.len() == batch_records() || !more) && !records.is_empty() {
                        let chunk = Chunk::Records(std::mem::replace(
                            &mut records,
                            Vec::with_capacity(batch_records()),
                        ));
                        if !emit(chunk, read.load(Ordering::Relaxed)) {
                            return Ok(());
                        }
                    }
                    if !more {
                        return Ok(());
                    }
                }
            })(),
            Self::Cram { mut reader, read, .. } => (|| {
                while let Some(container) =
                    reader.next_container().wrap_err("Failed to read the input CRAM file")?
                {
                    if !emit(Chunk::Container(container), read.load(Ordering::Relaxed)) {
                        break;
                    }
                }
                Ok(())
            })(),
            Self::Regions { mut reader, regions } => (|| {
                let mut records = Vec::with_capacity(batch_records());
                let mut done = 0;
                let mut previous: Option<Span> = None;
                for &region in &regions.0 {
                    if region.start >= region.end {
                        continue;
                    }
                    let pos = |at: u64| u32::try_from(at).ok().and_then(Pos0::new);
                    let (Some(first), Some(last)) =
                        (pos(region.start), region.end.checked_sub(1).and_then(pos))
                    else {
                        continue;
                    };
                    let mut query = reader
                        .query_raw(region.tid, (first..=last).into())
                        .wrap_err_with(|| format!("Failed to query the input for {region:?}"))?;
                    loop {
                        let mut record = RawBamRecord::default();
                        let more = match query.read_record(&mut record) {
                            Ok(more) => more,
                            Err(error) => {
                                // The records before it are written, as a prefix
                                if !records.is_empty() {
                                    emit(Chunk::Records(records), done);
                                }
                                return Err(error).wrap_err("Failed to read the input file");
                            }
                        };
                        let at = record.pos().map_or(region.start, |pos| pos.as_u64());
                        if more && !Regions::queried_before(previous, Span::of(&record)) {
                            records.push(record);
                        }
                        if (records.len() == batch_records() || !more) && !records.is_empty() {
                            let progress = done + at.clamp(region.start, region.end) - region.start;
                            let chunk = Chunk::Records(std::mem::replace(
                                &mut records,
                                Vec::with_capacity(batch_records()),
                            ));
                            if !emit(chunk, progress) {
                                return Ok(());
                            }
                        }
                        if !more {
                            break;
                        }
                    }
                    done += region.end - region.start;
                    previous = Some(region);
                }
                Ok(())
            })(),
        };
        if let Err(error) = result {
            // A hung-up receiver needs no error: the run already stopped
            send.send(Err(error)).ok();
        }
    }
}

/// Consecutive records of the input, the unit of work.
struct Batch {
    /// From 1, in input order
    number: u64,
    chunk: Chunk,
    /// The progress once it is written
    progress: usize,
}

/// A batch's records as the source hands them over.
enum Chunk {
    Records(Vec<RawBamRecord>),
    /// Decoded by the worker that takes it
    Container(CramContainer),
}

impl fmt::Display for Batch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.chunk {
            Chunk::Records(records) => {
                write!(f, "batch {} of the input ({} records)", self.number, records.len())
            }
            Chunk::Container(container) => {
                write!(f, "CRAM container {} of the input", container.number())
            }
        }
    }
}

/// What a worker pulls from the reader: a batch, or why there is none.
struct Pulled(Result<Batch>);

impl fmt::Display for Pulled {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.0 {
            Ok(batch) => batch.fmt(f),
            Err(_) => f.write_str("the rest of the input"),
        }
    }
}

/// A reference stretch, 0-based and half-open.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Span {
    tid: Tid,
    start: u64,
    end: u64,
}

impl Span {
    /// Where `record` aligns; `None` for an unplaced record.
    fn of(record: &RawBamRecord) -> Option<Self> {
        let (tid, pos) = (record.tid()?, record.pos()?);
        // A record that consumes no reference counts as covering its position,
        // as in htslib's region queries.
        let last = record.end_pos().unwrap_or(pos);
        Some(Self { tid, start: pos.as_u64(), end: last.as_u64().saturating_add(1) })
    }

    fn overlaps(self, other: Self) -> bool {
        self.tid == other.tid && self.start < other.end && other.start < self.end
    }

    /// The span covering both, if they share a contig and it stays short.
    fn joined(self, other: Self) -> Option<Self> {
        let joined = Self {
            tid: self.tid,
            start: self.start.min(other.start),
            end: self.end.max(other.end),
        };
        (self.tid == other.tid && joined.end - joined.start <= lookup_max_span()).then_some(joined)
    }
}

/// `--region`, as sorted, non-overlapping spans.
struct Regions(Vec<Span>);

impl Regions {
    fn new(input: Option<&CliRegionInput>, header: &BamHeader) -> Result<Self> {
        let mut spans = input
            .map(CliRegionInput::regions)
            .unwrap_or_default()
            .iter()
            .map(|region| {
                let tid = header.tid(&region.chromosome).wrap_err_with(|| {
                    format!("Chromosome {} of `--region` is not in the BAM file", region.chromosome)
                })?;
                let length = header.target_len(tid).wrap_err("Failed to fetch header length")?;
                let start = region.start.map_or(0, |p| p.to_zero_based().as_u64());
                // A closed 1-based end is the same number as a half-open 0-based one.
                let end = region.end.map_or(length, |p| p.as_u64());
                ensure!(
                    start <= length && end <= length,
                    "Region {region} extends past the end of chromosome {}",
                    region.chromosome
                );
                Ok(Span { tid, start, end })
            })
            .collect::<Result<Vec<_>>>()?;
        spans.sort_unstable();
        let mut merged: Vec<Span> = Vec::with_capacity(spans.len());
        for span in spans {
            match merged.last_mut() {
                Some(last) if last.tid == span.tid && span.start <= last.end => {
                    last.end = last.end.max(span.end);
                }
                _ => merged.push(span),
            }
        }
        Ok(Self(merged))
    }

    /// Whether a record a region's query returned came with the query of the
    /// region before, `previous`: regions are disjoint and sorted, so one
    /// reaching any earlier region reaches that one.
    fn queried_before(previous: Option<Span>, record: Option<Span>) -> bool {
        previous.zip(record).is_some_and(|(previous, record)| record.overlaps(previous))
    }
}

/// What each worker reads from, opened once on its own thread.
struct Worker {
    calls: RastairBedReader,
    /// Only legacy mode needs the reference
    fasta: Option<FastaReader>,
    /// For CRAM input
    decoder: Option<CramContainerDecoder>,
}

impl Worker {
    fn open(
        params: &BamRewriteArgs,
        mode: BamMode,
        decoder: Option<&CramContainerDecoder>,
    ) -> Result<Self> {
        let calls = RastairBedReader::new(params.calls_file.path())
            .wrap_err("Failed to open calls file in worker thread")?;
        let fasta = match mode {
            BamMode::Standard => None,
            BamMode::Legacy => Some(
                open_fasta(params.segments.fasta_file.path())
                    .wrap_err("Failed to open FASTA in worker thread")?,
            ),
        };
        let decoder = decoder
            .map(CramContainerDecoder::fork)
            .transpose()
            .wrap_err("Failed to open the CRAM decoder in worker thread")?;
        Ok(Self { calls, fasta, decoder })
    }

    #[instrument(level = "debug", skip_all, fields(%batch))]
    fn rewrite(
        &mut self,
        header: &BamHeader,
        batch: Batch,
        mode: BamMode,
    ) -> Result<Vec<RawBamRecord>> {
        let mut records = match batch.chunk {
            Chunk::Records(records) => records,
            Chunk::Container(container) => {
                let decoder = self.decoder.as_mut().wrap_err("no CRAM decoder").this_is_a_bug()?;
                let mut records = Vec::new();
                decoder.decode(&container, &mut records).wrap_err("Failed to decode CRAM")?;
                records
            }
        };
        // Runs of records whose span one lookup covers; an unplaced record
        // joins any run.
        let mut start = 0;
        let mut span: Option<Span> = None;
        for end in 0..=records.len() {
            let next = records.get(end).map(Span::of);
            let joined = match next {
                Some(None) => continue,
                Some(Some(this)) => span.map_or(Some(this), |span| span.joined(this)),
                None => None,
            };
            if let Some(joined) = joined {
                span = Some(joined);
                continue;
            }
            let run = records.get_mut(start..end).wrap_err("run out of range").this_is_a_bug()?;
            self.rewrite_run(header, run, span, mode)?;
            start = end;
            span = next.flatten();
        }
        Ok(records)
    }

    /// Rewrite `records`, which all lie in `span` or are unplaced.
    fn rewrite_run(
        &mut self,
        header: &BamHeader,
        records: &mut [RawBamRecord],
        span: Option<Span>,
        mode: BamMode,
    ) -> Result<()> {
        let Some(span) = span else {
            for record in records {
                rewrite_named(&FxHashMap::default(), record, mode, |_| None)?;
            }
            return Ok(());
        };
        let contig = header.target_name(span.tid).wrap_err("record of an unknown contig")?;
        let region = Region { contig: contig.into(), start: span.start, end: span.end };

        let calls: FxHashMap<u32, RastairCall> = self
            .calls
            .query(&region.clone().try_into()?)
            .wrap_err("failed to query calls file")?
            .into_iter()
            .map(|call| (call.pos, call.call))
            .collect();

        // Only a de-novo CpG's context needs the reference, 2 bp to either side.
        let pad = 2u64;
        let ref_start = span.start.saturating_sub(pad);
        let ref_seq = match &self.fasta {
            Some(fasta)
                if calls.values().any(|call| matches!(call, RastairCall::DeNovoCpg { .. })) =>
            {
                let length = header.target_len(span.tid).unwrap_or(span.end);
                let ref_end = span.end.saturating_add(pad).min(length);
                fasta
                    .fetch_seq(&region.contig, ref_start, ref_end)
                    .wrap_err("Failed to fetch FASTA region for context lookup")?
            }
            _ => Vec::new(),
        };
        let ref_base = |pos: u32| -> Option<Base> {
            let idx = usize::try_from(u64::from(pos).checked_sub(ref_start)?).ok()?;
            ref_seq.get(idx).map(|&b| Base::from(b))
        };

        for record in records {
            rewrite_named(&calls, record, mode, ref_base)?;
        }
        Ok(())
    }
}

/// [`rewrite_record`], failing with the read's name.
fn rewrite_named(
    calls: &FxHashMap<u32, RastairCall>,
    record: &mut RawBamRecord,
    mode: BamMode,
    ref_base: impl Fn(u32) -> Option<Base>,
) -> Result<()> {
    rewrite_record(calls, record, mode, ref_base).wrap_err_with(|| {
        format!("failed to rewrite record {}", String::from_utf8_lossy(record.qname()))
    })
}

/// The bases of `record` aligned to the reference, as read index and
/// reference position.
fn aligned_bases(record: &RawBamRecord) -> impl Iterator<Item = MatchPosition> + '_ {
    // An unplaced record has no aligned bases, whatever its CIGAR says (the
    // error `aligned_pairs` gives when it has one).
    record.aligned_pairs().ok().into_iter().flat_map(|pairs| pairs.matches_only())
}

#[instrument(level = "debug", skip_all, fields(pos = ?record.pos()))]
fn rewrite_record(
    calls: &FxHashMap<u32, RastairCall>,
    record: &mut RawBamRecord,
    mode: BamMode,
    ref_base: impl Fn(u32) -> Option<Base>,
) -> Result<()> {
    let flags = record.flags();
    let is_first_in_pair = flags.is_first_in_template();
    let Some(strand) = Strand::from_flags(flags) else {
        // Methylation evidence is a strand's converted base, so a read of
        // unknown strand carries none: no MM/ML, and an XM of all `.`.
        return match mode {
            BamMode::Standard => remove_modifications(record),
            BamMode::Legacy => {
                XrTags::unknown_strand(record.seq_len(), is_first_in_pair).apply_to_record(record)
            }
        };
    };

    match mode {
        BamMode::Standard => {
            // Methylation detection works in stored (+ strand) orientation:
            // T→C for OT, A→G for OB. Positions are stored-SEQ indices.
            let MethylatedInfo { seq, methylated_positions } =
                get_methylated_positions(calls, record, strand);

            // The MM tag spec requires positions relative to the original read
            // (5' to 3'). For forward reads this equals the stored SEQ. For
            // reverse reads the stored SEQ is the reverse complement of the
            // original read, so we must convert both the sequence and positions
            // to original-read orientation, and flip the base/strand qualifier.
            let mods = if flags.is_reverse() {
                let seq_len = seq.len();
                let original_positions: SmallVec<u32, 10> = methylated_positions
                    .iter()
                    .map(|&p| {
                        let flipped = seq_len - 1 - p.as_usize();
                        u32::try_from(flipped).wrap_err("flipped position does not fit in u32")
                    })
                    .collect::<Result<_>>()?;
                let flipped_strand = match strand {
                    Strand::OT => Strand::OB,
                    Strand::OB => Strand::OT,
                };
                MethylatedPositions::new(
                    flipped_strand,
                    &reverse_complement(&seq),
                    &original_positions,
                )
            } else {
                let positions: SmallVec<u32, 10> =
                    methylated_positions.iter().map(|p| p.get()).collect();
                MethylatedPositions::new(strand, &seq, &positions)
            };

            let unmod_base = CpgSide::converted_on(strand).unmod_base();
            for &qpos in &methylated_positions {
                set_base(record, qpos, unmod_base)?;
            }
            mods.apply_to_record(record)?;
        }
        BamMode::Legacy => {
            let annotations = build_legacy_annotations(calls, record, strand, ref_base);
            let xr_tags =
                XrTags::new_legacy(record.seq_len(), strand, is_first_in_pair, &annotations);
            xr_tags.apply_to_record(record)?;
        }
    }

    Ok(())
}

/// Set the base at `qpos`, which is only ever an unmodified C or G.
fn set_base(record: &mut RawBamRecord, qpos: QPos, base: Base) -> Result<()> {
    record.set_base(qpos, IupacBase::from(base)).wrap_err("failed to rewrite a methylated base")
}

/// Build per-read-position XM annotations for legacy mode.
///
/// For each position in the read that overlaps a CpG call (ref or de-novo),
/// determines whether it's methylated and what the reference context is
/// (CpG/CHG/CHH) so the XM tag uses the correct letter.
fn build_legacy_annotations(
    calls: &FxHashMap<u32, RastairCall>,
    record: &RawBamRecord,
    strand: Strand,
    ref_base: impl Fn(u32) -> Option<Base>,
) -> FxHashMap<usize, XmAnnotation> {
    let mut annotations = FxHashMap::default();
    let side = CpgSide::converted_on(strand);

    for MatchPosition { qpos, rpos, .. } in aligned_bases(record) {
        let pos_in_ref = rpos.as_u32();
        let is_denovo = match calls.get(&pos_in_ref) {
            Some(RastairCall::Cpg { .. }) => false,
            Some(RastairCall::DeNovoCpg { .. }) => true,
            _ => continue,
        };

        let Some(base) = record.base_at(qpos).map(IupacBase::to_base) else { continue };
        if let Some(methylated) = side.read_methylated(base) {
            let context = if is_denovo {
                determine_context(pos_in_ref, strand, &ref_base)
            } else {
                MethylationContext::CpG
            };
            annotations.insert(qpos.as_usize(), XmAnnotation { methylated, context });
        }
    }

    annotations
}

struct MethylatedInfo {
    /// The stored SEQ with the methylated bases unmodified, as ASCII
    seq: Vec<u8>,
    methylated_positions: SmallVec<QPos, 10>,
}

fn get_methylated_positions(
    calls: &FxHashMap<u32, RastairCall>,
    record: &RawBamRecord,
    strand: Strand,
) -> MethylatedInfo {
    // The stored SEQ is always in + strand (reference) orientation, regardless
    // of whether the read mapped to the reverse strand. Methylation evidence
    // bases (T for OT at C, A for OB at G) are also defined in + strand terms,
    // and the aligned pairs index the stored SEQ. So we work entirely in stored
    // orientation — no reverse-complement needed.
    let mut seq = Vec::with_capacity(record.seq_len());
    record.decode_seq_into(&mut seq);
    let mut methylated_positions = SmallVec::new();
    let side = CpgSide::converted_on(strand);

    for MatchPosition { qpos, rpos, .. } in aligned_bases(record) {
        match calls.get(&rpos.as_u32()) {
            Some(RastairCall::Cpg { .. } | RastairCall::DeNovoCpg { .. }) => {}
            _ => continue,
        }

        let Some(base) = seq.get_mut(qpos.as_usize()) else { continue };
        if side.read_methylated(Base::from(*base)) == Some(true) {
            *base = *side.unmod_base();
            methylated_positions.push(qpos);
        }
    }

    MethylatedInfo { seq, methylated_positions }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&base| match base {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect()
}

/// The `@PG` line for rastair; `add_pg` numbers its ID if the input already
/// has one (`rastair.rewrite_bam.1`, …), as on a rerun on rastair's output.
fn rastair_pg() -> PgRecord {
    PgRecord {
        id: "rastair.rewrite_bam".into(),
        pn: Some("rastair".into()),
        vn: Some(env!("CARGO_PKG_VERSION").into()),
        // A tab or newline in an argument would break the header line.
        cl: Some(
            std::env::args_os()
                .skip(1)
                .map(|arg| arg.to_string_lossy().into_owned())
                .collect::<Vec<_>>()
                .join(" ")
                .replace(char::is_control, " ")
                .into(),
        ),
        ds: Some("Rewrote BAM with methylation information".into()),
        pp: None,
    }
}

#[cfg(test)]
#[allow(clippy::cast_possible_truncation, reason = "lots of noise otherwise for small numbers")]
mod tests;
