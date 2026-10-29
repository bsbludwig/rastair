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
//! **The input is streamed front to back**, so every record comes out, in
//! input order, unplaced ones included, and no index is needed. The reading
//! side (`RecordBatches`) cuts the stream into `Batch`es of up to
//! `BATCH_RECORDS` records on one contig within `BATCH_MAX_SPAN` bases
//! (unplaced records join any batch). The workers pull batches one at a time
//! ([`segments::process_in_order`]), so reading is serial and runs on whichever
//! worker asks next; a worker looks up the calls (and, for de-novo CpGs in
//! legacy mode, the reference) over its batch's span once, rewrites the
//! records in place, and the writer writes the batches in order. The span cap
//! keeps that lookup small for unsorted input, which only makes batches
//! smaller.
//!
//! `--region` keeps the records overlapping it, each once; it filters the
//! stream, so the whole file is still read. `RecordBatches` is the one place
//! that knows where records come from: an indexed region query would be
//! another source of `Batch`es, and so would CRAM, where a batch would be one
//! container, decoded by the worker that takes it rather than by the reader.
//!
//! A record of unknown strand (paired, but neither first nor second) carries
//! no methylation; an unmapped or unplaced one has no aligned bases, so it gets
//! the tags of a read without methylation evidence.

use crate::{
    bed::reader::{RastairBedReader, RastairCall},
    sequence::{ReaderParams, Region},
    utils::{
        CliRegionInput, cli,
        file_helpers::{FastaReader, open_fasta},
    },
    vcf::CpgSide,
};
use clap::{Parser, value_parser};
use clio::ClioPath;
use color_eyre::eyre::{Context, ContextCompat as _, Report, Result, ensure};
use rustc_hash::FxHashMap;
use seqair::{
    bam::{
        BamHeader, CigarOp, RawBamRecord,
        aligned_pairs::{AlignedPairs, MatchPosition},
        header::PgRecord,
        writer::BamWriterBuilder,
    },
    reader::{AlignmentStream, Tid},
};
use seqair_types::{Base, QPos, SmallVec, Strand};
use std::{
    fmt,
    io::{self, Write as _},
    num::NonZeroU8,
    sync::{
        Arc,
        atomic::{AtomicU64, Ordering},
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
    /// Maximum length of a segment in bases
    ///
    /// Used for splitting work between threads. Tweak this to adjust memory
    /// usage. Between 1,000 and 1,000,000.
    #[arg(long, default_value_t = 10_000, value_parser = crate::sequence::segment_length())]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub segment_max_length: u64,

    /// Number of threads to use for processing the BAM file.
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
}

/// Records per [`Batch`]: about what a 10 kb segment held at 30x.
const BATCH_RECORDS: usize = 2_000;
/// The longest reference span of a [`Batch`], which bounds its calls and
/// reference lookup.
const BATCH_MAX_SPAN: u64 = 1_000_000;
/// Threads inflating the input's BGZF blocks.
const INFLATE_THREADS: NonZeroU8 = NonZeroU8::MIN.saturating_add(1);
/// Threads compressing the output, as for `samtools view -@ 3`.
const COMPRESS_THREADS: usize = 3;

#[tracing::instrument(level = "info", skip_all, fields(
    ?mode,
))]
pub fn rewrite(params: &BamRewriteArgs, mode: BamMode) -> Result<()> {
    progress::register_signal_handler();
    let panics = threads::PanicCheck::start();

    let batches = RecordBatches::open(&params.segments)?;
    let input_bytes = batches.input_bytes;
    let input_header = BamHeader::from_template(batches.stream.header());
    let mut header = BamHeader::from_template(batches.stream.header());
    header.add_pg(rastair_pg());

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
        .compression_level(1)
        .compression_threads(COMPRESS_THREADS)
        .build()
        .wrap_err("failed to create writer")?;

    let kib = |bytes: u64| usize::try_from(bytes / 1024).unwrap_or(usize::MAX);
    let mut progress = ProgressTracker::counting(kib(input_bytes), "KiB of input");
    let processed = segments::process_in_order(
        &pool,
        batches.map(Pulled),
        || Worker::open(params, mode),
        |worker, Pulled(batch)| {
            fault_injection::fault_point(fault_injection::FaultPoint::Worker)?;
            let batch = batch?;
            let read = batch.read;
            worker.rewrite(&input_header, batch, mode).map(|records| (records, read))
        },
        |(records, read)| {
            for record in &records {
                writer.write_raw(record).wrap_err("failed to write record to new BAM file")?;
            }
            progress.done_up_to(kib(read));
            fault_injection::fault_point(fault_injection::FaultPoint::Writer)
        },
    )
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

/// Consecutive records of the input, the unit of work.
struct Batch {
    /// The input index of the first record, from 0
    first: u64,
    records: Vec<RawBamRecord>,
    /// Where the placed records lie, if any
    span: Option<Span>,
    /// How many input bytes were read once the batch was complete
    read: u64,
}

impl fmt::Display for Batch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let last = self.first.saturating_add(self.records.len() as u64);
        write!(f, "records {}-{last} of the input", self.first.saturating_add(1))
    }
}

/// What a worker pulls from [`RecordBatches`]: a batch, or why there is none.
struct Pulled(Result<Batch>);

impl fmt::Display for Pulled {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.0 {
            Ok(batch) => batch.fmt(f),
            Err(_) => f.write_str("the next records of the input"),
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

    /// The span covering both, if they share a contig and it stays short.
    fn joined(self, other: Self) -> Option<Self> {
        let joined = Self {
            tid: self.tid,
            start: self.start.min(other.start),
            end: self.end.max(other.end),
        };
        (self.tid == other.tid && joined.end - joined.start <= BATCH_MAX_SPAN).then_some(joined)
    }
}

/// Reads the input front to back as [`Batch`]es: where the records come from.
struct RecordBatches {
    stream: AlignmentStream,
    regions: Option<RegionFilter>,
    /// The size of the input file
    input_bytes: u64,
    read: Arc<AtomicU64>,
    /// The record that started the next batch
    pending: Option<RawBamRecord>,
    /// A read error, returned after the batch of the records before it
    error: Option<Report>,
    /// The input index of the next record
    next: u64,
}

impl RecordBatches {
    fn open(params: &ReaderParams) -> Result<Self> {
        let path = params.bam_file.path();
        let file = std::fs::File::open(path)
            .wrap_err_with(|| format!("Failed to open `{}`", path.display()))?;
        let input_bytes =
            file.metadata().wrap_err_with(|| format!("Failed to read `{}`", path.display()))?.len();
        let read = Arc::new(AtomicU64::new(0));
        let inflate_pool = rayon::ThreadPoolBuilder::new()
            .thread_name(|idx| format!("bam-inflate-{idx}"))
            .num_threads(std::num::NonZeroUsize::from(INFLATE_THREADS).get())
            .build()
            .wrap_err("Failed to create thread pool for reading the BAM file")?;
        let stream =
            AlignmentStream::builder(CountingReader { inner: file, read: Arc::clone(&read) })
                .threads_in(INFLATE_THREADS, inflate_pool)
                .build()
                .wrap_err_with(|| format!("Failed to read `{}`", path.display()))?;
        let regions = params
            .regions
            .as_ref()
            .map(|regions| RegionFilter::new(regions, stream.header()))
            .transpose()?;
        Ok(Self { stream, regions, input_bytes, read, pending: None, error: None, next: 0 })
    }
}

impl Iterator for RecordBatches {
    type Item = Result<Batch>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(error) = self.error.take() {
            return Some(Err(error));
        }
        let first = self.next;
        let mut records = Vec::with_capacity(BATCH_RECORDS);
        let mut span = None;
        let mut next = self.pending.take();
        while records.len() < BATCH_RECORDS {
            let record = match next.take() {
                Some(record) => record,
                None => {
                    let mut record = RawBamRecord::default();
                    match self.stream.read_record(&mut record) {
                        Ok(true) => {}
                        Ok(false) => break,
                        Err(error) => {
                            self.error = Some(
                                Report::new(error).wrap_err("Failed to read the input BAM file"),
                            );
                            break;
                        }
                    }
                    if self.regions.as_ref().is_some_and(|regions| !regions.keeps(&record)) {
                        continue;
                    }
                    record
                }
            };
            if let Some(this) = Span::of(&record) {
                match span.map_or(Some(this), |span: Span| span.joined(this)) {
                    Some(joined) => span = Some(joined),
                    None => {
                        self.pending = Some(record);
                        break;
                    }
                }
            }
            records.push(record);
        }
        if records.is_empty() {
            return self.error.take().map(Err);
        }
        self.next = first.saturating_add(records.len() as u64);
        Some(Ok(Batch { first, records, span, read: self.read.load(Ordering::Relaxed) }))
    }
}

/// `--region`, as sorted, non-overlapping spans.
struct RegionFilter(Vec<Span>);

impl RegionFilter {
    fn new(input: &CliRegionInput, header: &BamHeader) -> Result<Self> {
        let mut spans = input
            .regions()
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

    /// Whether `record` overlaps a region; an unplaced one does not.
    fn keeps(&self, record: &RawBamRecord) -> bool {
        let Some(span) = Span::of(record) else { return false };
        // The regions do not overlap, so of those starting before the record
        // ends, only the last one can reach it.
        let before =
            self.0.partition_point(|region| (region.tid, region.start) < (span.tid, span.end));
        before
            .checked_sub(1)
            .and_then(|index| self.0.get(index))
            .is_some_and(|region| region.tid == span.tid && region.end > span.start)
    }
}

/// What each worker reads from, opened once on its own thread.
struct Worker {
    calls: RastairBedReader,
    /// Only legacy mode needs the reference
    fasta: Option<FastaReader>,
}

impl Worker {
    fn open(params: &BamRewriteArgs, mode: BamMode) -> Result<Self> {
        let calls = RastairBedReader::new(params.calls_file.path())
            .wrap_err("Failed to open calls file in worker thread")?;
        let fasta = match mode {
            BamMode::Standard => None,
            BamMode::Legacy => Some(
                open_fasta(params.segments.fasta_file.path())
                    .wrap_err("Failed to open FASTA in worker thread")?,
            ),
        };
        Ok(Self { calls, fasta })
    }

    #[instrument(level = "debug", skip_all, fields(%batch))]
    fn rewrite(
        &mut self,
        header: &BamHeader,
        batch: Batch,
        mode: BamMode,
    ) -> Result<Vec<RawBamRecord>> {
        let Batch { mut records, span, .. } = batch;
        let Some(span) = span else {
            for record in &mut records {
                rewrite_record(&FxHashMap::default(), record, mode, |_| None)?;
            }
            return Ok(records);
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

        for record in &mut records {
            rewrite_record(&calls, record, mode, ref_base).wrap_err_with(|| {
                format!("failed to rewrite record {}", String::from_utf8_lossy(record.qname()))
            })?;
        }
        Ok(records)
    }
}

/// The bases of `record` aligned to the reference, as read index and
/// reference position; none for an unplaced record.
fn aligned_bases(record: &RawBamRecord) -> Vec<MatchPosition> {
    let Some(pos) = record.pos() else { return Vec::new() };
    let cigar: SmallVec<CigarOp, 8> = record.cigar().collect();
    AlignedPairs::new(pos, &cigar).matches_only().collect()
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
    record
        .set_base(qpos, seqair_types::IupacBase::from(base))
        .wrap_err("failed to rewrite a methylated base")
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
    let seq = record.bases();
    let mut annotations = FxHashMap::default();
    let side = CpgSide::converted_on(strand);

    for MatchPosition { qpos, rpos, .. } in aligned_bases(record) {
        let pos_in_ref = rpos.as_u32();
        let is_denovo = match calls.get(&pos_in_ref) {
            Some(RastairCall::Cpg { .. }) => false,
            Some(RastairCall::DeNovoCpg { .. }) => true,
            _ => continue,
        };

        let Some(&base) = seq.get(qpos.as_usize()) else { continue };
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

/// The `@PG` line for rastair
fn rastair_pg() -> PgRecord {
    PgRecord {
        id: "rastair.rewrite_bam".into(),
        pn: Some("rastair".into()),
        vn: Some(env!("CARGO_PKG_VERSION").into()),
        cl: Some(std::env::args().skip(1).collect::<Vec<_>>().join(" ").into()),
        ds: Some("Rewrote BAM with methylation information".into()),
        pp: None,
    }
}

#[cfg(test)]
#[allow(clippy::cast_possible_truncation, reason = "lots of noise otherwise for small numbers")]
mod tests;
