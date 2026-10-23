//! `rastair bam`: rewrite a BAM with per-read methylation tags.
//!
//! Two modes: `legacy` writes XR/XG/XM and leaves SEQ alone; `standard` writes
//! MM/ML and rewrites SEQ (T→C on OT, A→G on OB). Tags are built in
//! [`base_modification`].
//!
//! **XM and MM/ML encode per-read methylation, not position-level calls.** A
//! CpG with beta 0.3 (`methylated: false` in [`RastairCall::Cpg`]) still has
//! methylated reads that must show `Z` in XM and appear in MM/ML; `methylated`
//! only decides whether the position is a CpG site. Per-read methylation
//! comes from the observed base alone: T at an OT C, A at an OB G.
//!
//! **MM/ML and XM count different reads.** MM/ML only encode modifications at
//! C bases of the stored SEQ, and of two mates overlapping a CpG only one has
//! a C there (the other has the G), while XM annotates both. So MM/ML readers
//! such as modkit see about half the reads XM counting does: fractions agree,
//! counts differ ~2:1. Compare the two modes by fraction, at a minimum
//! coverage.

use crate::{
    bed::reader::{RastairBedReader, RastairCall},
    sequence::{ChunkRegion, ReaderParams},
    utils::{
        cli,
        file_helpers::{FastaReader, open_fasta},
        logging::ThisIsABug,
    },
};
use clap::{Parser, value_parser};
use clio::ClioPath;
use color_eyre::eyre::{Context, Result};
use rust_htslib::bam::{
    self, FetchDefinition, Header, Read, Record, Writer, ext::BamRecordExtensions as _,
    header::HeaderRecord,
};
use rustc_hash::FxHashMap;
use seqair_types::SmallVec;
use seqair_types::{Base, Strand, StrandFromRecord};
use std::thread::available_parallelism;

mod base_modification;
use crate::runtime::{
    fault_injection, partial_output::PartialOutput, progress, progress::ProgressTracker, segments,
    threads,
};
pub use base_modification::{
    MethylatedPositions, MethylationContext, XmAnnotation, XrTags, determine_context,
};
use tracing::{instrument, warn};

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
    /// usage.
    #[arg(long, default_value_t = 10_000)]
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

#[tracing::instrument(level = "info", skip_all, fields(
    ?mode,
))]
pub fn rewrite(params: &BamRewriteArgs, mode: BamMode) -> Result<()> {
    progress::register_signal_handler();
    let panics = threads::PanicCheck::start();

    // Open BAM on the main thread to get header and compute regions
    let (header, regions) = {
        let readers = params.segments.readers().wrap_err("Failed to read BAM/FASTA files")?;
        let regions: Vec<ChunkRegion> = readers
            .segments(params.segment_max_length, 0)
            .wrap_err("Could not fetch segments from BAM file")?
            .collect();
        let header = {
            let mut h = Header::from_template(readers.bam.header());
            add_rastair_header(&mut h);
            h
        };
        (header, regions)
    };

    let worker_threads = params.threads.saturating_sub(1).max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .thread_name(|idx| format!("bam-worker-{idx}"))
        .num_threads(worker_threads)
        .build()
        .wrap_err("Failed to create thread pool for BAM rewrite")?;

    let output = PartialOutput::new("BAM output", &params.output);
    let output_file = output.write_path()?;
    let mut writer = if output_file.is_std() {
        Writer::from_stdout(&header, bam::Format::Bam)
    } else {
        Writer::from_path(output_file.path(), &header, bam::Format::Bam)
    }
    .wrap_err("failed to create writer")?;
    writer
        .set_compression_level(bam::CompressionLevel::Fastest)
        .wrap_err("failed to set compression level")?;
    writer.set_threads(3).wrap_err("failed to set threads")?;

    let mut progress = ProgressTracker::new(regions.len());
    let last = regions.last();
    let processed = segments::process_in_order(
        &pool,
        &regions,
        || WorkerReaders::open(params),
        |WorkerReaders { bam, calls, fasta }, segment| {
            fault_injection::fault_point(fault_injection::FaultPoint::Worker)?;
            let is_last = last.is_some_and(|last| std::ptr::eq(segment, last));
            rewrite_region(bam, calls, fasta, segment, mode, is_last)
        },
        |records| {
            for record in records {
                writer.write(&record).wrap_err("failed to write record to new BAM file")?;
            }
            progress.segment_done();
            fault_injection::fault_point(fault_injection::FaultPoint::Writer)
        },
    )
    .wrap_err("Failed to process BAM regions in parallel");

    // Close the output even if processing failed, so it holds a valid prefix
    drop(writer);
    let closed = if output_file.is_std() { Ok(()) } else { ensure_closed(output_file.path()) };
    segments::first_error(processed, closed)?;
    panics.finish()?;
    output.commit()?;

    Ok(())
}

/// What each worker reads from, opened once on its own thread.
struct WorkerReaders {
    bam: bam::IndexedReader,
    calls: RastairBedReader,
    fasta: FastaReader,
}

impl WorkerReaders {
    fn open(params: &BamRewriteArgs) -> Result<Self> {
        let mut bam = bam::IndexedReader::from_path(params.segments.bam_file.path())
            .wrap_err("Failed to open BAM/CRAM in worker thread")?;
        bam.set_reference(params.segments.fasta_file.path())
            .wrap_err("Failed to set FASTA reference for BAM/CRAM reader in worker thread")?;
        let calls = RastairBedReader::new(params.calls_file.path())
            .wrap_err("Failed to open calls file in worker thread")?;
        let fasta = open_fasta(params.segments.fasta_file.path())
            .wrap_err("Failed to open FASTA in worker thread")?;
        Ok(Self { bam, calls, fasta })
    }
}

#[instrument(level = "debug", skip_all, fields(region = %region.region))]
fn rewrite_region(
    bam: &mut bam::IndexedReader,
    calls_reader: &mut RastairBedReader,
    fasta: &mut FastaReader,
    region: &ChunkRegion,
    mode: BamMode,
    is_last_segment: bool,
) -> Result<Vec<Record>> {
    FetchDefinition::try_from(&region.region)
        .wrap_err("Could not convert region string")
        .and_then(|r| bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))
        .wrap_err_with(|| format!("Could not fetch region `{}` from BAM file", region.region))?;

    let noodle_region = region
        .region
        .clone()
        .try_into()
        .wrap_err("Failed to convert region representation")
        .this_is_a_bug()?;

    let calls: FxHashMap<u32, RastairCall> = calls_reader
        .query(&noodle_region)
        .wrap_err("failed to query calls file")?
        .iter()
        .map(|call| (call.pos, call.call.clone()))
        .collect();

    // Fetch the reference sequence for this region (with 2bp padding for context lookup).
    // Reads may extend beyond the region, so lookups that fall outside return None.
    let pad = 2u64;
    let ref_start = region.start.saturating_sub(pad);
    let ref_end = region.end.saturating_add(pad).min(region.last_position);
    let ref_seq = fasta
        .fetch_seq(&region.contig, ref_start, ref_end)
        .wrap_err("Failed to fetch FASTA region for context lookup")?;
    let ref_base = |pos: u32| -> Option<Base> {
        let idx = u64::from(pos).checked_sub(ref_start)? as usize;
        ref_seq.get(idx).map(|&b| Base::from(b))
    };

    let region_start = region.start.cast_signed();
    let region_end = region.end.cast_signed();
    let mut record = Record::new();
    let mut out = Vec::new();
    while let Some(result) = bam.read(&mut record) {
        if let Err(error) = result {
            warn!(%error, "Failed to read BAM record");
            continue;
        }

        // `bam.fetch` returns all reads overlapping the region. Adjacent
        // segments share boundary positions (segment N ends at X, segment N+1
        // starts at X), so reads at the boundary would be emitted twice.
        // Each segment owns reads in [region_start, region_end); the last
        // segment also includes region_end.
        let pos = record.pos();
        if pos < region_start || (!is_last_segment && pos >= region_end) {
            continue;
        }

        rewrite_record(&calls, &mut record, mode, ref_base).wrap_err("failed to rewrite record")?;
        out.push(record.clone());
    }

    Ok(out)
}

#[instrument(level = "debug", skip_all, fields(pos = record.pos()))]
fn rewrite_record(
    calls: &FxHashMap<u32, RastairCall>,
    record: &mut Record,
    mode: BamMode,
    ref_base: impl Fn(u32) -> Option<Base>,
) -> Result<()> {
    let strand = StrandFromRecord::strand(record);
    let is_first_in_pair = record.is_first_in_template();

    let is_reverse = record.is_reverse();

    match mode {
        BamMode::Standard => {
            // Methylation detection works in stored (+ strand) orientation:
            // T→C for OT, A→G for OB. Positions are stored-SEQ indices.
            let MethylatedInfo { seq, methylated_positions } =
                get_methylated_positions(calls, record)?;

            // The MM tag spec requires positions relative to the original read
            // (5' to 3'). For forward reads this equals the stored SEQ. For
            // reverse reads the stored SEQ is the reverse complement of the
            // original read, so we must convert both the sequence and positions
            // to original-read orientation, and flip the base/strand qualifier.
            let (mm_seq, mm_positions, mm_strand) = if is_reverse {
                let seq_len = seq.len();
                let original_positions: SmallVec<u32, 10> = methylated_positions
                    .iter()
                    .map(|&p| {
                        let flipped = seq_len - 1 - p as usize;
                        u32::try_from(flipped).wrap_err("flipped position does not fit in u32")
                    })
                    .collect::<Result<_>>()?;
                let original_seq = reverse_complement(&seq);
                let flipped_strand = match strand {
                    Strand::OT => Strand::OB,
                    Strand::OB => Strand::OT,
                    Strand::Unknown => Strand::Unknown,
                };
                (original_seq, original_positions, flipped_strand)
            } else {
                (seq.clone(), methylated_positions, strand)
            };

            let mods = MethylatedPositions::new(mm_strand, &mm_seq, &mm_positions);
            record.set_seq(&seq);
            mods.apply_to_record(record)?;
        }
        BamMode::Legacy => {
            let annotations = build_legacy_annotations(calls, record, ref_base)?;
            let xr_tags =
                XrTags::new_legacy(record.seq_len(), strand, is_first_in_pair, &annotations);
            xr_tags.apply_to_record(record)?;
        }
    }

    Ok(())
}

/// Build per-read-position XM annotations for legacy mode.
///
/// For each position in the read that overlaps a CpG call (ref or de-novo),
/// determines whether it's methylated and what the reference context is
/// (CpG/CHG/CHH) so the XM tag uses the correct letter.
fn build_legacy_annotations(
    calls: &FxHashMap<u32, RastairCall>,
    record: &Record,
    ref_base: impl Fn(u32) -> Option<Base>,
) -> Result<FxHashMap<usize, XmAnnotation>> {
    let strand = StrandFromRecord::strand(record);
    let seq = record.seq().as_bytes();
    let mut annotations = FxHashMap::default();

    let (target_base, evidence_base) = match strand {
        Strand::OT => (Base::C, Base::T),
        Strand::OB => (Base::G, Base::A),
        Strand::Unknown => return Ok(annotations),
    };

    for [pos_in_read, pos_in_ref] in record.aligned_pairs_full() {
        let Some(pos_in_read) = pos_in_read else { continue };
        let Some(pos_in_ref) = pos_in_ref else { continue };
        let pos_in_ref =
            u32::try_from(pos_in_ref).wrap_err("reference position does not fit in u32")?;
        let pos_in_read =
            usize::try_from(pos_in_read).wrap_err("read position does not fit in usize")?;

        let is_denovo = match calls.get(&pos_in_ref) {
            Some(RastairCall::Cpg { .. }) => false,
            Some(RastairCall::DeNovoCpg { .. }) => true,
            _ => continue,
        };

        let Some(&raw_base) = seq.get(pos_in_read) else { continue };
        let observed_base = Base::from(raw_base);

        // A CpG position in the read shows either the target base (unmethylated)
        // or the evidence base (methylated, converted by TAPS)
        if observed_base == target_base || observed_base == evidence_base {
            let context = if is_denovo {
                determine_context(pos_in_ref, strand, &ref_base)
            } else {
                MethylationContext::CpG
            };
            let methylated = observed_base == evidence_base;
            annotations.insert(pos_in_read, XmAnnotation { methylated, context });
        }
    }

    Ok(annotations)
}

struct MethylatedInfo {
    seq: Vec<u8>,
    methylated_positions: SmallVec<u32, 10>,
}

fn get_methylated_positions(
    calls: &FxHashMap<u32, RastairCall>,
    record: &Record,
) -> Result<MethylatedInfo> {
    use Base::*;

    let strand = StrandFromRecord::strand(record);

    // The stored SEQ is always in + strand (reference) orientation, regardless
    // of whether the read mapped to the reverse strand. Methylation evidence
    // bases (T for OT at C, A for OB at G) are also defined in + strand terms,
    // and `aligned_pairs_full` returns positions into the stored SEQ. So we work
    // entirely in stored orientation — no reverse-complement needed.
    let mut seq = record.seq().as_bytes();
    let mut methylated_positions: SmallVec<u32, 10> = SmallVec::new();

    for [pos_in_read, pos_in_ref] in record.aligned_pairs_full() {
        let Some(pos_in_read) = pos_in_read else {
            continue;
        };
        let Some(pos_in_ref) = pos_in_ref else {
            continue;
        };
        let pos_in_ref =
            u32::try_from(pos_in_ref).wrap_err("reference position does not fit in u32")?;
        let pos_in_read =
            usize::try_from(pos_in_read).wrap_err("read position does not fit in usize")?;

        // Only process positions that are called as CpG sites (ref or de-novo).
        // Per-read methylation is determined by the observed base, not the
        // position-level call: in TAPS, T at OT C = methylated, A at OB G = methylated.
        match calls.get(&pos_in_ref) {
            Some(RastairCall::Cpg { .. } | RastairCall::DeNovoCpg { .. }) => {}
            _ => continue,
        };

        {
            let Some(&raw_base) = seq.get(pos_in_read) else { continue };
            let observed_base = Base::from(raw_base);

            match strand {
                Strand::OT => {
                    if observed_base == T {
                        seq[pos_in_read] = *C;
                        methylated_positions.push(
                            u32::try_from(pos_in_read)
                                .wrap_err("read position does not fit in u32")?,
                        );
                    }
                }
                Strand::OB => {
                    if observed_base == A {
                        seq[pos_in_read] = *G;
                        methylated_positions.push(
                            u32::try_from(pos_in_read)
                                .wrap_err("read position does not fit in u32")?,
                        );
                    }
                }
                Strand::Unknown => continue,
            }
        }
    }

    Ok(MethylatedInfo { seq, methylated_positions })
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

/// Add a PG header line for rastair
fn add_rastair_header(header: &mut Header) {
    header.push_record(
        HeaderRecord::new(b"PG")
            .push_tag(b"ID", "rastair.rewrite_bam")
            .push_tag(b"PN", "rastair")
            .push_tag(b"VN", env!("CARGO_PKG_VERSION"))
            .push_tag(b"CL", std::env::args().skip(1).collect::<Vec<_>>().join(" "))
            .push_tag(b"DS", "Rewrote BAM with methylation information"),
    );
}

/// rust-htslib closes a BAM writer in `Drop` and ignores whether that worked,
/// and with writer threads, a full disk may only show then. htslib writes the
/// BGZF EOF block last, after flushing everything else, so check for it.
// ponytail: misses a failure after the EOF block (the final close(2)); use an
// explicit close once rust-htslib has one.
fn ensure_closed(path: &std::path::Path) -> Result<()> {
    use std::io::{Read as _, Seek as _, SeekFrom};
    const BGZF_EOF: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];
    let mut end = [0; 28];
    let read = std::fs::File::open(path).and_then(|mut file| {
        file.seek(SeekFrom::End(-28))?;
        file.read_exact(&mut end)
    });
    color_eyre::eyre::ensure!(
        read.is_ok() && end == BGZF_EOF,
        "Failed to write the BAM file `{}`: it does not end in a BGZF EOF block, the disk may be full",
        path.display()
    );
    Ok(())
}

#[cfg(test)]
#[allow(clippy::cast_possible_truncation, reason = "lots of noise otherwise for small numbers")]
mod tests;
