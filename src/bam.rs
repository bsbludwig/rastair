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
use color_eyre::eyre::{Context, ContextCompat, Result, eyre};
use rastair_types::SmallVec;
use rastair_types::{Base, Strand, StrandFromRecord};
use rayon::prelude::*;
use rust_htslib::bam::{
    self, FetchDefinition, Header, Read, Record, Writer, ext::BamRecordExtensions as _,
    header::HeaderRecord,
};
use rustc_hash::FxHashMap;
use std::thread::available_parallelism;

mod base_modification;
use crate::progress::ProgressTracker;
pub use base_modification::{
    MethylatedPositions, MethylationContext, XmAnnotation, XrTags, determine_context,
};
use tracing::{instrument, trace, warn};

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
    #[arg(long, default_value_t = 100_000)]
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
    output = %params.output.path().display(),
    ?mode,
))]
pub fn rewrite(params: &BamRewriteArgs, mode: BamMode) -> Result<()> {
    crate::progress::register_signal_handler();

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
    let (bam_sender, bam_receiver) = ordered_channel::bounded::<Vec<Record>>(worker_threads * 10);

    let total_segments = regions.len();
    let output_file = params.output.clone();
    let writer_thread = std::thread::Builder::new()
        .name("bam-writer".to_string())
        .spawn(move || -> Result<()> {
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

            let mut progress = ProgressTracker::new(total_segments);
            for records in bam_receiver {
                for record in records {
                    writer.write(&record).wrap_err("failed to write record to new BAM file")?;
                }
                progress.segment_done();
            }
            Ok(())
        })
        .wrap_err("failed to spawn BAM writer thread")?;

    let n_regions = regions.len();
    rayon::ThreadPoolBuilder::new()
        .thread_name(|idx| format!("bam-worker-{idx}"))
        .num_threads(worker_threads)
        .build()
        .wrap_err("Failed to create thread pool for BAM rewrite")?
        .install(move || {
            regions.iter().enumerate().par_bridge().try_for_each_with(
                bam_sender,
                |sender, (index, segment)| {
                    let is_last = index == n_regions - 1;
                    rewrite_region_parallel(index, segment, is_last, sender, params, mode)
                },
            )
        })
        .wrap_err("Failed to process BAM regions in parallel")?;

    writer_thread
        .join()
        .map_err(|_| eyre!("BAM writer thread panicked"))
        .this_is_a_bug()?
        .wrap_err("Error in BAM writer thread")?;

    Ok(())
}

/// Wrapper for parallel BAM region processing with thread-local readers.
fn rewrite_region_parallel(
    index: usize,
    segment: &ChunkRegion,
    is_last: bool,
    sender: &mut ordered_channel::Sender<Vec<Record>>,
    params: &BamRewriteArgs,
    mode: BamMode,
) -> Result<()> {
    thread_local! {
        static BAM_READER: std::cell::RefCell<Option<bam::IndexedReader>> =
            const { std::cell::RefCell::new(None) };
        static BED_READER: std::cell::RefCell<Option<RastairBedReader>> =
            const { std::cell::RefCell::new(None) };
        static FASTA_READER: std::cell::RefCell<Option<FastaReader>> =
            const { std::cell::RefCell::new(None) };
    }

    let records = BAM_READER.with(|bam_cell| {
        BED_READER.with(|bed_cell| {
            FASTA_READER.with(|fasta_cell| -> Result<Vec<Record>> {
                let mut bam_opt = bam_cell.borrow_mut();
                let mut bed_opt = bed_cell.borrow_mut();
                let mut fasta_opt = fasta_cell.borrow_mut();

                if bam_opt.is_none() {
                    *bam_opt = Some(
                        bam::IndexedReader::from_path(params.segments.bam_file.path())
                            .wrap_err("Failed to open BAM in worker thread")?,
                    );
                }
                if bed_opt.is_none() {
                    *bed_opt = Some(
                        RastairBedReader::new(params.calls_file.path())
                            .wrap_err("Failed to open calls file in worker thread")?,
                    );
                }
                if fasta_opt.is_none() {
                    *fasta_opt = Some(
                        open_fasta(params.segments.fasta_file.path())
                            .wrap_err("Failed to open FASTA in worker thread")?,
                    );
                }

                let bam = bam_opt
                    .as_mut()
                    .wrap_err("thread-local BAM reader not initialized")
                    .this_is_a_bug()?;
                let bed = bed_opt
                    .as_mut()
                    .wrap_err("thread-local BED reader not initialized")
                    .this_is_a_bug()?;
                let fasta = fasta_opt
                    .as_mut()
                    .wrap_err("thread-local FASTA reader not initialized")
                    .this_is_a_bug()?;

                rewrite_region(bam, bed, fasta, segment, mode, is_last)
            })
        })
    })?;

    if let Err(err) = sender.send(index, records) {
        trace!(error = format!("{err:#}"), "Failed to send BAM records, channel probably closed");
    }

    Ok(())
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
    let ref_seq = {
        let mut seq = Vec::new();
        fasta
            .fetch(&region.contig, ref_start, ref_end)
            .wrap_err("Failed to fetch FASTA region for context lookup")?;
        fasta.read(&mut seq).wrap_err("Failed to read FASTA sequence")?;
        seq
    };
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

#[cfg(test)]
#[allow(clippy::cast_possible_truncation, reason = "lots of noise otherwise for small numbers")]
mod tests {
    use color_eyre::eyre::{ContextCompat as _, bail, ensure};
    use insta::{assert_compact_debug_snapshot, assert_snapshot};
    use rust_htslib::bam::record::Aux;

    use super::*;

    /// No-op reference lookup for tests that don't involve de-novo CpGs
    fn no_ref(_pos: u32) -> Option<Base> {
        None
    }

    /// Golden test: rewrite the same read through both Legacy (XM) and Standard (MM)
    /// modes, verify exact output strings and cross-check that both agree on which
    /// positions are methylated.
    ///
    /// Uses a real read from the test BAM with manually constructed calls that have
    /// a known expected outcome.
    #[test]
    fn golden_legacy_and_standard_agree() -> Result<()> {
        use Base::*;

        // Read 2 in the test region: flag 163, OB strand, NOT reversed, second in pair
        // SEQ: AAGGCATGCACCACCACGCCTGGCTTGGTTTGGTTTTTGATTGGTTGGTTGGTCTTTTGAGACAGGGTTTCTCTGTGT
        // BAM pos: 6103076 (0-based)
        // Ref: aaggcgtgcaccaccacgcctggcttggtttggtttttgattggttgGttggtcttttgagacagggtttctctgtgt
        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
        bam.fetch(FetchDefinition::All)?;
        let mut record = Record::new();
        while let Some(result) = bam.read(&mut record) {
            result?;
            if record.flags() == 163 && record.pos() == 6103076 {
                break;
            }
        }
        ensure!(record.flags() == 163, "could not find flag-163 record at 6103076");

        let strand = StrandFromRecord::strand(&record);
        assert_eq!(strand, Strand::OB, "expected OB strand");
        assert!(!record.is_reverse(), "expected forward read");
        let seq = record.seq().as_bytes();
        assert_snapshot!(as_base_string(&seq), @"AAGGCATGCACCACCACGCCTGGCTTGGTTTGGTTTTTGATTGGTTGGTTGGTCTTTTGAGACAGGGTTTCTCTGTGT");

        // Construct calls for two CpG positions (G side, matching OB strand):
        //
        // Ref pos 6103081 (read offset 5): methylated CpG
        //   Read[5] = A → methylation evidence on OB (A at G ref position)
        //
        // Ref pos 6103093 (read offset 17): unmethylated CpG
        //   Read[17] = G → target base, no methylation evidence
        let calls = FxHashMap::from_iter([
            (6103081, RastairCall::Cpg { base: G, methylated: true }),
            (6103093, RastairCall::Cpg { base: G, methylated: false }),
        ]);

        // Reference lookup for context determination (both positions are real CpGs)
        // Ref around 6103081: ...cgtgc... (pos-1=C → CpG context on bottom strand)
        // Ref around 6103093: ...acgcc... (pos-1=C → CpG context on bottom strand)
        let ref_bases = FxHashMap::from_iter([
            (6103079, G),
            (6103080, C),
            (6103081, G),
            (6103082, T),
            (6103091, A),
            (6103092, C),
            (6103093, G),
            (6103094, C),
        ]);
        let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };

        // === Legacy mode ===
        let mut legacy_record = record.clone();
        rewrite_record(&calls, &mut legacy_record, BamMode::Legacy, &ref_lookup)?;

        let Aux::String(xm_tag) = legacy_record.aux(b"XM")? else {
            bail!("XM not a string");
        };
        let Aux::String(xr_tag) = legacy_record.aux(b"XR")? else {
            bail!("XR not a string");
        };
        let Aux::String(xg_tag) = legacy_record.aux(b"XG")? else {
            bail!("XG not a string");
        };

        // Flag 163: second in pair → XR=GA; OB → XG=GA
        assert_snapshot!(xr_tag, @"GA");
        assert_snapshot!(xg_tag, @"GA");

        // XM: position 5 = Z (methylated CpG), position 17 = z (unmethylated CpG), rest = '.'
        assert_eq!(xm_tag.len(), seq.len());
        assert_snapshot!(xm_tag, @".....Z...........z............................................................");

        // Verify no stray annotations
        let xm_annotations: Vec<(usize, char)> =
            xm_tag.chars().enumerate().filter(|(_, c)| *c != '.').collect();
        assert_compact_debug_snapshot!(xm_annotations, @"[(5, 'Z'), (17, 'z')]");

        // SEQ should NOT be rewritten in legacy mode
        assert_eq!(legacy_record.seq().as_bytes(), seq, "Legacy mode should not modify SEQ");

        // === Standard mode ===
        let mut standard_record = record.clone();
        rewrite_record(&calls, &mut standard_record, BamMode::Standard, &ref_lookup)?;

        let Aux::String(mm_tag) = standard_record.aux(b"MM")? else {
            bail!("MM not a string");
        };

        // OB strand, not reversed → MM uses G-m format with stored-SEQ positions
        // Rewritten SEQ: position 5 changed from A→G (un-modify methylation evidence)
        let new_seq = standard_record.seq().as_bytes();
        assert_snapshot!(as_base_string(&new_seq), @"AAGGCGTGCACCACCACGCCTGGCTTGGTTTGGTTTTTGATTGGTTGGTTGGTCTTTTGAGACAGGGTTTCTCTGTGT");

        // Verify position 5 was rewritten (A→G) and position 17 unchanged (G stays G)
        assert_eq!(Base::from(seq[5]), A, "original should have A at pos 5");
        assert_eq!(Base::from(new_seq[5]), G, "rewritten should have G at pos 5");
        assert_eq!(Base::from(new_seq[17]), G, "pos 17 should stay G (unmethylated)");

        // MM tag: G-m encoding with skip list
        assert_snapshot!(mm_tag, @"G-m,2;");

        // === Cross-check: both modes agree on methylated positions ===
        let xm_methylated = decode_xm_to_positions(xm_tag);

        let fwd_seq = seq_for_mm_tag(&standard_record);
        let (mm_base, mm_methylated) = decode_mm_to_positions(mm_tag, &fwd_seq)?;
        assert_eq!(mm_base, G, "MM should target G for OB strand");

        // Both should identify position 5 as the only methylated position
        assert_eq!(xm_methylated, vec![5], "XM methylated positions");
        assert_eq!(mm_methylated, vec![5], "MM methylated positions");

        // Standard mode should NOT have XR/XG/XM tags
        assert!(standard_record.aux(b"XR").is_err());
        // Legacy mode should NOT have MM/ML tags
        assert!(legacy_record.aux(b"MM").is_err());

        Ok(())
    }

    /// XM and MM/ML tags must reflect per-read methylation evidence, not the
    /// position-level call. A CpG with beta ≤ 0.5 (`methylated: false`) still
    /// has individually methylated reads that must show Z (not z) in XM and
    /// appear in MM/ML.
    #[test]
    fn per_read_methylation_independent_of_position_call() -> Result<()> {
        use Base::*;

        // Same OB flag-163 read as golden_legacy_and_standard_agree
        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
        bam.fetch(FetchDefinition::All)?;
        let mut record = Record::new();
        while let Some(result) = bam.read(&mut record) {
            result?;
            if record.flags() == 163 && record.pos() == 6103076 {
                break;
            }
        }
        ensure!(record.flags() == 163, "could not find flag-163 record at 6103076");

        let strand = StrandFromRecord::strand(&record);
        assert_eq!(strand, Strand::OB);

        // Position 6103081 (read offset 5): Read[5] = A → evidence base for OB
        // Position 6103093 (read offset 17): Read[17] = G → target base
        //
        // Mark BOTH as `methylated: false` (simulating beta ≤ 0.5).
        // Despite the position-level call, read offset 5 shows methylation
        // evidence (A at G ref) and must be annotated as Z, not z.
        let calls = FxHashMap::from_iter([
            (6103081, RastairCall::Cpg { base: G, methylated: false }),
            (6103093, RastairCall::Cpg { base: G, methylated: false }),
        ]);

        let ref_bases =
            FxHashMap::from_iter([(6103080, C), (6103081, G), (6103092, C), (6103093, G)]);
        let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };

        // === Legacy mode ===
        let mut legacy_record = record.clone();
        rewrite_record(&calls, &mut legacy_record, BamMode::Legacy, &ref_lookup)?;

        let Aux::String(xm_tag) = legacy_record.aux(b"XM")? else {
            bail!("XM not a string");
        };

        // Offset 5 must be Z (methylated per-read), offset 17 must be z (unmethylated per-read)
        let xm_annotations: Vec<(usize, char)> =
            xm_tag.chars().enumerate().filter(|(_, c)| *c != '.').collect();
        assert_compact_debug_snapshot!(xm_annotations, @"[(5, 'Z'), (17, 'z')]");

        // === Standard mode ===
        let mut standard_record = record.clone();
        rewrite_record(&calls, &mut standard_record, BamMode::Standard, &ref_lookup)?;

        let new_seq = standard_record.seq().as_bytes();
        // Position 5 must be rewritten A→G (undo methylation evidence)
        assert_eq!(Base::from(new_seq[5]), G, "position 5 should be rewritten A→G");

        let Aux::String(mm_tag) = standard_record.aux(b"MM")? else {
            bail!("MM not a string");
        };
        // MM tag should report the methylation at position 5
        let fwd_seq = seq_for_mm_tag(&standard_record);
        let (mm_base, mm_methylated) = decode_mm_to_positions(mm_tag, &fwd_seq)?;
        assert_eq!(mm_base, G);
        assert_eq!(mm_methylated, vec![5], "MM should report methylation at offset 5");

        Ok(())
    }

    /// Golden test for an OT strand read (flag 147, reversed, second in pair).
    /// Verifies both modes produce correct output for the reverse-strand case
    /// where MM tag positions must be flipped to original-read orientation.
    #[test]
    fn golden_ot_reversed_legacy_and_standard_agree() -> Result<()> {
        use Base::*;

        // Read 1: flag 83, OB strand, first in pair, REVERSED
        // BAM pos: 6103075 (0-based)
        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
        bam.fetch((0, 6103075, 6103100))?;
        let mut record = Record::new();
        bam.read(&mut record).wrap_err("no records")?.wrap_err("read")?;

        assert_eq!(record.flags(), 83);
        let strand = StrandFromRecord::strand(&record);
        assert_eq!(strand, Strand::OB);
        assert!(record.is_reverse());
        let seq = record.seq().as_bytes();

        // Two CpG calls (G side for OB):
        // Pos 6103081 (read offset 6): Read[6]=G → target, unmethylated in this read
        // Pos 6103093 (read offset 18): Read[18]=G → target, unmethylated in this read
        //
        // Also add a methylated call where read actually shows evidence:
        // Let's find a position where the read has A (OB evidence)
        let mut evidence_pos = None;
        for [pos_in_read, pos_in_ref] in record.aligned_pairs() {
            if Base::from(seq[pos_in_read as usize]) == A {
                evidence_pos = Some((pos_in_read as usize, pos_in_ref as u32));
                break;
            }
        }
        let (ev_read_pos, ev_ref_pos) = evidence_pos.wrap_err("no A base found in read")?;

        let calls = FxHashMap::from_iter([
            (6103081, RastairCall::Cpg { base: G, methylated: false }),
            (6103093, RastairCall::Cpg { base: G, methylated: false }),
            (ev_ref_pos, RastairCall::Cpg { base: G, methylated: true }),
        ]);

        let ref_bases = FxHashMap::from_iter([
            (6103080, C),
            (6103081, G),
            (6103092, C),
            (6103093, G),
            (ev_ref_pos.saturating_sub(1), C),
            (ev_ref_pos, G),
        ]);
        let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };

        // === Legacy mode ===
        let mut legacy_record = record.clone();
        rewrite_record(&calls, &mut legacy_record, BamMode::Legacy, &ref_lookup)?;

        let Aux::String(xm_tag) = legacy_record.aux(b"XM")? else {
            bail!("XM not a string");
        };
        assert_eq!(xm_tag.len(), seq.len());

        // Verify the three annotated positions
        let xm_annotations: Vec<(usize, char)> =
            xm_tag.chars().enumerate().filter(|(_, c)| *c != '.').collect();
        // ev_read_pos should be Z (methylated), positions 6 and 18 should be z (unmethylated)
        let mut expected = vec![(6, 'z'), (18, 'z'), (ev_read_pos, 'Z')];
        expected.sort();
        assert_eq!(xm_annotations, expected, "XM annotations mismatch");

        // === Standard mode ===
        let mut standard_record = record.clone();
        rewrite_record(&calls, &mut standard_record, BamMode::Standard, &ref_lookup)?;

        let Aux::String(mm_tag) = standard_record.aux(b"MM")? else {
            bail!("MM not a string");
        };

        // This is a reversed read, so MM positions are in original-read orientation
        // (reverse complement of stored SEQ). The strand qualifier flips OB→OT.
        assert_snapshot!(mm_tag.chars().take(3).collect::<String>(), @"C+m");

        // === Cross-check ===
        let xm_methylated = decode_xm_to_positions(xm_tag);
        assert_eq!(xm_methylated, vec![ev_read_pos], "XM: only evidence position should be Z");

        let fwd_seq = seq_for_mm_tag(&standard_record);
        let (mm_base, mm_methylated) = decode_mm_to_positions(mm_tag, &fwd_seq)?;
        // For reversed OB→OT, MM uses C base
        assert_eq!(mm_base, C);

        // The methylated position in stored-SEQ is ev_read_pos.
        // In original-read (RC) orientation: seq_len - 1 - ev_read_pos
        let expected_mm_pos = seq.len() - 1 - ev_read_pos;
        assert_eq!(
            mm_methylated,
            vec![expected_mm_pos],
            "MM: methylated position in original-read coords"
        );

        Ok(())
    }

    #[test]
    fn sequence_transform() -> Result<()> {
        use Base::*;

        // it's easiest to take a record from a real BAM file instead of
        // constructing one manually
        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")
            .wrap_err("failed to open test BAM")?;
        bam.fetch((0, 6103075, 6103100))?;
        let mut record = Record::new();
        bam.read(&mut record).wrap_err("no records")?.wrap_err("failed to read record")?;

        // just to be sure we have the expected record
        assert_snapshot!(record.pos(), @"6103075");
        assert_snapshot!(StrandFromRecord::strand(&record), @"OB");

        let ol_seq = record.seq().as_bytes();

        // We want to create some test calls that are setting some CpG positions
        // to be methylated.
        //
        // In the sequence these are not encoded as `CG` but as either `TG` (top
        // strand) or `CA` (bottom strand). our test record is on the bottom
        // strand, so we'll look for `CA` in the sequence and add a call entry
        // for it. For good measure, we'll also add calls for CG since the real
        // calls will likely contain them.
        let calls = {
            let mut calls = FxHashMap::default();
            let seq = &ol_seq;

            for [pos_in_read, pos_in_ref] in record.aligned_pairs() {
                let current_base = Base::from(seq[pos_in_read as usize]);
                let Some(next_base) = seq.get((pos_in_read + 1) as usize).map(Base::from) else {
                    continue;
                };

                match (current_base, next_base) {
                    (C, G) => {
                        calls.insert(
                            pos_in_ref as u32,
                            RastairCall::Cpg { methylated: true, base: C },
                        );
                    }
                    (C, A) => {
                        calls.insert(
                            pos_in_ref as u32 + 1,
                            RastairCall::Cpg { methylated: true, base: G },
                        );
                    }
                    _ => {}
                }
            }
            calls
        };

        let data = get_methylated_positions(&calls, &record)?;
        let new_seq = &data.seq;

        // for human comparison
        assert_snapshot!(as_base_string(&ol_seq), @"AAAGGCGTGCACCACCACGCCTGGCTTGGTTTGGTTTTTGATTGGTTGGTTGGTCTTTTGAGACAGGGTTTCTCTGTGTA");
        assert_snapshot!(as_base_string(new_seq), @"AAAGGCGTGCGCCGCCGCGCCTGGCTTGGTTTGGTTTTTGATTGGTTGGTTGGTCTTTTGAGACGGGGTTTCTCTGTGTA");

        assert_compact_debug_snapshot!(data.methylated_positions, @"[10, 13, 16, 64]");
        // verify these are all positions where A (methylation evidence) was
        // rewritten to G
        for &pos in &data.methylated_positions {
            let base = Base::from(ol_seq[pos as usize]);
            if base != A {
                bail!("expected A at methylated position {}, found {}", pos, base.as_str());
            }
            let new_base = Base::from(new_seq[pos as usize]);
            if new_base != G {
                bail!(
                    "expected G at rewritten methylated position {}, found {}",
                    pos,
                    new_base.as_str()
                );
            }
        }

        // This is a reverse OB read (flag=83). The MM tag must be in
        // original-read orientation (RC of stored SEQ), with flipped
        // base/strand (OB→OT = C+m).
        let is_reverse = record.is_reverse();
        assert!(is_reverse, "test record should be reverse-strand");

        let seq_len = new_seq.len();
        let original_positions: SmallVec<u32, 10> = data
            .methylated_positions
            .iter()
            .map(|&p| u32::try_from(seq_len - 1 - p as usize).expect("fits"))
            .collect();
        let original_seq = reverse_complement(new_seq);

        let methylated_positions =
            MethylatedPositions::new(Strand::OT, &original_seq, &original_positions);

        assert_compact_debug_snapshot!(methylated_positions, @"MethylatedPositions { base: C, strand: OT, positions: [5, 16, 0, 0] }");
        methylated_positions
            .apply_to_record(&mut record)
            .wrap_err("failed to apply modifications to record")?;

        let Aux::String(mod_string) = record.aux(b"MM").wrap_err("missing MM tag")? else {
            bail!("MM tag is not a string");
        };
        assert_snapshot!(mod_string, @"C+m,5,16,0,0;");

        Ok(())
    }

    #[test]
    fn test_xr_xg_xm_tags() -> Result<()> {
        use Base::*;

        // Test with flag 83 (OB strand, first in pair reverse)
        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")
            .wrap_err("failed to open test BAM")?;
        bam.fetch((0, 6103075, 6103100))?;
        let mut record = Record::new();
        bam.read(&mut record).wrap_err("no records")?.wrap_err("failed to read record")?;

        // Verify this is flag 83
        assert_eq!(record.flags(), 83, "Expected flag 83");
        assert_eq!(StrandFromRecord::strand(&record), Strand::OB);

        let ol_seq = record.seq().as_bytes();

        // Create test calls with some methylated CpGs
        let calls = {
            let mut calls = FxHashMap::default();
            let seq = &ol_seq;

            for [pos_in_read, pos_in_ref] in record.aligned_pairs() {
                let current_base = Base::from(seq[pos_in_read as usize]);
                let Some(next_base) = seq.get((pos_in_read + 1) as usize).map(Base::from) else {
                    continue;
                };

                match (current_base, next_base) {
                    (C, G) => {
                        calls.insert(
                            pos_in_ref as u32,
                            RastairCall::Cpg { methylated: true, base: C },
                        );
                    }
                    (C, A) => {
                        calls.insert(
                            pos_in_ref as u32 + 1,
                            RastairCall::Cpg { methylated: true, base: G },
                        );
                    }
                    _ => {}
                }
            }
            calls
        };

        rewrite_record(&calls, &mut record, BamMode::Legacy, no_ref)?;

        // Check XR tag (flag 83: first in pair => CT)
        let Aux::String(xr_tag) = record.aux(b"XR").wrap_err("missing XR tag")? else {
            bail!("XR tag is not a string");
        };
        assert_snapshot!(xr_tag, @"CT");

        // Check XG tag (OB strand => GA)
        let Aux::String(xg_tag) = record.aux(b"XG").wrap_err("missing XG tag")? else {
            bail!("XG tag is not a string");
        };
        assert_snapshot!(xg_tag, @"GA");

        // Check XM tag exists and has correct length
        let Aux::String(xm_tag) = record.aux(b"XM").wrap_err("missing XM tag")? else {
            bail!("XM tag is not a string");
        };
        assert_eq!(xm_tag.len(), record.seq_len());

        Ok(())
    }

    #[test]
    fn test_xr_xg_tags_all_flags() -> Result<()> {
        let test_cases = vec![
            (99, "CT", "CT"),  // flag 99: OT, first in pair => XR:CT, XG:CT
            (147, "GA", "CT"), // flag 147: OT, second in pair => XR:GA, XG:CT
            (83, "CT", "GA"),  // flag 83: OB, first in pair => XR:CT, XG:GA
            (163, "GA", "GA"), // flag 163: OB, second in pair => XR:GA, XG:GA
        ];

        for (flag, expected_xr, expected_xg) in test_cases {
            let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")
                .wrap_err("failed to open test BAM")?;

            bam.fetch(FetchDefinition::All)?;
            let mut record = Record::new();
            let mut found = false;

            while let Some(result) = bam.read(&mut record) {
                result?;
                if record.flags() == flag {
                    found = true;
                    break;
                }
            }

            ensure!(found, "No record with flag {flag} found in test BAM");

            let calls = FxHashMap::default();
            rewrite_record(&calls, &mut record, BamMode::Legacy, no_ref)?;

            let Aux::String(xr_tag) = record.aux(b"XR").wrap_err("missing XR tag")? else {
                bail!("XR tag is not a string for flag {}", flag);
            };
            assert_eq!(xr_tag, expected_xr, "XR tag mismatch for flag {}", flag);

            let Aux::String(xg_tag) = record.aux(b"XG").wrap_err("missing XG tag")? else {
                bail!("XG tag is not a string for flag {}", flag);
            };
            assert_eq!(xg_tag, expected_xg, "XG tag mismatch for flag {}", flag);
        }

        Ok(())
    }

    /// Decode an MM tag string back to absolute sequence positions.
    ///
    /// MM format: `BASE[+-]m,skip1,skip2,...;`
    /// The skips are deltas between consecutive occurrences of the target base.
    ///
    /// Both MM and XM tags in this codebase are generated from
    /// `seq_for_mm_tag` (the forward-oriented sequence). To decode MM
    /// consistently, pass the same forward-oriented sequence.
    fn decode_mm_to_positions(mm_tag: &str, seq: &[u8]) -> Result<(Base, Vec<usize>)> {
        let mm_tag = mm_tag.trim_end_matches(';');
        ensure!(mm_tag.len() >= 3, "MM tag too short: {mm_tag}");

        let base = Base::from(mm_tag.as_bytes()[0]);

        // Format: "C+m,skip1,skip2,..." or "G-m;" (no positions)
        // Skip past "C+m" prefix; if there's a comma, skip that too
        if mm_tag.len() <= 3 {
            return Ok((base, vec![]));
        }
        let skips_str = &mm_tag[4..]; // skip "C+m," or "G-m,"

        if skips_str.is_empty() {
            return Ok((base, vec![]));
        }

        let base_positions: Vec<usize> = seq
            .iter()
            .enumerate()
            .filter(|&(_, b)| Base::from(*b) == base)
            .map(|(i, _)| i)
            .collect();

        let mut positions = Vec::new();
        let mut base_idx: usize = 0;

        for skip in skips_str.split(',') {
            let skip: usize = skip.parse().wrap_err_with(|| format!("bad skip: {skip}"))?;
            base_idx += skip;
            let seq_pos = base_positions.get(base_idx).ok_or_else(|| {
                color_eyre::eyre::eyre!("MM skip out of bounds: base_idx={base_idx}")
            })?;
            positions.push(*seq_pos);
            base_idx += 1;
        }

        Ok((base, positions))
    }

    /// Decode an XM tag string to absolute sequence positions of methylated bases.
    fn decode_xm_to_positions(xm_tag: &str) -> Vec<usize> {
        xm_tag.chars().enumerate().filter(|(_, c)| *c == 'Z').map(|(i, _)| i).collect()
    }

    /// Get the sequence in MM tag orientation. For forward reads this is the
    /// stored SEQ; for reverse reads it is the reverse complement (original
    /// read orientation), since MM positions count from the original 5' end.
    fn seq_for_mm_tag(record: &Record) -> Vec<u8> {
        let stored = record.seq().as_bytes();
        if record.is_reverse() { reverse_complement(&stored) } else { stored }
    }

    #[test]
    fn roundtrip_standard_mm_tags() -> Result<()> {
        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")
            .wrap_err("failed to open test BAM")?;

        bam.fetch((0, 6103075, 6103200))?;
        let mut record = Record::new();
        let mut tested = 0u32;

        while let Some(result) = bam.read(&mut record) {
            result?;
            let ol_seq = record.seq().as_bytes();
            let calls = build_cpg_calls(&record, &ol_seq);
            rewrite_record(&calls, &mut record, BamMode::Standard, no_ref)?;

            // Standard mode should NOT produce XR/XG/XM tags
            assert!(record.aux(b"XR").is_err(), "Standard mode should not produce XR tag");

            // MM/ML are omitted when there are no methylated positions (absent = "no data").
            if let Ok(aux) = record.aux(b"MM") {
                let Aux::String(mm_tag) = aux else {
                    bail!("MM tag is not a string at pos={} flag={}", record.pos(), record.flags());
                };

                let fwd_seq = seq_for_mm_tag(&record);
                let (_, mm_positions) =
                    decode_mm_to_positions(mm_tag, &fwd_seq).wrap_err_with(|| {
                        format!(
                            "decode MM {:?} pos={} flag={}",
                            mm_tag,
                            record.pos(),
                            record.flags()
                        )
                    })?;
                assert!(mm_positions.len() <= fwd_seq.len(), "MM positions out of range");
            }

            tested += 1;
        }

        ensure!(tested > 0, "no records tested");
        Ok(())
    }

    #[test]
    fn roundtrip_write_read_bam_standard() -> Result<()> {
        let temp_dir = tempfile::TempDir::new()?;
        let temp_bam = temp_dir.path().join("roundtrip.bam");

        let header = {
            let bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
            Header::from_template(bam.header())
        };

        let mut writer =
            Writer::from_path(&temp_bam, &header, bam::Format::Bam).wrap_err("create writer")?;

        let mut bam =
            bam::IndexedReader::from_path("tests/data/test.bam").wrap_err("open test BAM")?;
        bam.fetch((0, 6103075, 6103100))?;

        let mut records_written = 0u32;
        let mut record = Record::new();
        while let Some(result) = bam.read(&mut record) {
            result?;
            let ol_seq = record.seq().as_bytes();
            let calls = build_cpg_calls(&record, &ol_seq);
            rewrite_record(&calls, &mut record, BamMode::Standard, no_ref)?;
            writer.write(&record)?;
            records_written += 1;
        }
        drop(writer);
        ensure!(records_written > 0, "no records written");

        let mut reader =
            bam::Reader::from_path(&temp_bam).wrap_err("open roundtrip BAM for reading")?;
        let mut records_read = 0u32;

        while let Some(result) = reader.read(&mut record) {
            result?;
            records_read += 1;

            let fwd_seq = seq_for_mm_tag(&record);

            // MM/ML are omitted when there are no methylated positions (absent = "no data").
            if let Ok(aux) = record.aux(b"MM") {
                let Aux::String(mm_tag) = aux else {
                    bail!("MM not a string");
                };

                let (_, mm_positions) = decode_mm_to_positions(mm_tag, &fwd_seq)?;

                let Aux::ArrayU8(ml_data) = record.aux(b"ML").wrap_err("missing ML")? else {
                    bail!("ML not an array");
                };
                assert_eq!(
                    ml_data.iter().count(),
                    mm_positions.len(),
                    "ML length should match number of methylated positions"
                );
            }

            assert!(record.aux(b"XR").is_err(), "Standard mode should not produce XR");
        }

        assert_eq!(records_written, records_read, "record count mismatch after roundtrip");

        Ok(())
    }

    #[test]
    fn roundtrip_write_read_bam_legacy() -> Result<()> {
        let temp_dir = tempfile::TempDir::new()?;
        let temp_bam = temp_dir.path().join("roundtrip_legacy.bam");

        let header = {
            let bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
            Header::from_template(bam.header())
        };

        let mut writer =
            Writer::from_path(&temp_bam, &header, bam::Format::Bam).wrap_err("create writer")?;

        let mut bam =
            bam::IndexedReader::from_path("tests/data/test.bam").wrap_err("open test BAM")?;
        bam.fetch((0, 6103075, 6103100))?;

        let mut records_written = 0u32;
        let mut record = Record::new();
        while let Some(result) = bam.read(&mut record) {
            result?;
            let ol_seq = record.seq().as_bytes();
            let calls = build_cpg_calls(&record, &ol_seq);
            rewrite_record(&calls, &mut record, BamMode::Legacy, no_ref)?;
            writer.write(&record)?;
            records_written += 1;
        }
        drop(writer);
        ensure!(records_written > 0, "no records written");

        let mut reader =
            bam::Reader::from_path(&temp_bam).wrap_err("open roundtrip BAM for reading")?;
        let mut records_read = 0u32;

        while let Some(result) = reader.read(&mut record) {
            result?;
            records_read += 1;

            let Aux::String(_) = record.aux(b"XR").wrap_err("missing XR")? else {
                bail!("XR not a string");
            };
            let Aux::String(_) = record.aux(b"XG").wrap_err("missing XG")? else {
                bail!("XG not a string");
            };
            let Aux::String(xm_tag) = record.aux(b"XM").wrap_err("missing XM")? else {
                bail!("XM not a string");
            };
            assert_eq!(xm_tag.len(), record.seq_len(), "XM length mismatch");

            assert!(record.aux(b"MM").is_err(), "Legacy mode should not produce MM");
        }

        assert_eq!(records_written, records_read, "record count mismatch after roundtrip");

        Ok(())
    }

    #[test]
    fn roundtrip_no_methylation_standard() -> Result<()> {
        let temp_dir = tempfile::TempDir::new()?;
        let temp_bam = temp_dir.path().join("roundtrip_empty.bam");

        let header = {
            let bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
            Header::from_template(bam.header())
        };

        let mut writer = Writer::from_path(&temp_bam, &header, bam::Format::Bam)?;

        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
        bam.fetch((0, 6103075, 6103100))?;

        let mut record = Record::new();
        bam.read(&mut record).wrap_err("no records")?.wrap_err("read")?;

        let calls = FxHashMap::default();
        rewrite_record(&calls, &mut record, BamMode::Standard, no_ref)?;
        writer.write(&record)?;
        drop(writer);

        let mut reader = bam::Reader::from_path(&temp_bam)?;
        reader.read(&mut record).wrap_err("no records")?.wrap_err("read back")?;

        // With no methylated positions, MM/ML tags are not written at all.
        // An absent MM means "no modification data" per SAM spec.
        // Writing empty C+m; or empty ML:B:C crashes tools like modbedtools.
        assert!(
            record.aux(b"MM").is_err(),
            "MM tag should not be present when there are no methylated positions"
        );
        assert!(
            record.aux(b"ML").is_err(),
            "ML tag should not be present when there are no methylated positions"
        );

        Ok(())
    }

    #[test]
    fn roundtrip_no_methylation_legacy() -> Result<()> {
        let temp_dir = tempfile::TempDir::new()?;
        let temp_bam = temp_dir.path().join("roundtrip_empty_legacy.bam");

        let header = {
            let bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
            Header::from_template(bam.header())
        };

        let mut writer = Writer::from_path(&temp_bam, &header, bam::Format::Bam)?;

        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
        bam.fetch((0, 6103075, 6103100))?;

        let mut record = Record::new();
        bam.read(&mut record).wrap_err("no records")?.wrap_err("read")?;

        let calls = FxHashMap::default();
        rewrite_record(&calls, &mut record, BamMode::Legacy, no_ref)?;
        writer.write(&record)?;
        drop(writer);

        let mut reader = bam::Reader::from_path(&temp_bam)?;
        reader.read(&mut record).wrap_err("no records")?.wrap_err("read back")?;

        let Aux::String(xm_tag) = record.aux(b"XM").wrap_err("missing XM")? else {
            bail!("XM not a string");
        };
        let xm_positions = decode_xm_to_positions(xm_tag);
        assert!(xm_positions.is_empty(), "expected no XM methylation");

        ensure!(!xm_tag.contains('Z'), "XM should have no uppercase Z without methylation calls");

        Ok(())
    }

    /// Build CpG methylation calls from a record's sequence by looking for
    /// CA (OB strand evidence) and CG dinucleotides.
    fn build_cpg_calls(record: &Record, seq: &[u8]) -> FxHashMap<u32, RastairCall> {
        use Base::*;
        let mut calls = FxHashMap::default();
        for [pos_in_read, pos_in_ref] in record.aligned_pairs() {
            let current_base = Base::from(seq[pos_in_read as usize]);
            let Some(next_base) = seq.get((pos_in_read + 1) as usize).map(Base::from) else {
                continue;
            };
            match (current_base, next_base) {
                (C, G) => {
                    calls.insert(pos_in_ref as u32, RastairCall::Cpg { methylated: true, base: C });
                }
                (C, A) => {
                    calls.insert(
                        pos_in_ref as u32 + 1,
                        RastairCall::Cpg { methylated: true, base: G },
                    );
                }
                _ => {}
            }
        }
        calls
    }

    fn as_base_string(seq: &[u8]) -> String {
        seq.iter().map(|b| Base::from(*b).as_str()).collect()
    }

    /// Only positions that are actual CpG sites in the calls should be marked
    /// z/Z in the XM tag. Every other C (OT) or G (OB) that is NOT a CpG in
    /// the reference must be `.`.
    #[test]
    fn xm_only_marks_called_cpg_positions() -> Result<()> {
        use Base::*;

        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
        bam.fetch((0, 6103075, 6103100))?;
        let mut record = Record::new();
        bam.read(&mut record).wrap_err("no records")?.wrap_err("read")?;

        let strand = StrandFromRecord::strand(&record);
        let seq = record.seq().as_bytes();

        // Count how many G bases (OB target) are in the read
        let total_target_bases = seq
            .iter()
            .filter(|&&b| Base::from(b) == if strand == Strand::OB { G } else { C })
            .count();

        // Create calls for only ONE CpG position — pick the first CG or CA
        // dinucleotide we find
        let mut calls = FxHashMap::default();
        let mut cpg_count = 0u32;
        for [pos_in_read, pos_in_ref] in record.aligned_pairs() {
            let current_base = Base::from(seq[pos_in_read as usize]);
            let Some(next_base) = seq.get((pos_in_read + 1) as usize).map(Base::from) else {
                continue;
            };

            // For OB: look for CA (methylation evidence for G on next pos)
            if strand == Strand::OB && current_base == C && next_base == A && cpg_count == 0 {
                calls.insert(pos_in_ref as u32 + 1, RastairCall::Cpg { methylated: true, base: G });
                cpg_count += 1;
                break;
            }
            // For OT: look for TG (methylation evidence for C)
            if strand == Strand::OT && current_base == T && next_base == G && cpg_count == 0 {
                calls.insert(pos_in_ref as u32, RastairCall::Cpg { methylated: true, base: C });
                cpg_count += 1;
                break;
            }
        }
        ensure!(cpg_count == 1, "need at least one CpG call for this test");

        rewrite_record(&calls, &mut record, BamMode::Legacy, no_ref)?;

        let Aux::String(xm_tag) = record.aux(b"XM")? else {
            bail!("XM not a string");
        };

        let z_count = xm_tag.chars().filter(|c| *c == 'Z' || *c == 'z').count();

        assert!(
            z_count < total_target_bases,
            "XM has {z_count} CpG annotations but only {cpg_count} CpG call(s) exist. \
             There are {total_target_bases} target bases in the read — non-CpG positions \
             should be '.' not 'z'. XM: {xm_tag}"
        );

        Ok(())
    }

    /// De-novo CpG calls (DeNovoCpg variant) must also appear in legacy XM tags.
    #[test]
    fn xm_includes_denovo_cpg_methylation() -> Result<()> {
        use Base::*;

        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
        bam.fetch((0, 6103075, 6103100))?;
        let mut record = Record::new();
        bam.read(&mut record).wrap_err("no records")?.wrap_err("read")?;

        let strand = StrandFromRecord::strand(&record);
        let seq = record.seq().as_bytes();

        // Create a single DeNovoCpg call at the first methylation evidence base
        let mut calls = FxHashMap::default();
        let mut found = false;
        for [pos_in_read, pos_in_ref] in record.aligned_pairs() {
            let observed = Base::from(seq[pos_in_read as usize]);
            let is_evidence = match strand {
                Strand::OT => observed == T,
                Strand::OB => observed == A,
                Strand::Unknown => false,
            };
            if is_evidence {
                let base = match strand {
                    Strand::OT => C,
                    Strand::OB => G,
                    Strand::Unknown => unreachable!(),
                };
                calls.insert(pos_in_ref as u32, RastairCall::DeNovoCpg { methylated: true, base });
                found = true;
                break;
            }
        }
        ensure!(found, "need a methylation-evidence base for this test");

        rewrite_record(&calls, &mut record, BamMode::Legacy, no_ref)?;

        let Aux::String(xm_tag) = record.aux(b"XM")? else {
            bail!("XM not a string");
        };

        let has_methylated = xm_tag.chars().any(|c| c == 'Z' || c == 'X' || c == 'H');
        assert!(
            has_methylated,
            "XM tag should contain a methylation mark for the DeNovoCpg call, \
             but got: {xm_tag}"
        );

        Ok(())
    }

    /// De-novo CpG calls must also be found by find_methylated_positions (standard mode).
    #[test]
    fn standard_mode_includes_denovo_cpg() -> Result<()> {
        use Base::*;

        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
        bam.fetch((0, 6103075, 6103100))?;
        let mut record = Record::new();
        bam.read(&mut record).wrap_err("no records")?.wrap_err("read")?;

        let strand = StrandFromRecord::strand(&record);
        let seq = record.seq().as_bytes();

        let mut calls = FxHashMap::default();
        let mut expected_pos = None;
        for [pos_in_read, pos_in_ref] in record.aligned_pairs() {
            let observed = Base::from(seq[pos_in_read as usize]);
            let is_evidence = match strand {
                Strand::OT => observed == T,
                Strand::OB => observed == A,
                Strand::Unknown => false,
            };
            if is_evidence {
                let base = match strand {
                    Strand::OT => C,
                    Strand::OB => G,
                    Strand::Unknown => unreachable!(),
                };
                calls.insert(pos_in_ref as u32, RastairCall::DeNovoCpg { methylated: true, base });
                expected_pos = Some(pos_in_read as u32);
                break;
            }
        }
        ensure!(expected_pos.is_some(), "need a methylation-evidence base");

        let data = get_methylated_positions(&calls, &record)?;
        assert!(
            !data.methylated_positions.is_empty(),
            "get_methylated_positions should find DeNovoCpg calls, but found none"
        );

        Ok(())
    }

    /// De-novo CpGs should use context-dependent letters (x/X for CHG, h/H for CHH)
    /// rather than z/Z (which is reserved for reference CpG context).
    #[test]
    fn xm_denovo_uses_context_letters() -> Result<()> {
        use Base::*;

        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
        bam.fetch((0, 6103075, 6103100))?;
        let mut record = Record::new();
        bam.read(&mut record).wrap_err("no records")?.wrap_err("read")?;

        let strand = StrandFromRecord::strand(&record);
        let seq = record.seq().as_bytes();

        // Create a DeNovoCpg call and a synthetic ref lookup that returns
        // non-CpG context (CHH: no G after C for OT, no C before G for OB)
        let mut calls = FxHashMap::default();
        let mut denovo_ref_pos = None;
        for [pos_in_read, pos_in_ref] in record.aligned_pairs() {
            let observed = Base::from(seq[pos_in_read as usize]);
            let is_evidence = match strand {
                Strand::OT => observed == T,
                Strand::OB => observed == A,
                Strand::Unknown => false,
            };
            if is_evidence {
                let base = match strand {
                    Strand::OT => C,
                    Strand::OB => G,
                    Strand::Unknown => unreachable!(),
                };
                calls.insert(pos_in_ref as u32, RastairCall::DeNovoCpg { methylated: true, base });
                denovo_ref_pos = Some(pos_in_ref as u32);
                break;
            }
        }
        let denovo_ref_pos = denovo_ref_pos.wrap_err("no evidence base found")?;

        // Synthetic reference that returns A at all neighboring positions → CHH context
        let ref_lookup = |pos: u32| -> Option<Base> {
            if pos == denovo_ref_pos {
                match strand {
                    Strand::OT => Some(C),
                    Strand::OB => Some(G),
                    Strand::Unknown => None,
                }
            } else {
                Some(A) // non-CpG neighbor → CHH
            }
        };

        rewrite_record(&calls, &mut record, BamMode::Legacy, ref_lookup)?;

        let Aux::String(xm_tag) = record.aux(b"XM")? else {
            bail!("XM not a string");
        };

        // De-novo with CHH context should use H (methylated CHH), not Z
        assert!(
            xm_tag.contains('H'),
            "De-novo CpG in CHH context should produce 'H', got: {xm_tag}"
        );
        assert!(
            !xm_tag.contains('Z') && !xm_tag.contains('z'),
            "De-novo CpG should not use z/Z (CpG context letters), got: {xm_tag}"
        );

        Ok(())
    }

    /// Golden test for an OT strand read (flag 99, forward, first in pair).
    /// OT targets C positions with T as methylation evidence.
    #[test]
    fn golden_ot_forward_legacy_and_standard_agree() -> Result<()> {
        use Base::*;

        // Flag 99: OT, first in pair, NOT reversed
        // SAM pos 6103079 (1-based) → BAM pos 6103078 (0-based)
        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
        bam.fetch(FetchDefinition::All)?;
        let mut record = Record::new();
        while let Some(result) = bam.read(&mut record) {
            result?;
            if record.flags() == 99 && record.pos() == 6103078 {
                break;
            }
        }
        ensure!(record.flags() == 99, "could not find flag-99 record at 6103078");

        let strand = StrandFromRecord::strand(&record);
        assert_eq!(strand, Strand::OT);
        assert!(!record.is_reverse());
        let seq = record.seq().as_bytes();
        assert_snapshot!(as_base_string(&seq), @"GGTGTGCACCACCATGCCTGGCTTGGTTTGGTTTTTGATTGGTTGGTTGGTCTTTTGAGACAGGGTTTCTCTGTGTAGCT");

        // CpG calls (C side for OT):
        // Ref pos 6103080 (read offset 2): methylated
        //   Read[2] = T → methylation evidence on OT
        // Ref pos 6103092 (read offset 14): methylated
        //   Read[14] = T → methylation evidence on OT
        let calls = FxHashMap::from_iter([
            (6103080, RastairCall::Cpg { base: C, methylated: true }),
            (6103092, RastairCall::Cpg { base: C, methylated: true }),
        ]);

        let ref_bases =
            FxHashMap::from_iter([(6103080, C), (6103081, G), (6103092, C), (6103093, G)]);
        let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };

        // === Legacy mode ===
        let mut legacy_record = record.clone();
        rewrite_record(&calls, &mut legacy_record, BamMode::Legacy, &ref_lookup)?;

        let Aux::String(xm_tag) = legacy_record.aux(b"XM")? else { bail!("XM") };
        let Aux::String(xr_tag) = legacy_record.aux(b"XR")? else { bail!("XR") };
        let Aux::String(xg_tag) = legacy_record.aux(b"XG")? else { bail!("XG") };

        // Flag 99: first in pair → XR=CT; OT → XG=CT
        assert_snapshot!(xr_tag, @"CT");
        assert_snapshot!(xg_tag, @"CT");

        let xm_annotations: Vec<(usize, char)> =
            xm_tag.chars().enumerate().filter(|(_, c)| *c != '.').collect();
        assert_compact_debug_snapshot!(xm_annotations, @"[(2, 'Z'), (14, 'Z')]");

        // SEQ unchanged in legacy mode
        assert_eq!(legacy_record.seq().as_bytes(), seq);

        // === Standard mode ===
        let mut standard_record = record.clone();
        rewrite_record(&calls, &mut standard_record, BamMode::Standard, &ref_lookup)?;

        let Aux::String(mm_tag) = standard_record.aux(b"MM")? else { bail!("MM") };
        // OT, forward → C+m
        assert_snapshot!(mm_tag.chars().take(3).collect::<String>(), @"C+m");
        assert_snapshot!(mm_tag, @"C+m,0,5;");

        // Rewritten SEQ: T→C at methylated positions 2 and 14
        let new_seq = standard_record.seq().as_bytes();
        assert_eq!(Base::from(seq[2]), T);
        assert_eq!(Base::from(new_seq[2]), C);
        assert_eq!(Base::from(seq[14]), T);
        assert_eq!(Base::from(new_seq[14]), C);

        // === Cross-check ===
        let xm_methylated = decode_xm_to_positions(xm_tag);
        let fwd_seq = seq_for_mm_tag(&standard_record);
        let (mm_base, mm_methylated) = decode_mm_to_positions(mm_tag, &fwd_seq)?;
        assert_eq!(mm_base, C);
        assert_eq!(xm_methylated, vec![2, 14]);
        assert_eq!(mm_methylated, vec![2, 14]);

        Ok(())
    }

    /// Golden test for a de-novo CpG: verify both modes produce correct output
    /// and that the XM tag uses context-dependent letters (not z/Z).
    #[test]
    fn golden_denovo_legacy_and_standard_agree() -> Result<()> {
        use Base::*;

        // Use the flag-163 OB read (same as golden_legacy_and_standard_agree)
        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
        bam.fetch(FetchDefinition::All)?;
        let mut record = Record::new();
        while let Some(result) = bam.read(&mut record) {
            result?;
            if record.flags() == 163 && record.pos() == 6103076 {
                break;
            }
        }
        ensure!(record.flags() == 163, "could not find record");

        let strand = StrandFromRecord::strand(&record);
        assert_eq!(strand, Strand::OB);
        let seq = record.seq().as_bytes();

        // Read[5] = A → methylation evidence at ref pos 6103081 on OB
        assert_eq!(Base::from(seq[5]), A);

        // Make it a DeNovoCpg with CHG context:
        // OB context looks at ref[pos-1] and ref[pos-2] (complement perspective)
        // ref[pos-1]=T, ref[pos-2]=C → CHG (second base back is C)
        let calls = FxHashMap::from_iter([(
            6103081_u32,
            RastairCall::DeNovoCpg { base: G, methylated: true },
        )]);
        let ref_bases = FxHashMap::from_iter([
            (6103079_u32, C), // pos-2 → C → CHG
            (6103080, T),     // pos-1 → T (complement=A, not G, so not CpG)
            (6103081, G),
        ]);
        let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };

        // === Legacy mode ===
        let mut legacy_record = record.clone();
        rewrite_record(&calls, &mut legacy_record, BamMode::Legacy, &ref_lookup)?;

        let Aux::String(xm_tag) = legacy_record.aux(b"XM")? else { bail!("XM") };
        let xm_annotations: Vec<(usize, char)> =
            xm_tag.chars().enumerate().filter(|(_, c)| *c != '.').collect();
        // CHG methylated → 'X'
        assert_compact_debug_snapshot!(xm_annotations, @"[(5, 'X')]");

        // === Standard mode ===
        let mut standard_record = record.clone();
        rewrite_record(&calls, &mut standard_record, BamMode::Standard, &ref_lookup)?;

        let Aux::String(mm_tag) = standard_record.aux(b"MM")? else { bail!("MM") };
        // Should still produce MM tag (standard mode doesn't distinguish ref vs de-novo)
        assert_snapshot!(mm_tag, @"G-m,2;");

        // === Cross-check: both agree position 5 is methylated ===
        let xm_methylated: Vec<usize> = xm_tag
            .chars()
            .enumerate()
            .filter(|(_, c)| c.is_ascii_uppercase() && *c != '.')
            .map(|(i, _)| i)
            .collect();
        let fwd_seq = seq_for_mm_tag(&standard_record);
        let (_, mm_methylated) = decode_mm_to_positions(mm_tag, &fwd_seq)?;
        assert_eq!(xm_methylated, vec![5]);
        assert_eq!(mm_methylated, vec![5]);

        Ok(())
    }

    /// De-novo CpGs must produce context-dependent XM letters (z/Z for CpG,
    /// x/X for CHG, h/H for CHH) through the full rewrite_record path, for
    /// both OT and OB strands, and both methylated and unmethylated reads.
    #[test]
    fn denovo_all_contexts_both_strands() -> Result<()> {
        use Base::*;

        // === OB strand (flag-163 read at 6103076) ===
        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
        bam.fetch(FetchDefinition::All)?;
        let mut ob_record = Record::new();
        while let Some(result) = bam.read(&mut ob_record) {
            result?;
            if ob_record.flags() == 163 && ob_record.pos() == 6103076 {
                break;
            }
        }
        ensure!(ob_record.flags() == 163, "could not find OB record");
        let ob_seq = ob_record.seq().as_bytes();
        // Read[5]=A (evidence), Read[17]=G (target) on OB
        assert_eq!(Base::from(ob_seq[5]), A);
        assert_eq!(Base::from(ob_seq[17]), G);

        // OB de-novo CpG context: ref[pos-1]==C → CpG
        {
            let calls = FxHashMap::from_iter([
                (6103081_u32, RastairCall::DeNovoCpg { base: G, methylated: true }),
                (6103093, RastairCall::DeNovoCpg { base: G, methylated: false }),
            ]);
            let ref_bases = FxHashMap::from_iter([
                (6103080_u32, C), // pos-1 for 6103081 → CpG
                (6103081, G),
                (6103092, C), // pos-1 for 6103093 → CpG
                (6103093, G),
            ]);
            let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };
            let mut record = ob_record.clone();
            rewrite_record(&calls, &mut record, BamMode::Legacy, &ref_lookup)?;
            let Aux::String(xm) = record.aux(b"XM")? else { bail!("XM") };
            let annotations: Vec<(usize, char)> =
                xm.chars().enumerate().filter(|(_, c)| *c != '.').collect();
            assert_compact_debug_snapshot!(annotations, @"[(5, 'Z'), (17, 'z')]");
        }

        // OB de-novo CHG context: ref[pos-1]!=C, ref[pos-2]==C → CHG
        {
            let calls = FxHashMap::from_iter([
                (6103081_u32, RastairCall::DeNovoCpg { base: G, methylated: true }),
                (6103093, RastairCall::DeNovoCpg { base: G, methylated: false }),
            ]);
            let ref_bases = FxHashMap::from_iter([
                (6103079_u32, C), // pos-2 for 6103081 → CHG
                (6103080, T),     // pos-1 ≠ C
                (6103081, G),
                (6103091, C), // pos-2 for 6103093 → CHG
                (6103092, T), // pos-1 ≠ C
                (6103093, G),
            ]);
            let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };
            let mut record = ob_record.clone();
            rewrite_record(&calls, &mut record, BamMode::Legacy, &ref_lookup)?;
            let Aux::String(xm) = record.aux(b"XM")? else { bail!("XM") };
            let annotations: Vec<(usize, char)> =
                xm.chars().enumerate().filter(|(_, c)| *c != '.').collect();
            assert_compact_debug_snapshot!(annotations, @"[(5, 'X'), (17, 'x')]");
        }

        // OB de-novo CHH context: ref[pos-1]!=C, ref[pos-2]!=C → CHH
        {
            let calls = FxHashMap::from_iter([
                (6103081_u32, RastairCall::DeNovoCpg { base: G, methylated: true }),
                (6103093, RastairCall::DeNovoCpg { base: G, methylated: false }),
            ]);
            let ref_bases = FxHashMap::from_iter([
                (6103079_u32, A), // pos-2 ≠ C
                (6103080, T),     // pos-1 ≠ C
                (6103081, G),
                (6103091, A), // pos-2 ≠ C
                (6103092, T), // pos-1 ≠ C
                (6103093, G),
            ]);
            let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };
            let mut record = ob_record.clone();
            rewrite_record(&calls, &mut record, BamMode::Legacy, &ref_lookup)?;
            let Aux::String(xm) = record.aux(b"XM")? else { bail!("XM") };
            let annotations: Vec<(usize, char)> =
                xm.chars().enumerate().filter(|(_, c)| *c != '.').collect();
            assert_compact_debug_snapshot!(annotations, @"[(5, 'H'), (17, 'h')]");
        }

        // === OT strand (flag-99 read at 6103078) ===
        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
        bam.fetch(FetchDefinition::All)?;
        let mut ot_record = Record::new();
        while let Some(result) = bam.read(&mut ot_record) {
            result?;
            if ot_record.flags() == 99 && ot_record.pos() == 6103078 {
                break;
            }
        }
        ensure!(ot_record.flags() == 99, "could not find OT record");
        let ot_seq = ot_record.seq().as_bytes();
        // Read[2]=T (evidence), Read[14]=T (evidence) on OT
        assert_eq!(Base::from(ot_seq[2]), T);
        assert_eq!(Base::from(ot_seq[14]), T);

        // OT de-novo CpG context: ref[pos+1]==G → CpG
        {
            let calls = FxHashMap::from_iter([
                (6103080_u32, RastairCall::DeNovoCpg { base: C, methylated: true }),
                (6103092, RastairCall::DeNovoCpg { base: C, methylated: true }),
            ]);
            let ref_bases = FxHashMap::from_iter([
                (6103080_u32, C),
                (6103081, G), // pos+1 → CpG
                (6103092, C),
                (6103093, G), // pos+1 → CpG
            ]);
            let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };
            let mut record = ot_record.clone();
            rewrite_record(&calls, &mut record, BamMode::Legacy, &ref_lookup)?;
            let Aux::String(xm) = record.aux(b"XM")? else { bail!("XM") };
            let annotations: Vec<(usize, char)> =
                xm.chars().enumerate().filter(|(_, c)| *c != '.').collect();
            assert_compact_debug_snapshot!(annotations, @"[(2, 'Z'), (14, 'Z')]");
        }

        // OT de-novo CHG context: ref[pos+1]!=G, ref[pos+2]==G → CHG
        {
            let calls = FxHashMap::from_iter([
                (6103080_u32, RastairCall::DeNovoCpg { base: C, methylated: true }),
                (6103092, RastairCall::DeNovoCpg { base: C, methylated: true }),
            ]);
            let ref_bases = FxHashMap::from_iter([
                (6103080_u32, C),
                (6103081, A), // pos+1 ≠ G
                (6103082, G), // pos+2 → CHG
                (6103092, C),
                (6103093, A), // pos+1 ≠ G
                (6103094, G), // pos+2 → CHG
            ]);
            let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };
            let mut record = ot_record.clone();
            rewrite_record(&calls, &mut record, BamMode::Legacy, &ref_lookup)?;
            let Aux::String(xm) = record.aux(b"XM")? else { bail!("XM") };
            let annotations: Vec<(usize, char)> =
                xm.chars().enumerate().filter(|(_, c)| *c != '.').collect();
            assert_compact_debug_snapshot!(annotations, @"[(2, 'X'), (14, 'X')]");
        }

        // OT de-novo CHH context: ref[pos+1]!=G, ref[pos+2]!=G → CHH
        {
            let calls = FxHashMap::from_iter([
                (6103080_u32, RastairCall::DeNovoCpg { base: C, methylated: true }),
                (6103092, RastairCall::DeNovoCpg { base: C, methylated: true }),
            ]);
            let ref_bases = FxHashMap::from_iter([
                (6103080_u32, C),
                (6103081, A), // pos+1 ≠ G
                (6103082, A), // pos+2 ≠ G → CHH
                (6103092, C),
                (6103093, A), // pos+1 ≠ G
                (6103094, A), // pos+2 ≠ G → CHH
            ]);
            let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };
            let mut record = ot_record.clone();
            rewrite_record(&calls, &mut record, BamMode::Legacy, &ref_lookup)?;
            let Aux::String(xm) = record.aux(b"XM")? else { bail!("XM") };
            let annotations: Vec<(usize, char)> =
                xm.chars().enumerate().filter(|(_, c)| *c != '.').collect();
            assert_compact_debug_snapshot!(annotations, @"[(2, 'H'), (14, 'H')]");
        }

        Ok(())
    }

    /// Reads with indels (insertions/deletions) must have correct position
    /// mapping through aligned_pairs_full. Verify a read with a 1bp deletion
    /// (61M1D19M) produces correct XM and MM output.
    #[test]
    fn golden_indel_read_legacy_and_standard_agree() -> Result<()> {
        use Base::*;

        // Flag 99 read at 6106220 (BAM 0-based) with CIGAR 61M1D19M
        let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
        bam.fetch(FetchDefinition::All)?;
        let mut record = Record::new();
        while let Some(result) = bam.read(&mut record) {
            result?;
            if record.flags() == 99 && record.pos() == 6106220 {
                break;
            }
        }
        ensure!(record.flags() == 99 && record.pos() == 6106220, "could not find indel read");

        let strand = StrandFromRecord::strand(&record);
        assert_eq!(strand, Strand::OT);
        let seq = record.seq().as_bytes();

        // CIGAR: 61M1D19M
        // Read[0-60] → ref[6106220-6106280] (61M)
        // ref[6106281] deleted (1D)
        // Read[61-79] → ref[6106282-6106300] (19M)
        //
        // CpG calls:
        // Ref pos 6106229 (read offset 9): OT C side
        //   Read[9] = T → methylation evidence → methylated
        // Ref pos 6106300 (read offset 79): OT C side, AFTER the deletion
        //   Read[79] = C → target base → unmethylated
        assert_eq!(Base::from(seq[9]), T, "read[9] should be T (methylation evidence)");
        assert_eq!(Base::from(seq[79]), C, "read[79] should be C (unmethylated, after deletion)");

        let calls = FxHashMap::from_iter([
            (6106229_u32, RastairCall::Cpg { base: C, methylated: true }),
            (6106300, RastairCall::Cpg { base: C, methylated: false }),
        ]);
        let ref_bases =
            FxHashMap::from_iter([(6106229_u32, C), (6106230, G), (6106300, C), (6106301, G)]);
        let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };

        // === Legacy mode ===
        let mut legacy_record = record.clone();
        rewrite_record(&calls, &mut legacy_record, BamMode::Legacy, &ref_lookup)?;

        let Aux::String(xm_tag) = legacy_record.aux(b"XM")? else { bail!("XM") };
        assert_eq!(xm_tag.len(), seq.len());
        let xm_annotations: Vec<(usize, char)> =
            xm_tag.chars().enumerate().filter(|(_, c)| *c != '.').collect();
        // Position 9 = Z (methylated), position 79 = z (unmethylated, AFTER deletion)
        assert_compact_debug_snapshot!(xm_annotations, @"[(9, 'Z'), (79, 'z')]");

        // === Standard mode ===
        let mut standard_record = record.clone();
        rewrite_record(&calls, &mut standard_record, BamMode::Standard, &ref_lookup)?;

        let Aux::String(mm_tag) = standard_record.aux(b"MM")? else { bail!("MM") };

        // Rewritten SEQ: T→C at position 9
        let new_seq = standard_record.seq().as_bytes();
        assert_eq!(Base::from(new_seq[9]), C, "methylated T should be rewritten to C");
        assert_eq!(Base::from(new_seq[79]), C, "unmethylated C should stay C");

        // === Cross-check ===
        let xm_methylated = decode_xm_to_positions(xm_tag);
        let fwd_seq = seq_for_mm_tag(&standard_record);
        let (mm_base, mm_methylated) = decode_mm_to_positions(mm_tag, &fwd_seq)?;
        assert_eq!(mm_base, C);
        // Only position 9 is methylated
        assert_eq!(xm_methylated, vec![9]);
        assert_eq!(mm_methylated, vec![9]);

        Ok(())
    }

    /// All four standard flag combinations (99, 147, 83, 163) with the same
    /// genomic CpG must produce consistent annotations. OT reads annotate the
    /// C side, OB reads annotate the G side.
    #[test]
    fn all_four_flags_consistent_annotations() -> Result<()> {
        use Base::*;

        // The CpG at ref 6103080(C)/6103081(G) is covered by all four flags
        // in the test region.
        let ref_bases = FxHashMap::from_iter([(6103080_u32, C), (6103081, G)]);
        let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };

        let test_cases: [(u16, Strand, &str, &str); 4] = [
            (99, Strand::OT, "CT", "CT"),
            (147, Strand::OT, "GA", "CT"),
            (83, Strand::OB, "CT", "GA"),
            (163, Strand::OB, "GA", "GA"),
        ];

        for (flag, expected_strand, expected_xr, expected_xg) in test_cases {
            let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
            bam.fetch(FetchDefinition::All)?;
            let mut record = Record::new();
            let mut found = false;
            while let Some(result) = bam.read(&mut record) {
                result?;
                if record.flags() == flag {
                    // Make sure this read covers our CpG position
                    let start = record.pos() as u32;
                    let end = start + record.seq_len() as u32;
                    if start <= 6103081 && end > 6103081 {
                        found = true;
                        break;
                    }
                }
            }
            ensure!(found, "No record with flag {flag} covering CpG at 6103080/81");

            let strand = StrandFromRecord::strand(&record);
            assert_eq!(strand, expected_strand, "strand mismatch for flag {flag}");

            // Build calls appropriate for this strand
            let calls = if strand == Strand::OT {
                FxHashMap::from_iter([(
                    6103080_u32,
                    RastairCall::Cpg { base: C, methylated: true },
                )])
            } else {
                FxHashMap::from_iter([(
                    6103081_u32,
                    RastairCall::Cpg { base: G, methylated: true },
                )])
            };

            // === Legacy ===
            let mut legacy_record = record.clone();
            rewrite_record(&calls, &mut legacy_record, BamMode::Legacy, &ref_lookup)?;

            let Aux::String(xr_tag) = legacy_record.aux(b"XR")? else { bail!("XR for {flag}") };
            let Aux::String(xg_tag) = legacy_record.aux(b"XG")? else { bail!("XG for {flag}") };
            let Aux::String(xm_tag) = legacy_record.aux(b"XM")? else { bail!("XM for {flag}") };

            assert_eq!(xr_tag, expected_xr, "XR mismatch for flag {flag}");
            assert_eq!(xg_tag, expected_xg, "XG mismatch for flag {flag}");
            assert_eq!(xm_tag.len(), record.seq_len(), "XM length mismatch for flag {flag}");

            // Every flag should have at least one annotation (Z or z)
            let has_annotation = xm_tag.chars().any(|c| c != '.');
            assert!(has_annotation, "XM has no annotations for flag {flag}: {xm_tag}");

            // === Standard ===
            let mut standard_record = record.clone();
            rewrite_record(&calls, &mut standard_record, BamMode::Standard, &ref_lookup)?;

            // Cross-check: both modes find the same number of methylated positions.
            // MM/ML are absent when no methylation evidence is found (absent = 0 positions).
            let xm_methylated = decode_xm_to_positions(xm_tag);
            let mm_methylated_count = match standard_record.aux(b"MM") {
                Ok(Aux::String(mm_tag)) => {
                    let fwd_seq = seq_for_mm_tag(&standard_record);
                    let (_, positions) = decode_mm_to_positions(mm_tag, &fwd_seq)?;
                    positions.len()
                }
                Ok(_) => bail!("MM not a string for flag {flag}"),
                Err(_) => 0,
            };
            assert_eq!(
                xm_methylated.len(),
                mm_methylated_count,
                "Methylated position count mismatch for flag {flag}: \
                 XM={xm_methylated:?}, MM_count={mm_methylated_count}"
            );
        }

        Ok(())
    }

    #[test]
    fn test_reverse_complement() {
        // Simple test
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"CCCC"), b"GGGG");

        // More complex sequence
        assert_eq!(reverse_complement(b"ATCGATCG"), b"CGATCGAT");

        // With unknown bases
        assert_eq!(reverse_complement(b"ANTN"), b"NANT");

        // Empty sequence
        assert_eq!(reverse_complement(b""), b"");

        // Real-world example: forward sequence with methylation evidence
        let forward = b"ATTGT"; // Original with T at position 2 (methylated C)
        let reverse = reverse_complement(forward);
        assert_eq!(reverse, b"ACAAT");
    }
}
