use crate::{
    bed::reader::{RastairBedReader, RastairCall},
    sequence::{ChunkRegion, ReaderParams, Region},
    utils::{cli, logging::ThisIsABug},
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
use rustc_hash::{FxHashMap, FxHashSet};
use std::thread::available_parallelism;

mod base_modification;
pub use base_modification::{MethylatedPositions, XrTags};
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

            for records in bam_receiver {
                for record in records {
                    writer.write(&record).wrap_err("failed to write record to new BAM file")?;
                }
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
    }

    let records = BAM_READER.with(|bam_cell| {
        BED_READER.with(|bed_cell| -> Result<Vec<Record>> {
            let mut bam_opt = bam_cell.borrow_mut();
            let mut bed_opt = bed_cell.borrow_mut();

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

            let bam = bam_opt
                .as_mut()
                .wrap_err("thread-local BAM reader not initialized")
                .this_is_a_bug()?;
            let bed = bed_opt
                .as_mut()
                .wrap_err("thread-local BED reader not initialized")
                .this_is_a_bug()?;

            rewrite_region(bam, bed, &segment.region, mode, is_last)
        })
    })?;

    if let Err(err) = sender.send(index, records) {
        trace!(error = format!("{err:#}"), "Failed to send BAM records, channel probably closed");
    }

    Ok(())
}

#[instrument(level = "debug", skip_all, fields(region = %region))]
fn rewrite_region(
    bam: &mut bam::IndexedReader,
    calls_reader: &mut RastairBedReader,
    region: &Region,
    mode: BamMode,
    is_last_segment: bool,
) -> Result<Vec<Record>> {
    FetchDefinition::try_from(region)
        .wrap_err("Could not convert region string")
        .and_then(|r| bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))
        .wrap_err_with(|| format!("Could not fetch region `{}` from BAM file", region))?;

    let noodle_region = region
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

        rewrite_record(&calls, &mut record, mode).wrap_err("failed to rewrite record")?;
        out.push(record.clone());
    }

    Ok(out)
}

#[instrument(level = "debug", skip_all, fields(pos = record.pos()))]
fn rewrite_record(
    calls: &FxHashMap<u32, RastairCall>,
    record: &mut Record,
    mode: BamMode,
) -> Result<()> {
    let strand = StrandFromRecord::strand(record);
    let is_first_in_pair = record.is_first_in_template();

    let is_reverse = record.is_reverse();

    match mode {
        BamMode::Standard => {
            // Methylation detection works in stored (+ strand) orientation:
            // T→C for OT, A→G for OB. Positions are stored-SEQ indices.
            let MethylatedInfo { seq, methylated_positions } =
                get_methylated_positions(calls, record);

            // The MM tag spec requires positions relative to the original read
            // (5' to 3'). For forward reads this equals the stored SEQ. For
            // reverse reads the stored SEQ is the reverse complement of the
            // original read, so we must convert both the sequence and positions
            // to original-read orientation, and flip the base/strand qualifier.
            let (mm_seq, mm_positions, mm_strand) = if is_reverse {
                let seq_len = seq.len();
                let original_positions: SmallVec<u32, 10> = methylated_positions
                    .iter()
                    .map(|&p| u32::try_from(seq_len - 1 - p as usize).expect("position fits"))
                    .collect();
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
            let cpg = find_legacy_cpg_positions(calls, record);
            let seq = record.seq().as_bytes();

            let xr_tags =
                XrTags::new_legacy(&seq, strand, is_first_in_pair, &cpg.methylated, &cpg.all_cpg);
            xr_tags.apply_to_record(record)?;
        }
    }

    Ok(())
}

struct LegacyCpgPositions {
    /// Read positions where a methylated CpG was observed
    methylated: FxHashSet<usize>,
    /// All read positions that overlap a CpG call (methylated or not)
    all_cpg: FxHashSet<usize>,
}

/// Find CpG positions for legacy XM tag generation.
///
/// Returns both the set of methylated positions and the set of all CpG
/// positions (methylated + unmethylated) so the XM generator can distinguish
/// "CpG but unmethylated" (`z`) from "not a CpG at all" (`.`).
fn find_legacy_cpg_positions(
    calls: &FxHashMap<u32, RastairCall>,
    record: &Record,
) -> LegacyCpgPositions {
    let strand = StrandFromRecord::strand(record);
    let seq = record.seq().as_bytes();
    let mut methylated = FxHashSet::default();
    let mut all_cpg = FxHashSet::default();

    let target_base = match strand {
        Strand::OT => Base::C,
        Strand::OB => Base::G,
        Strand::Unknown => return LegacyCpgPositions { methylated, all_cpg },
    };
    let evidence_base = match strand {
        Strand::OT => Base::T,
        Strand::OB => Base::A,
        Strand::Unknown => unreachable!(),
    };

    for [pos_in_read, pos_in_ref] in record.aligned_pairs_full() {
        let Some(pos_in_read) = pos_in_read else { continue };
        let Some(pos_in_ref) = pos_in_ref else { continue };
        let pos_in_ref = u32::try_from(pos_in_ref).expect("position fits in u32");
        let pos_in_read = pos_in_read as usize;

        let is_methylated = match calls.get(&pos_in_ref) {
            Some(
                RastairCall::Cpg { methylated, .. } | RastairCall::DeNovoCpg { methylated, .. },
            ) => *methylated,
            _ => continue,
        };

        let observed_base = Base::from(seq[pos_in_read]);

        // A CpG position in the read shows either the target base (unmethylated)
        // or the evidence base (methylated, converted by TAPS)
        if observed_base == target_base || observed_base == evidence_base {
            all_cpg.insert(pos_in_read);
            if is_methylated && observed_base == evidence_base {
                methylated.insert(pos_in_read);
            }
        }
    }

    LegacyCpgPositions { methylated, all_cpg }
}

struct MethylatedInfo {
    seq: Vec<u8>,
    methylated_positions: SmallVec<u32, 10>,
}

fn get_methylated_positions(
    calls: &FxHashMap<u32, RastairCall>,
    record: &Record,
) -> MethylatedInfo {
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
        let pos_in_ref = u32::try_from(pos_in_ref).expect("position fits in u32");
        let pos_in_read = pos_in_read as usize;

        let methylated = match calls.get(&pos_in_ref) {
            Some(RastairCall::Cpg { methylated, .. } | RastairCall::DeNovoCpg { methylated, .. }) => {
                *methylated
            }
            _ => continue,
        };

        if methylated {
            let observed_base = Base::from(seq[pos_in_read]);

            match strand {
                Strand::OT => {
                    if observed_base == T {
                        seq[pos_in_read] = *C;
                        methylated_positions.push(pos_in_read as u32);
                    }
                }
                Strand::OB => {
                    if observed_base == A {
                        seq[pos_in_read] = *G;
                        methylated_positions.push(pos_in_read as u32);
                    }
                }
                Strand::Unknown => continue,
            }
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

        let data = get_methylated_positions(&calls, &record);
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

        rewrite_record(&calls, &mut record, BamMode::Legacy)?;

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
            rewrite_record(&calls, &mut record, BamMode::Legacy)?;

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
            rewrite_record(&calls, &mut record, BamMode::Standard)?;

            let Aux::String(mm_tag) = record.aux(b"MM").wrap_err("missing MM tag")? else {
                bail!("MM tag is not a string");
            };

            let fwd_seq = seq_for_mm_tag(&record);
            let (_, mm_positions) =
                decode_mm_to_positions(mm_tag, &fwd_seq).wrap_err_with(|| {
                    format!("decode MM {:?} pos={} flag={}", mm_tag, record.pos(), record.flags())
                })?;

            // Standard mode should NOT produce XR/XG/XM tags
            assert!(record.aux(b"XR").is_err(), "Standard mode should not produce XR tag");

            // MM should produce valid positions
            assert!(mm_positions.len() <= fwd_seq.len(), "MM positions out of range");

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
            rewrite_record(&calls, &mut record, BamMode::Standard)?;
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

            let Aux::String(mm_tag) = record.aux(b"MM").wrap_err("missing MM tag")? else {
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
            rewrite_record(&calls, &mut record, BamMode::Legacy)?;
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
        rewrite_record(&calls, &mut record, BamMode::Standard)?;
        writer.write(&record)?;
        drop(writer);

        let mut reader = bam::Reader::from_path(&temp_bam)?;
        reader.read(&mut record).wrap_err("no records")?.wrap_err("read back")?;

        let Aux::String(mm_tag) = record.aux(b"MM").wrap_err("missing MM")? else {
            bail!("MM not a string");
        };

        let fwd_seq = seq_for_mm_tag(&record);
        let (_, mm_positions) = decode_mm_to_positions(mm_tag, &fwd_seq)?;
        assert!(mm_positions.is_empty(), "expected no MM methylation");

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
        rewrite_record(&calls, &mut record, BamMode::Legacy)?;
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
                calls.insert(
                    pos_in_ref as u32 + 1,
                    RastairCall::Cpg { methylated: true, base: G },
                );
                cpg_count += 1;
                break;
            }
            // For OT: look for TG (methylation evidence for C)
            if strand == Strand::OT && current_base == T && next_base == G && cpg_count == 0 {
                calls.insert(
                    pos_in_ref as u32,
                    RastairCall::Cpg { methylated: true, base: C },
                );
                cpg_count += 1;
                break;
            }
        }
        ensure!(cpg_count == 1, "need at least one CpG call for this test");

        rewrite_record(&calls, &mut record, BamMode::Legacy)?;

        let Aux::String(xm_tag) = record.aux(b"XM")? else {
            bail!("XM not a string");
        };

        let z_count = xm_tag.chars().filter(|c| *c == 'Z' || *c == 'z').count();

        // BUG: currently every C or G in the read is marked as z/Z, but only
        // the one CpG position we called should be annotated.
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
                calls.insert(
                    pos_in_ref as u32,
                    RastairCall::DeNovoCpg { methylated: true, base },
                );
                found = true;
                break;
            }
        }
        ensure!(found, "need a methylation-evidence base for this test");

        rewrite_record(&calls, &mut record, BamMode::Legacy)?;

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
                calls.insert(
                    pos_in_ref as u32,
                    RastairCall::DeNovoCpg { methylated: true, base },
                );
                expected_pos = Some(pos_in_read as u32);
                break;
            }
        }
        ensure!(expected_pos.is_some(), "need a methylation-evidence base");

        let data = get_methylated_positions(&calls, &record);
        assert!(
            !data.methylated_positions.is_empty(),
            "get_methylated_positions should find DeNovoCpg calls, but found none"
        );

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
