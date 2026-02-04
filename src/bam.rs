use crate::{
    bed::reader::{RastairBedReader, RastairCall},
    sequence::{ChunkRegion, ReaderParams, Region},
    utils::{cli, logging::ThisIsABug},
};
use clap::{Parser, value_parser};
use clio::ClioPath;
use color_eyre::eyre::{Context, Result};
use rastair_types::SmallVec;
use rastair_types::{Base, Strand, StrandFromRecord};
use rust_htslib::bam::{
    self, FetchDefinition, Header, Read, Record, Writer, ext::BamRecordExtensions as _,
    header::HeaderRecord,
};
use rustc_hash::{FxHashMap, FxHashSet};

mod base_modification;
pub use base_modification::{MethylatedPositions, XrTags};
use tracing::{instrument, warn};

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
))]
pub fn rewrite(params: &BamRewriteArgs) -> Result<()> {
    let mut readers = params.segments.readers().wrap_err("Failed to read BAM/FASTA files")?;
    let regions: Vec<ChunkRegion> = readers
        .segments(params.segment_max_length, 0)
        .wrap_err("Could not fetch segments from BAM file")?
        .collect();

    let mut calls_reader =
        RastairBedReader::new(params.calls_file.path()).wrap_err("failed to open calls file")?;

    let output_file = &params.output;
    let mut writer = {
        let header = {
            let mut header = Header::from_template(readers.bam.header());
            add_rastair_header(&mut header);
            header
        };

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
        writer
    };

    for segment in &regions {
        rewrite_region(&mut readers.bam, &mut calls_reader, &mut writer, &segment.region)
            .wrap_err_with(|| format!("Failed to rewrite region {}", segment.region))?;
    }

    Ok(())
}

#[instrument(level = "debug", skip_all, fields(region = %region))]
fn rewrite_region(
    bam: &mut bam::IndexedReader,
    calls_reader: &mut RastairBedReader,
    writer: &mut rust_htslib::bam::Writer,
    region: &Region,
) -> Result<()> {
    FetchDefinition::try_from(region)
        .wrap_err("Could not convert region string")
        .and_then(|r| bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))
        .wrap_err_with(|| format!("Could not fetch region `{}` from BAM file", region))?;

    let noodle_region = region
        .clone()
        .try_into()
        .wrap_err("Failed to convert region representation")
        .this_is_a_bug()?;

    // Load all calls in this region into a map of position -> call
    // Since we're going region by region, all calls are from the same contig
    let calls: FxHashMap<u32, RastairCall> = calls_reader
        .query(&noodle_region)
        .wrap_err("failed to query calls file")?
        .iter()
        .map(|call| (call.pos, call.call.clone()))
        .collect();

    let mut record = Record::new();
    while let Some(result) = bam.read(&mut record) {
        if let Err(error) = result {
            warn!(%error, "Failed to read BAM record");
        }

        rewrite_record(&calls, &mut record).wrap_err("failed to rewrite record")?;

        writer.write(&record).wrap_err("failed to write record to new BAM file")?;
    }

    Ok(())
}

/// Rewrite a single BAM record with methylation information from calls
///
/// 1. Find all calls that are covered by our record
/// 2. Find the positions in the read that correspond to the called positions
/// 3. Build a `MethylatedPositions` struct for MM/ML tags
/// 4. Build `XrTags` for XR/XG/XM tags
/// 5. Rewrite the sequence to un-modify the bases of methylated positions
///    (i.e., change T or A that is methylation evidence back to C or G)
/// 6. Apply the modifications to the record
#[instrument(level = "debug", skip_all, fields(pos = record.pos()))]
fn rewrite_record(
    // All calls in the region of the record (they are all the same contig)
    calls: &FxHashMap<u32, RastairCall>,
    record: &mut Record,
) -> Result<()> {
    let is_reverse = record.is_reverse();
    let MethylatedInfo { seq, methylated_positions } =
        get_methylated_positions(calls, record, is_reverse);

    let strand = StrandFromRecord::strand(record);

    // For reverse reads, we need to work with the original sequence for MM tag generation
    // but keep the SEQ in reverse complement form
    let seq_for_mm_tag = if is_reverse { reverse_complement(&seq) } else { seq.clone() };

    let mods = MethylatedPositions::new(strand, &seq_for_mm_tag, &methylated_positions);

    // Determine if this read is first or second in pair for XR tag
    let is_first_in_pair = record.is_first_in_template();

    // Create XR/XG/XM tags (reusing a HashSet for performance)
    let mut methylated_set = FxHashSet::default();
    let xr_tags = XrTags::new(
        &seq_for_mm_tag,
        strand,
        &methylated_positions,
        is_first_in_pair,
        &mut methylated_set,
    );

    record.set_seq(&seq);

    mods.apply_to_record(record)?;
    xr_tags.apply_to_record(record)?;
    Ok(())
}

struct MethylatedInfo {
    seq: Vec<u8>,
    methylated_positions: SmallVec<u32, 10>,
}

fn get_methylated_positions(
    // All calls in the region of the record (they are all the same contig)
    calls: &FxHashMap<u32, RastairCall>,
    record: &Record,
    is_reverse: bool,
) -> MethylatedInfo {
    use Base::*;

    let strand = StrandFromRecord::strand(record);

    // For reverse reads, we need to work with the original sequence
    let current_seq = record.seq().as_bytes();
    let mut seq = if is_reverse { reverse_complement(&current_seq) } else { current_seq };

    let seq_len = seq.len();
    let mut methylated_positions: SmallVec<u32, 10> = SmallVec::new();

    for [pos_in_read, pos_in_ref] in record.aligned_pairs_full() {
        let Some(pos_in_read) = pos_in_read else {
            continue;
        };
        let Some(pos_in_ref) = pos_in_ref else {
            continue;
        };
        let pos_in_ref = u32::try_from(pos_in_ref).expect("position fits in u32");
        let pos_in_read = u32::try_from(pos_in_read).expect("position fits in u32");

        if let Some(call) = calls.get(&pos_in_ref)
            && let RastairCall::Cpg { methylated, .. } = call
            && *methylated
        {
            // For reverse reads, map current position to original position
            let original_pos = if is_reverse {
                u32::try_from(seq_len - 1 - pos_in_read as usize).expect("position fits in u32")
            } else {
                pos_in_read
            };

            let pos = original_pos as usize;
            let observed_base = Base::from(seq[pos]);

            // Let's rewrite the sequence to un-modify the bases and collect
            // methylated positions
            match strand {
                // If we are on the top strand, we expect to see T at methylated
                // C positions
                Strand::OT => {
                    if observed_base == T {
                        seq[pos] = *C;
                        methylated_positions.push(original_pos);
                    }
                }
                // If we are on the bottom strand, we expect to see A at
                // methylated G
                Strand::OB => {
                    if observed_base == A {
                        seq[pos] = *G;
                        methylated_positions.push(original_pos);
                    }
                }
                Strand::Unknown => continue,
            }
        }
    }

    // For reverse reads, reverse complement back to get the modified SEQ
    if is_reverse {
        seq = reverse_complement(&seq);
    }

    MethylatedInfo { seq, methylated_positions }
}

/// Reverse complement a DNA sequence
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

        let data = get_methylated_positions(&calls, &record, false);
        let new_seq = &data.seq;

        // for human comparison
        assert_snapshot!(as_base_string(&ol_seq), @"AAAGGCGTGCACCACCACGCCTGGCTTGGTTTGGTTTTTGATTGGTTGGTTGGTCTTTTGAGACAGGGTTTCTCTGTGTA");
        assert_snapshot!(as_base_string(new_seq), @"AAAGGCGTGCGCCGCCGCGCCTGGCTTGGTTTGGTTTTTGATTGGTTGGTTGGTCTTTTGAGACGGGGTTTCTCTGTGTA");
        // changed positions:                          __ _ _ ^  ^  ^ _   __   __   __     _   __  __  __      _ _  ^___       _ _
        // positions in list of all Gs:                       4  0  0                                               16

        assert_compact_debug_snapshot!(data.methylated_positions, @"[10, 13, 16, 64]");
        // okay but are these really all As that were methylated Gs?
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

        let methylated_positions =
            MethylatedPositions::new(Strand::OB, new_seq, &data.methylated_positions);

        assert_compact_debug_snapshot!(methylated_positions, @"MethylatedPositions { base: G, strand: OB, positions: [4, 0, 0, 16] }");
        methylated_positions
            .apply_to_record(&mut record)
            .wrap_err("failed to apply modifications to record")?;

        let Aux::String(mod_string) = record.aux(b"MM").wrap_err("missing MM tag")? else {
            bail!("MM tag is not a string");
        };
        assert_snapshot!(mod_string, @"G-m,4,0,0,16;");

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

        rewrite_record(&calls, &mut record)?;

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
        // Test all four flag combinations
        let test_cases = vec![
            (99, "CT", "CT"),  // flag 99: OT, first in pair => XR:CT, XG:CT
            (147, "GA", "CT"), // flag 147: OT, second in pair => XR:GA, XG:CT
            (83, "CT", "GA"),  // flag 83: OB, first in pair => XR:CT, XG:GA
            (163, "GA", "GA"), // flag 163: OB, second in pair => XR:GA, XG:GA
        ];

        for (flag, expected_xr, expected_xg) in test_cases {
            let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")
                .wrap_err("failed to open test BAM")?;

            // Find a record with the desired flag
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

            let calls = FxHashMap::default(); // Empty calls for this test
            rewrite_record(&calls, &mut record)?;

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

    fn as_base_string(seq: &[u8]) -> String {
        seq.iter().map(|b| Base::from(*b).as_str()).collect()
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
