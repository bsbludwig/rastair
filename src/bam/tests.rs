use color_eyre::eyre::{bail, ensure};
use insta::{assert_compact_debug_snapshot, assert_snapshot};
use seqair::bam::aux::{AuxValue, find_tag};
use seqair_types::BamFlags;

use super::*;

const TEST_BAM: &str = "tests/data/test.bam";

fn open_test_bam() -> Result<AlignmentStream> {
    Ok(AlignmentStream::from_reader(std::fs::File::open(TEST_BAM)?)?)
}

/// The first record of the test BAM, in file order, that `pred` accepts.
fn find_record(pred: impl Fn(&RawBamRecord) -> bool) -> Result<RawBamRecord> {
    let mut stream = open_test_bam()?;
    let mut record = RawBamRecord::default();
    while stream.read_record(&mut record)? {
        if pred(&record) {
            return Ok(record);
        }
    }
    bail!("no record in the test BAM matches")
}

/// The records on the first contig that overlap `[start, end)`, in file order.
fn records_overlapping(start: u32, end: u32) -> Result<Vec<RawBamRecord>> {
    let mut stream = open_test_bam()?;
    let mut record = RawBamRecord::default();
    let mut found = Vec::new();
    while stream.read_record(&mut record)? {
        let (Some(tid), Some(pos), Some(last)) = (record.tid(), record.pos(), record.end_pos())
        else {
            continue;
        };
        // The file is sorted, so nothing after `end` on this contig can overlap.
        if tid.as_u32() > 0 || pos.as_u32() >= end {
            break;
        }
        if last.as_u32() >= start {
            found.push(record.clone());
        }
    }
    Ok(found)
}

/// What `fetch((0, start, end))` followed by one `read` returned.
fn first_overlapping(start: u32, end: u32) -> Result<RawBamRecord> {
    records_overlapping(start, end)?.into_iter().next().wrap_err("no records")
}

fn find_flag_at(flag_value: u16, pos_value: u32) -> Result<RawBamRecord> {
    find_record(|r| flag(r) == flag_value && pos(r) == pos_value)
        .wrap_err_with(|| format!("could not find flag-{flag_value} record at {pos_value}"))
}

fn flag(record: &RawBamRecord) -> u16 {
    record.flags().raw()
}

fn pos(record: &RawBamRecord) -> u32 {
    record.pos().map_or(0, |pos| pos.as_u32())
}

fn stored_seq(record: &RawBamRecord) -> Vec<u8> {
    let mut seq = Vec::new();
    record.decode_seq_into(&mut seq);
    seq
}

/// A string tag; missing or of another type, an error.
fn string_tag(record: &RawBamRecord, tag: &[u8; 2]) -> Result<String> {
    match find_tag(record.aux().as_bytes(), *tag) {
        None => bail!("missing {} tag", tag.escape_ascii()),
        Some(AuxValue::String(value)) => Ok(String::from_utf8(value.to_vec())?),
        Some(other) => bail!("{} is not a string tag: {other:?}", tag.escape_ascii()),
    }
}

fn has_tag(record: &RawBamRecord, tag: &[u8; 2]) -> bool {
    find_tag(record.aux().as_bytes(), *tag).is_some()
}

/// The read bases that sit on the reference, as read index and reference position.
fn matched(record: &RawBamRecord) -> Vec<(usize, u32)> {
    let cigar = cigar_of(record);
    aligned_bases(record.pos(), &cigar).map(|m| (m.qpos.as_usize(), m.rpos.as_u32())).collect()
}

fn write_bam(path: &std::path::Path, records: &[RawBamRecord]) -> Result<()> {
    let header = BamHeader::from_template(open_test_bam()?.header());
    let mut writer = BamWriterBuilder::to_path(path, &header).build().wrap_err("create writer")?;
    for record in records {
        writer.write_raw(record)?;
    }
    writer.finish()?;
    Ok(())
}

fn read_bam(path: &std::path::Path) -> Result<Vec<RawBamRecord>> {
    let mut stream = AlignmentStream::from_reader(std::fs::File::open(path)?)
        .wrap_err("open roundtrip BAM for reading")?;
    let mut records = Vec::new();
    let mut record = RawBamRecord::default();
    while stream.read_record(&mut record)? {
        records.push(record.clone());
    }
    Ok(records)
}

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
    let record = find_flag_at(163, 6103076)?;

    let strand = Strand::from_flags(record.flags()).wrap_err("unknown strand")?;
    assert_eq!(strand, Strand::OB, "expected OB strand");
    assert!(!record.flags().is_reverse(), "expected forward read");
    let seq = stored_seq(&record);
    assert_snapshot!(as_base_string(&seq), @"AAGGCATGCACCACCACGCCTGGCTTGGTTTGGTTTTTGATTGGTTGGTTGGTCTTTTGAGACAGGGTTTCTCTGTGT");

    // Construct calls for two CpG positions (G side, matching OB strand):
    //
    // Ref pos 6103081 (read offset 5): methylated CpG
    //   Read[5] = A → methylation evidence on OB (A at G ref position)
    //
    // Ref pos 6103093 (read offset 17): unmethylated CpG
    //   Read[17] = G → target base, no methylation evidence
    let calls = FxHashMap::from_iter([
        (6103081, RastairCall::Cpg { base: G, methylated: Some(true) }),
        (6103093, RastairCall::Cpg { base: G, methylated: Some(false) }),
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
    rewrite_record(&calls, &mut legacy_record, BamMode::Legacy, ref_lookup)?;

    let xm_tag = string_tag(&legacy_record, b"XM")?;
    let xr_tag = string_tag(&legacy_record, b"XR")?;
    let xg_tag = string_tag(&legacy_record, b"XG")?;

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
    assert_eq!(stored_seq(&legacy_record), seq, "Legacy mode should not modify SEQ");

    // === Standard mode ===
    let mut standard_record = record.clone();
    rewrite_record(&calls, &mut standard_record, BamMode::Standard, ref_lookup)?;

    let mm_tag = string_tag(&standard_record, b"MM")?;

    // OB strand, not reversed → MM uses G-m format with stored-SEQ positions
    // Rewritten SEQ: position 5 changed from A→G (un-modify methylation evidence)
    let new_seq = stored_seq(&standard_record);
    assert_snapshot!(as_base_string(&new_seq), @"AAGGCGTGCACCACCACGCCTGGCTTGGTTTGGTTTTTGATTGGTTGGTTGGTCTTTTGAGACAGGGTTTCTCTGTGT");

    // Verify position 5 was rewritten (A→G) and position 17 unchanged (G stays G)
    assert_eq!(Base::from(seq[5]), A, "original should have A at pos 5");
    assert_eq!(Base::from(new_seq[5]), G, "rewritten should have G at pos 5");
    assert_eq!(Base::from(new_seq[17]), G, "pos 17 should stay G (unmethylated)");

    // MM tag: G-m encoding with skip list
    assert_snapshot!(mm_tag, @"G-m,2;");

    // === Cross-check: both modes agree on methylated positions ===
    let xm_methylated = decode_xm_to_positions(&xm_tag);

    let fwd_seq = seq_for_mm_tag(&standard_record);
    let (mm_base, mm_methylated) = decode_mm_to_positions(&mm_tag, &fwd_seq)?;
    assert_eq!(mm_base, G, "MM should target G for OB strand");

    // Both should identify position 5 as the only methylated position
    assert_eq!(xm_methylated, vec![5], "XM methylated positions");
    assert_eq!(mm_methylated, vec![5], "MM methylated positions");

    // Standard mode should NOT have XR/XG/XM tags
    assert!(!has_tag(&standard_record, b"XR"));
    // Legacy mode should NOT have MM/ML tags
    assert!(!has_tag(&legacy_record, b"MM"));

    Ok(())
}

/// XM and MM/ML tags must reflect per-read methylation evidence, not the
/// position-level call. A CpG with beta ≤ 0.5 (`methylated: Some(false)`) still
/// has individually methylated reads that must show Z (not z) in XM and
/// appear in MM/ML.
#[test]
fn per_read_methylation_independent_of_position_call() -> Result<()> {
    use Base::*;

    // Same OB flag-163 read as golden_legacy_and_standard_agree
    let record = find_flag_at(163, 6103076)?;

    let strand = Strand::from_flags(record.flags()).wrap_err("unknown strand")?;
    assert_eq!(strand, Strand::OB);

    // Position 6103081 (read offset 5): Read[5] = A → evidence base for OB
    // Position 6103093 (read offset 17): Read[17] = G → target base
    //
    // Mark BOTH as `methylated: Some(false)` (simulating beta ≤ 0.5).
    // Despite the position-level call, read offset 5 shows methylation
    // evidence (A at G ref) and must be annotated as Z, not z.
    let calls = FxHashMap::from_iter([
        (6103081, RastairCall::Cpg { base: G, methylated: Some(false) }),
        (6103093, RastairCall::Cpg { base: G, methylated: Some(false) }),
    ]);

    let ref_bases = FxHashMap::from_iter([(6103080, C), (6103081, G), (6103092, C), (6103093, G)]);
    let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };

    // === Legacy mode ===
    let mut legacy_record = record.clone();
    rewrite_record(&calls, &mut legacy_record, BamMode::Legacy, ref_lookup)?;

    let xm_tag = string_tag(&legacy_record, b"XM")?;

    // Offset 5 must be Z (methylated per-read), offset 17 must be z (unmethylated per-read)
    let xm_annotations: Vec<(usize, char)> =
        xm_tag.chars().enumerate().filter(|(_, c)| *c != '.').collect();
    assert_compact_debug_snapshot!(xm_annotations, @"[(5, 'Z'), (17, 'z')]");

    // === Standard mode ===
    let mut standard_record = record.clone();
    rewrite_record(&calls, &mut standard_record, BamMode::Standard, ref_lookup)?;

    let new_seq = stored_seq(&standard_record);
    // Position 5 must be rewritten A→G (undo methylation evidence)
    assert_eq!(Base::from(new_seq[5]), G, "position 5 should be rewritten A→G");

    let mm_tag = string_tag(&standard_record, b"MM")?;
    // MM tag should report the methylation at position 5
    let fwd_seq = seq_for_mm_tag(&standard_record);
    let (mm_base, mm_methylated) = decode_mm_to_positions(&mm_tag, &fwd_seq)?;
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
    let record = first_overlapping(6103075, 6103100)?;

    assert_eq!(flag(&record), 83);
    let strand = Strand::from_flags(record.flags()).wrap_err("unknown strand")?;
    assert_eq!(strand, Strand::OB);
    assert!(record.flags().is_reverse());
    let seq = stored_seq(&record);

    // Two CpG calls (G side for OB):
    // Pos 6103081 (read offset 6): Read[6]=G → target, unmethylated in this read
    // Pos 6103093 (read offset 18): Read[18]=G → target, unmethylated in this read
    //
    // Also add a methylated call where read actually shows evidence:
    // Let's find a position where the read has A (OB evidence)
    let mut evidence_pos = None;
    for (pos_in_read, pos_in_ref) in matched(&record) {
        if Base::from(seq[pos_in_read]) == A {
            evidence_pos = Some((pos_in_read, pos_in_ref));
            break;
        }
    }
    let (ev_read_pos, ev_ref_pos) = evidence_pos.wrap_err("no A base found in read")?;

    let calls = FxHashMap::from_iter([
        (6103081, RastairCall::Cpg { base: G, methylated: Some(false) }),
        (6103093, RastairCall::Cpg { base: G, methylated: Some(false) }),
        (ev_ref_pos, RastairCall::Cpg { base: G, methylated: Some(true) }),
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
    rewrite_record(&calls, &mut legacy_record, BamMode::Legacy, ref_lookup)?;

    let xm_tag = string_tag(&legacy_record, b"XM")?;
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
    rewrite_record(&calls, &mut standard_record, BamMode::Standard, ref_lookup)?;

    let mm_tag = string_tag(&standard_record, b"MM")?;

    // This is a reversed read, so MM positions are in original-read orientation
    // (reverse complement of stored SEQ). The strand qualifier flips OB→OT.
    assert_snapshot!(mm_tag.chars().take(3).collect::<String>(), @"C+m");

    // === Cross-check ===
    let xm_methylated = decode_xm_to_positions(&xm_tag);
    assert_eq!(xm_methylated, vec![ev_read_pos], "XM: only evidence position should be Z");

    let fwd_seq = seq_for_mm_tag(&standard_record);
    let (mm_base, mm_methylated) = decode_mm_to_positions(&mm_tag, &fwd_seq)?;
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
    let mut record = first_overlapping(6103075, 6103100)?;

    // just to be sure we have the expected record
    assert_snapshot!(pos(&record), @"6103075");
    let strand = Strand::from_flags(record.flags()).wrap_err("unknown strand")?;
    assert_snapshot!(strand, @"OB");

    let ol_seq = stored_seq(&record);

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

        for (pos_in_read, pos_in_ref) in matched(&record) {
            let current_base = Base::from(seq[pos_in_read]);
            let Some(next_base) = seq.get(pos_in_read + 1).map(Base::from) else {
                continue;
            };

            match (current_base, next_base) {
                (C, G) => {
                    calls.insert(pos_in_ref, RastairCall::Cpg { methylated: Some(true), base: C });
                }
                (C, A) => {
                    calls.insert(
                        pos_in_ref + 1,
                        RastairCall::Cpg { methylated: Some(true), base: G },
                    );
                }
                _ => {}
            }
        }
        calls
    };

    let data = get_methylated_positions(&calls, &record, strand);
    let new_seq = &data.seq;

    // for human comparison
    assert_snapshot!(as_base_string(&ol_seq), @"AAAGGCGTGCACCACCACGCCTGGCTTGGTTTGGTTTTTGATTGGTTGGTTGGTCTTTTGAGACAGGGTTTCTCTGTGTA");
    assert_snapshot!(as_base_string(new_seq), @"AAAGGCGTGCGCCGCCGCGCCTGGCTTGGTTTGGTTTTTGATTGGTTGGTTGGTCTTTTGAGACGGGGTTTCTCTGTGTA");

    assert_compact_debug_snapshot!(
        data.methylated_positions.iter().map(|p| p.get()).collect::<Vec<_>>(),
        @"[10, 13, 16, 64]"
    );
    // verify these are all positions where A (methylation evidence) was
    // rewritten to G
    for &pos in &data.methylated_positions {
        let pos = pos.as_usize();
        let base = Base::from(ol_seq[pos]);
        if base != A {
            bail!("expected A at methylated position {}, found {}", pos, base.as_str());
        }
        let new_base = Base::from(new_seq[pos]);
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
    let is_reverse = record.flags().is_reverse();
    assert!(is_reverse, "test record should be reverse-strand");

    let seq_len = new_seq.len();
    let original_positions: SmallVec<u32, 10> = data
        .methylated_positions
        .iter()
        .map(|&p| u32::try_from(seq_len - 1 - p.as_usize()).expect("fits"))
        .collect();
    let original_seq = reverse_complement(new_seq);

    let methylated_positions =
        MethylatedPositions::new(Strand::OT, &original_seq, &original_positions);

    assert_compact_debug_snapshot!(methylated_positions, @"MethylatedPositions { base: C, strand: OT, positions: [5, 16, 0, 0] }");
    methylated_positions
        .apply_to_record(&mut record)
        .wrap_err("failed to apply modifications to record")?;

    let mod_string = string_tag(&record, b"MM")?;
    assert_snapshot!(mod_string, @"C+m,5,16,0,0;");

    Ok(())
}

#[test]
fn test_xr_xg_xm_tags() -> Result<()> {
    use Base::*;

    // Test with flag 83 (OB strand, first in pair reverse)
    let mut record = first_overlapping(6103075, 6103100)?;

    // Verify this is flag 83
    assert_eq!(flag(&record), 83, "Expected flag 83");
    assert_eq!(Strand::from_flags(record.flags()), Some(Strand::OB));

    let ol_seq = stored_seq(&record);

    // Create test calls with some methylated CpGs
    let calls = {
        let mut calls = FxHashMap::default();
        let seq = &ol_seq;

        for (pos_in_read, pos_in_ref) in matched(&record) {
            let current_base = Base::from(seq[pos_in_read]);
            let Some(next_base) = seq.get(pos_in_read + 1).map(Base::from) else {
                continue;
            };

            match (current_base, next_base) {
                (C, G) => {
                    calls.insert(pos_in_ref, RastairCall::Cpg { methylated: Some(true), base: C });
                }
                (C, A) => {
                    calls.insert(
                        pos_in_ref + 1,
                        RastairCall::Cpg { methylated: Some(true), base: G },
                    );
                }
                _ => {}
            }
        }
        calls
    };

    rewrite_record(&calls, &mut record, BamMode::Legacy, no_ref)?;

    // Check XR tag (flag 83: first in pair => CT)
    let xr_tag = string_tag(&record, b"XR")?;
    assert_snapshot!(xr_tag, @"CT");

    // Check XG tag (OB strand => GA)
    let xg_tag = string_tag(&record, b"XG")?;
    assert_snapshot!(xg_tag, @"GA");

    // Check XM tag exists and has correct length
    let xm_tag = string_tag(&record, b"XM")?;
    assert_eq!(xm_tag.len(), record.seq_len());

    Ok(())
}

/// A mate flagged neither first nor second has no strand to read methylation
/// from, so neither mode calls any of its bases, even at called CpGs.
#[test]
fn a_read_of_unknown_strand_carries_no_methylation() -> Result<()> {
    let mut original = first_overlapping(6103075, 6103100)?;
    original.set_flags(BamFlags::from(flag(&original) & !0xC0))?;
    ensure!(Strand::from_flags(original.flags()).is_none(), "the strand should be unknown");
    let calls: FxHashMap<u32, RastairCall> = matched(&original)
        .into_iter()
        .flat_map(|(_, pos)| {
            [
                (pos, RastairCall::Cpg { methylated: Some(true), base: Base::C }),
                (pos, RastairCall::Cpg { methylated: Some(true), base: Base::G }),
            ]
        })
        .collect();

    let mut legacy = original.clone();
    rewrite_record(&calls, &mut legacy, BamMode::Legacy, no_ref)?;
    assert_eq!(string_tag(&legacy, b"XR")?, "GA");
    assert_eq!(string_tag(&legacy, b"XG")?, "CT");
    assert_eq!(string_tag(&legacy, b"XM")?, ".".repeat(original.seq_len()));

    let mut standard = original.clone();
    rewrite_record(&calls, &mut standard, BamMode::Standard, no_ref)?;
    assert!(!has_tag(&standard, b"MM") && !has_tag(&standard, b"ML"));
    assert_eq!(stored_seq(&standard), stored_seq(&original));
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

    for (wanted, expected_xr, expected_xg) in test_cases {
        let mut record = find_record(|r| flag(r) == wanted)
            .wrap_err_with(|| format!("No record with flag {wanted} found in test BAM"))?;

        let calls = FxHashMap::default();
        rewrite_record(&calls, &mut record, BamMode::Legacy, no_ref)?;

        let xr_tag = string_tag(&record, b"XR")?;
        assert_eq!(xr_tag, expected_xr, "XR tag mismatch for flag {}", wanted);

        let xg_tag = string_tag(&record, b"XG")?;
        assert_eq!(xg_tag, expected_xg, "XG tag mismatch for flag {}", wanted);
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

    let base_positions: Vec<usize> =
        seq.iter().enumerate().filter(|&(_, b)| Base::from(*b) == base).map(|(i, _)| i).collect();

    let mut positions = Vec::new();
    let mut base_idx: usize = 0;

    for skip in skips_str.split(',') {
        let skip: usize = skip.parse().wrap_err_with(|| format!("bad skip: {skip}"))?;
        base_idx += skip;
        let seq_pos = base_positions
            .get(base_idx)
            .ok_or_else(|| color_eyre::eyre::eyre!("MM skip out of bounds: base_idx={base_idx}"))?;
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
fn seq_for_mm_tag(record: &RawBamRecord) -> Vec<u8> {
    let stored = stored_seq(record);
    if record.flags().is_reverse() { reverse_complement(&stored) } else { stored }
}

#[test]
fn roundtrip_standard_mm_tags() -> Result<()> {
    let mut tested = 0u32;

    for mut record in records_overlapping(6103075, 6103200)? {
        let ol_seq = stored_seq(&record);
        let calls = build_cpg_calls(&record, &ol_seq);
        rewrite_record(&calls, &mut record, BamMode::Standard, no_ref)?;

        // Standard mode should NOT produce XR/XG/XM tags
        assert!(!has_tag(&record, b"XR"), "Standard mode should not produce XR tag");

        // MM/ML are omitted when there are no methylated positions (absent = "no data").
        if has_tag(&record, b"MM") {
            let mm_tag = string_tag(&record, b"MM")?;

            let fwd_seq = seq_for_mm_tag(&record);
            let (_, mm_positions) =
                decode_mm_to_positions(&mm_tag, &fwd_seq).wrap_err_with(|| {
                    format!("decode MM {:?} pos={} flag={}", mm_tag, pos(&record), flag(&record))
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

    let mut rewritten = records_overlapping(6103075, 6103100)?;
    for record in &mut rewritten {
        let calls = build_cpg_calls(record, &stored_seq(record));
        rewrite_record(&calls, record, BamMode::Standard, no_ref)?;
    }
    ensure!(!rewritten.is_empty(), "no records written");
    write_bam(&temp_bam, &rewritten)?;

    let read_back = read_bam(&temp_bam)?;
    for record in &read_back {
        let fwd_seq = seq_for_mm_tag(record);

        // MM/ML are omitted when there are no methylated positions (absent = "no data").
        if has_tag(record, b"MM") {
            let mm_tag = string_tag(record, b"MM")?;
            let (_, mm_positions) = decode_mm_to_positions(&mm_tag, &fwd_seq)?;

            let Some(AuxValue::ArrayU8(ml_data)) = find_tag(record.aux().as_bytes(), *b"ML") else {
                bail!("ML missing or not an array");
            };
            assert_eq!(
                ml_data.len(),
                mm_positions.len(),
                "ML length should match number of methylated positions"
            );
        }

        assert!(!has_tag(record, b"XR"), "Standard mode should not produce XR");
    }

    assert_eq!(rewritten.len(), read_back.len(), "record count mismatch after roundtrip");

    Ok(())
}

#[test]
fn roundtrip_write_read_bam_legacy() -> Result<()> {
    let temp_dir = tempfile::TempDir::new()?;
    let temp_bam = temp_dir.path().join("roundtrip_legacy.bam");

    let mut rewritten = records_overlapping(6103075, 6103100)?;
    for record in &mut rewritten {
        let calls = build_cpg_calls(record, &stored_seq(record));
        rewrite_record(&calls, record, BamMode::Legacy, no_ref)?;
    }
    ensure!(!rewritten.is_empty(), "no records written");
    write_bam(&temp_bam, &rewritten)?;

    let read_back = read_bam(&temp_bam)?;
    for record in &read_back {
        string_tag(record, b"XR")?;
        string_tag(record, b"XG")?;
        let xm_tag = string_tag(record, b"XM")?;
        assert_eq!(xm_tag.len(), record.seq_len(), "XM length mismatch");

        assert!(!has_tag(record, b"MM"), "Legacy mode should not produce MM");
    }

    assert_eq!(rewritten.len(), read_back.len(), "record count mismatch after roundtrip");

    Ok(())
}

#[test]
fn roundtrip_no_methylation_standard() -> Result<()> {
    let temp_dir = tempfile::TempDir::new()?;
    let temp_bam = temp_dir.path().join("roundtrip_empty.bam");

    let mut record = first_overlapping(6103075, 6103100)?;
    let calls = FxHashMap::default();
    rewrite_record(&calls, &mut record, BamMode::Standard, no_ref)?;
    write_bam(&temp_bam, &[record])?;

    let record = read_bam(&temp_bam)?.into_iter().next().wrap_err("no records")?;

    // With no methylated positions, MM/ML tags are not written at all.
    // An absent MM means "no modification data" per SAM spec.
    // Writing empty C+m; or empty ML:B:C crashes tools like modbedtools.
    assert!(
        !has_tag(&record, b"MM"),
        "MM tag should not be present when there are no methylated positions"
    );
    assert!(
        !has_tag(&record, b"ML"),
        "ML tag should not be present when there are no methylated positions"
    );

    Ok(())
}

#[test]
fn roundtrip_no_methylation_legacy() -> Result<()> {
    let temp_dir = tempfile::TempDir::new()?;
    let temp_bam = temp_dir.path().join("roundtrip_empty_legacy.bam");

    let mut record = first_overlapping(6103075, 6103100)?;
    let calls = FxHashMap::default();
    rewrite_record(&calls, &mut record, BamMode::Legacy, no_ref)?;
    write_bam(&temp_bam, &[record])?;

    let record = read_bam(&temp_bam)?.into_iter().next().wrap_err("no records")?;

    let xm_tag = string_tag(&record, b"XM")?;
    let xm_positions = decode_xm_to_positions(&xm_tag);
    assert!(xm_positions.is_empty(), "expected no XM methylation");

    ensure!(!xm_tag.contains('Z'), "XM should have no uppercase Z without methylation calls");

    Ok(())
}

/// Build CpG methylation calls from a record's sequence by looking for
/// CA (OB strand evidence) and CG dinucleotides.
fn build_cpg_calls(record: &RawBamRecord, seq: &[u8]) -> FxHashMap<u32, RastairCall> {
    use Base::*;
    let mut calls = FxHashMap::default();
    for (pos_in_read, pos_in_ref) in matched(record) {
        let current_base = Base::from(seq[pos_in_read]);
        let Some(next_base) = seq.get(pos_in_read + 1).map(Base::from) else {
            continue;
        };
        match (current_base, next_base) {
            (C, G) => {
                calls.insert(pos_in_ref, RastairCall::Cpg { methylated: Some(true), base: C });
            }
            (C, A) => {
                calls.insert(pos_in_ref + 1, RastairCall::Cpg { methylated: Some(true), base: G });
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

    let mut record = first_overlapping(6103075, 6103100)?;

    let strand = Strand::from_flags(record.flags()).wrap_err("unknown strand")?;
    let seq = stored_seq(&record);

    // Count how many G bases (OB target) are in the read
    let total_target_bases =
        seq.iter().filter(|&&b| Base::from(b) == if strand == Strand::OB { G } else { C }).count();

    // Create calls for only ONE CpG position — pick the first CG or CA
    // dinucleotide we find
    let mut calls = FxHashMap::default();
    let mut cpg_count = 0u32;
    for (pos_in_read, pos_in_ref) in matched(&record) {
        let current_base = Base::from(seq[pos_in_read]);
        let Some(next_base) = seq.get(pos_in_read + 1).map(Base::from) else {
            continue;
        };

        // For OB: look for CA (methylation evidence for G on next pos)
        if strand == Strand::OB && current_base == C && next_base == A && cpg_count == 0 {
            calls.insert(pos_in_ref + 1, RastairCall::Cpg { methylated: Some(true), base: G });
            cpg_count += 1;
            break;
        }
        // For OT: look for TG (methylation evidence for C)
        if strand == Strand::OT && current_base == T && next_base == G && cpg_count == 0 {
            calls.insert(pos_in_ref, RastairCall::Cpg { methylated: Some(true), base: C });
            cpg_count += 1;
            break;
        }
    }
    ensure!(cpg_count == 1, "need at least one CpG call for this test");

    rewrite_record(&calls, &mut record, BamMode::Legacy, no_ref)?;

    let xm_tag = string_tag(&record, b"XM")?;

    let z_count = xm_tag.chars().filter(|c| *c == 'Z' || *c == 'z').count();

    assert!(
        z_count < total_target_bases,
        "XM has {z_count} CpG annotations but only {cpg_count} CpG call(s) exist. \
         There are {total_target_bases} target bases in the read — non-CpG positions \
         should be '.' not 'z'. XM: {xm_tag}"
    );

    Ok(())
}

/// De-novo CpG calls (`DeNovoCpg` variant) must also appear in legacy XM tags.
#[test]
fn xm_includes_denovo_cpg_methylation() -> Result<()> {
    use Base::*;

    let mut record = first_overlapping(6103075, 6103100)?;

    let strand = Strand::from_flags(record.flags()).wrap_err("unknown strand")?;
    let seq = stored_seq(&record);

    // Create a single DeNovoCpg call at the first methylation evidence base
    let mut calls = FxHashMap::default();
    let mut found = false;
    for (pos_in_read, pos_in_ref) in matched(&record) {
        let observed = Base::from(seq[pos_in_read]);
        let is_evidence = match strand {
            Strand::OT => observed == T,
            Strand::OB => observed == A,
        };
        if is_evidence {
            let base = match strand {
                Strand::OT => C,
                Strand::OB => G,
            };
            calls.insert(pos_in_ref, RastairCall::DeNovoCpg { methylated: Some(true), base });
            found = true;
            break;
        }
    }
    ensure!(found, "need a methylation-evidence base for this test");

    rewrite_record(&calls, &mut record, BamMode::Legacy, no_ref)?;

    let xm_tag = string_tag(&record, b"XM")?;

    let has_methylated = xm_tag.chars().any(|c| c == 'Z' || c == 'X' || c == 'H');
    assert!(
        has_methylated,
        "XM tag should contain a methylation mark for the DeNovoCpg call, \
         but got: {xm_tag}"
    );

    Ok(())
}

/// De-novo CpG calls must also be found by `find_methylated_positions` (standard mode).
#[test]
fn standard_mode_includes_denovo_cpg() -> Result<()> {
    use Base::*;

    let record = first_overlapping(6103075, 6103100)?;

    let strand = Strand::from_flags(record.flags()).wrap_err("unknown strand")?;
    let seq = stored_seq(&record);

    let mut calls = FxHashMap::default();
    let mut expected_pos = None;
    for (pos_in_read, pos_in_ref) in matched(&record) {
        let observed = Base::from(seq[pos_in_read]);
        let is_evidence = match strand {
            Strand::OT => observed == T,
            Strand::OB => observed == A,
        };
        if is_evidence {
            let base = match strand {
                Strand::OT => C,
                Strand::OB => G,
            };
            calls.insert(pos_in_ref, RastairCall::DeNovoCpg { methylated: Some(true), base });
            expected_pos = Some(u32::try_from(pos_in_read)?);
            break;
        }
    }
    ensure!(expected_pos.is_some(), "need a methylation-evidence base");

    let data = get_methylated_positions(&calls, &record, strand);
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

    let mut record = first_overlapping(6103075, 6103100)?;

    let strand = Strand::from_flags(record.flags()).wrap_err("unknown strand")?;
    let seq = stored_seq(&record);

    // Create a DeNovoCpg call and a synthetic ref lookup that returns
    // non-CpG context (CHH: no G after C for OT, no C before G for OB)
    let mut calls = FxHashMap::default();
    let mut denovo_ref_pos = None;
    for (pos_in_read, pos_in_ref) in matched(&record) {
        let observed = Base::from(seq[pos_in_read]);
        let is_evidence = match strand {
            Strand::OT => observed == T,
            Strand::OB => observed == A,
        };
        if is_evidence {
            let base = match strand {
                Strand::OT => C,
                Strand::OB => G,
            };
            calls.insert(pos_in_ref, RastairCall::DeNovoCpg { methylated: Some(true), base });
            denovo_ref_pos = Some(pos_in_ref);
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
            }
        } else {
            Some(A) // non-CpG neighbor → CHH
        }
    };

    rewrite_record(&calls, &mut record, BamMode::Legacy, ref_lookup)?;

    let xm_tag = string_tag(&record, b"XM")?;

    // De-novo with CHH context should use H (methylated CHH), not Z
    assert!(xm_tag.contains('H'), "De-novo CpG in CHH context should produce 'H', got: {xm_tag}");
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
    let record = find_flag_at(99, 6103078)?;

    let strand = Strand::from_flags(record.flags()).wrap_err("unknown strand")?;
    assert_eq!(strand, Strand::OT);
    assert!(!record.flags().is_reverse());
    let seq = stored_seq(&record);
    assert_snapshot!(as_base_string(&seq), @"GGTGTGCACCACCATGCCTGGCTTGGTTTGGTTTTTGATTGGTTGGTTGGTCTTTTGAGACAGGGTTTCTCTGTGTAGCT");

    // CpG calls (C side for OT):
    // Ref pos 6103080 (read offset 2): methylated
    //   Read[2] = T → methylation evidence on OT
    // Ref pos 6103092 (read offset 14): methylated
    //   Read[14] = T → methylation evidence on OT
    let calls = FxHashMap::from_iter([
        (6103080, RastairCall::Cpg { base: C, methylated: Some(true) }),
        (6103092, RastairCall::Cpg { base: C, methylated: Some(true) }),
    ]);

    let ref_bases = FxHashMap::from_iter([(6103080, C), (6103081, G), (6103092, C), (6103093, G)]);
    let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };

    // === Legacy mode ===
    let mut legacy_record = record.clone();
    rewrite_record(&calls, &mut legacy_record, BamMode::Legacy, ref_lookup)?;

    let xm_tag = string_tag(&legacy_record, b"XM")?;
    let xr_tag = string_tag(&legacy_record, b"XR")?;
    let xg_tag = string_tag(&legacy_record, b"XG")?;

    // Flag 99: first in pair → XR=CT; OT → XG=CT
    assert_snapshot!(xr_tag, @"CT");
    assert_snapshot!(xg_tag, @"CT");

    let xm_annotations: Vec<(usize, char)> =
        xm_tag.chars().enumerate().filter(|(_, c)| *c != '.').collect();
    assert_compact_debug_snapshot!(xm_annotations, @"[(2, 'Z'), (14, 'Z')]");

    // SEQ unchanged in legacy mode
    assert_eq!(stored_seq(&legacy_record), seq);

    // === Standard mode ===
    let mut standard_record = record.clone();
    rewrite_record(&calls, &mut standard_record, BamMode::Standard, ref_lookup)?;

    let mm_tag = string_tag(&standard_record, b"MM")?;
    // OT, forward → C+m
    assert_snapshot!(mm_tag.chars().take(3).collect::<String>(), @"C+m");
    assert_snapshot!(mm_tag, @"C+m,0,5;");

    // Rewritten SEQ: T→C at methylated positions 2 and 14
    let new_seq = stored_seq(&standard_record);
    assert_eq!(Base::from(seq[2]), T);
    assert_eq!(Base::from(new_seq[2]), C);
    assert_eq!(Base::from(seq[14]), T);
    assert_eq!(Base::from(new_seq[14]), C);

    // === Cross-check ===
    let xm_methylated = decode_xm_to_positions(&xm_tag);
    let fwd_seq = seq_for_mm_tag(&standard_record);
    let (mm_base, mm_methylated) = decode_mm_to_positions(&mm_tag, &fwd_seq)?;
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
    let record = find_flag_at(163, 6103076)?;

    let strand = Strand::from_flags(record.flags()).wrap_err("unknown strand")?;
    assert_eq!(strand, Strand::OB);
    let seq = stored_seq(&record);

    // Read[5] = A → methylation evidence at ref pos 6103081 on OB
    assert_eq!(Base::from(seq[5]), A);

    // Make it a DeNovoCpg with CHG context:
    // OB context looks at ref[pos-1] and ref[pos-2] (complement perspective)
    // ref[pos-1]=T, ref[pos-2]=C → CHG (second base back is C)
    let calls = FxHashMap::from_iter([(
        6103081_u32,
        RastairCall::DeNovoCpg { base: G, methylated: Some(true) },
    )]);
    let ref_bases = FxHashMap::from_iter([
        (6103079_u32, C), // pos-2 → C → CHG
        (6103080, T),     // pos-1 → T (complement=A, not G, so not CpG)
        (6103081, G),
    ]);
    let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };

    // === Legacy mode ===
    let mut legacy_record = record.clone();
    rewrite_record(&calls, &mut legacy_record, BamMode::Legacy, ref_lookup)?;

    let xm_tag = string_tag(&legacy_record, b"XM")?;
    let xm_annotations: Vec<(usize, char)> =
        xm_tag.chars().enumerate().filter(|(_, c)| *c != '.').collect();
    // CHG methylated → 'X'
    assert_compact_debug_snapshot!(xm_annotations, @"[(5, 'X')]");

    // === Standard mode ===
    let mut standard_record = record.clone();
    rewrite_record(&calls, &mut standard_record, BamMode::Standard, ref_lookup)?;

    let mm_tag = string_tag(&standard_record, b"MM")?;
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
    let (_, mm_methylated) = decode_mm_to_positions(&mm_tag, &fwd_seq)?;
    assert_eq!(xm_methylated, vec![5]);
    assert_eq!(mm_methylated, vec![5]);

    Ok(())
}

/// De-novo CpGs must produce context-dependent XM letters (z/Z for CpG,
/// x/X for CHG, h/H for CHH) through the full `rewrite_record` path, for
/// both OT and OB strands, and both methylated and unmethylated reads.
#[test]
fn denovo_all_contexts_both_strands() -> Result<()> {
    use Base::*;

    // === OB strand (flag-163 read at 6103076) ===
    let ob_record = find_flag_at(163, 6103076)?;
    let ob_seq = stored_seq(&ob_record);
    // Read[5]=A (evidence), Read[17]=G (target) on OB
    assert_eq!(Base::from(ob_seq[5]), A);
    assert_eq!(Base::from(ob_seq[17]), G);

    // OB de-novo CpG context: ref[pos-1]==C → CpG
    {
        let calls = FxHashMap::from_iter([
            (6103081_u32, RastairCall::DeNovoCpg { base: G, methylated: Some(true) }),
            (6103093, RastairCall::DeNovoCpg { base: G, methylated: Some(false) }),
        ]);
        let ref_bases = FxHashMap::from_iter([
            (6103080_u32, C), // pos-1 for 6103081 → CpG
            (6103081, G),
            (6103092, C), // pos-1 for 6103093 → CpG
            (6103093, G),
        ]);
        let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };
        let mut record = ob_record.clone();
        rewrite_record(&calls, &mut record, BamMode::Legacy, ref_lookup)?;
        let xm = string_tag(&record, b"XM")?;
        let annotations: Vec<(usize, char)> =
            xm.chars().enumerate().filter(|(_, c)| *c != '.').collect();
        assert_compact_debug_snapshot!(annotations, @"[(5, 'Z'), (17, 'z')]");
    }

    // OB de-novo CHG context: ref[pos-1]!=C, ref[pos-2]==C → CHG
    {
        let calls = FxHashMap::from_iter([
            (6103081_u32, RastairCall::DeNovoCpg { base: G, methylated: Some(true) }),
            (6103093, RastairCall::DeNovoCpg { base: G, methylated: Some(false) }),
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
        rewrite_record(&calls, &mut record, BamMode::Legacy, ref_lookup)?;
        let xm = string_tag(&record, b"XM")?;
        let annotations: Vec<(usize, char)> =
            xm.chars().enumerate().filter(|(_, c)| *c != '.').collect();
        assert_compact_debug_snapshot!(annotations, @"[(5, 'X'), (17, 'x')]");
    }

    // OB de-novo CHH context: ref[pos-1]!=C, ref[pos-2]!=C → CHH
    {
        let calls = FxHashMap::from_iter([
            (6103081_u32, RastairCall::DeNovoCpg { base: G, methylated: Some(true) }),
            (6103093, RastairCall::DeNovoCpg { base: G, methylated: Some(false) }),
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
        rewrite_record(&calls, &mut record, BamMode::Legacy, ref_lookup)?;
        let xm = string_tag(&record, b"XM")?;
        let annotations: Vec<(usize, char)> =
            xm.chars().enumerate().filter(|(_, c)| *c != '.').collect();
        assert_compact_debug_snapshot!(annotations, @"[(5, 'H'), (17, 'h')]");
    }

    // === OT strand (flag-99 read at 6103078) ===
    let ot_record = find_flag_at(99, 6103078)?;
    let ot_seq = stored_seq(&ot_record);
    // Read[2]=T (evidence), Read[14]=T (evidence) on OT
    assert_eq!(Base::from(ot_seq[2]), T);
    assert_eq!(Base::from(ot_seq[14]), T);

    // OT de-novo CpG context: ref[pos+1]==G → CpG
    {
        let calls = FxHashMap::from_iter([
            (6103080_u32, RastairCall::DeNovoCpg { base: C, methylated: Some(true) }),
            (6103092, RastairCall::DeNovoCpg { base: C, methylated: Some(true) }),
        ]);
        let ref_bases = FxHashMap::from_iter([
            (6103080_u32, C),
            (6103081, G), // pos+1 → CpG
            (6103092, C),
            (6103093, G), // pos+1 → CpG
        ]);
        let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };
        let mut record = ot_record.clone();
        rewrite_record(&calls, &mut record, BamMode::Legacy, ref_lookup)?;
        let xm = string_tag(&record, b"XM")?;
        let annotations: Vec<(usize, char)> =
            xm.chars().enumerate().filter(|(_, c)| *c != '.').collect();
        assert_compact_debug_snapshot!(annotations, @"[(2, 'Z'), (14, 'Z')]");
    }

    // OT de-novo CHG context: ref[pos+1]!=G, ref[pos+2]==G → CHG
    {
        let calls = FxHashMap::from_iter([
            (6103080_u32, RastairCall::DeNovoCpg { base: C, methylated: Some(true) }),
            (6103092, RastairCall::DeNovoCpg { base: C, methylated: Some(true) }),
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
        rewrite_record(&calls, &mut record, BamMode::Legacy, ref_lookup)?;
        let xm = string_tag(&record, b"XM")?;
        let annotations: Vec<(usize, char)> =
            xm.chars().enumerate().filter(|(_, c)| *c != '.').collect();
        assert_compact_debug_snapshot!(annotations, @"[(2, 'X'), (14, 'X')]");
    }

    // OT de-novo CHH context: ref[pos+1]!=G, ref[pos+2]!=G → CHH
    {
        let calls = FxHashMap::from_iter([
            (6103080_u32, RastairCall::DeNovoCpg { base: C, methylated: Some(true) }),
            (6103092, RastairCall::DeNovoCpg { base: C, methylated: Some(true) }),
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
        rewrite_record(&calls, &mut record, BamMode::Legacy, ref_lookup)?;
        let xm = string_tag(&record, b"XM")?;
        let annotations: Vec<(usize, char)> =
            xm.chars().enumerate().filter(|(_, c)| *c != '.').collect();
        assert_compact_debug_snapshot!(annotations, @"[(2, 'H'), (14, 'H')]");
    }

    Ok(())
}

/// Reads with indels (insertions/deletions) must have correct position
/// mapping through `aligned_pairs_full`. Verify a read with a 1bp deletion
/// (61M1D19M) produces correct XM and MM output.
#[test]
fn golden_indel_read_legacy_and_standard_agree() -> Result<()> {
    use Base::*;

    // Flag 99 read at 6106220 (BAM 0-based) with CIGAR 61M1D19M
    let record = find_flag_at(99, 6106220)?;

    let strand = Strand::from_flags(record.flags()).wrap_err("unknown strand")?;
    assert_eq!(strand, Strand::OT);
    let seq = stored_seq(&record);

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
        (6106229_u32, RastairCall::Cpg { base: C, methylated: Some(true) }),
        (6106300, RastairCall::Cpg { base: C, methylated: Some(false) }),
    ]);
    let ref_bases =
        FxHashMap::from_iter([(6106229_u32, C), (6106230, G), (6106300, C), (6106301, G)]);
    let ref_lookup = |pos: u32| -> Option<Base> { ref_bases.get(&pos).copied() };

    // === Legacy mode ===
    let mut legacy_record = record.clone();
    rewrite_record(&calls, &mut legacy_record, BamMode::Legacy, ref_lookup)?;

    let xm_tag = string_tag(&legacy_record, b"XM")?;
    assert_eq!(xm_tag.len(), seq.len());
    let xm_annotations: Vec<(usize, char)> =
        xm_tag.chars().enumerate().filter(|(_, c)| *c != '.').collect();
    // Position 9 = Z (methylated), position 79 = z (unmethylated, AFTER deletion)
    assert_compact_debug_snapshot!(xm_annotations, @"[(9, 'Z'), (79, 'z')]");

    // === Standard mode ===
    let mut standard_record = record.clone();
    rewrite_record(&calls, &mut standard_record, BamMode::Standard, ref_lookup)?;

    let mm_tag = string_tag(&standard_record, b"MM")?;

    // Rewritten SEQ: T→C at position 9
    let new_seq = stored_seq(&standard_record);
    assert_eq!(Base::from(new_seq[9]), C, "methylated T should be rewritten to C");
    assert_eq!(Base::from(new_seq[79]), C, "unmethylated C should stay C");

    // === Cross-check ===
    let xm_methylated = decode_xm_to_positions(&xm_tag);
    let fwd_seq = seq_for_mm_tag(&standard_record);
    let (mm_base, mm_methylated) = decode_mm_to_positions(&mm_tag, &fwd_seq)?;
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

    for (wanted, expected_strand, expected_xr, expected_xg) in test_cases {
        // Make sure this read covers our CpG position
        let record = find_record(|r| {
            let start = pos(r);
            flag(r) == wanted && start <= 6103081 && start + r.seq_len() as u32 > 6103081
        })
        .wrap_err_with(|| format!("No record with flag {wanted} covering CpG at 6103080/81"))?;

        let strand = Strand::from_flags(record.flags()).wrap_err("unknown strand")?;
        assert_eq!(strand, expected_strand, "strand mismatch for flag {wanted}");

        // Build calls appropriate for this strand
        let calls = if strand == Strand::OT {
            FxHashMap::from_iter([(
                6103080_u32,
                RastairCall::Cpg { base: C, methylated: Some(true) },
            )])
        } else {
            FxHashMap::from_iter([(
                6103081_u32,
                RastairCall::Cpg { base: G, methylated: Some(true) },
            )])
        };

        // === Legacy ===
        let mut legacy_record = record.clone();
        rewrite_record(&calls, &mut legacy_record, BamMode::Legacy, ref_lookup)?;

        let xr_tag = string_tag(&legacy_record, b"XR")?;
        let xg_tag = string_tag(&legacy_record, b"XG")?;
        let xm_tag = string_tag(&legacy_record, b"XM")?;

        assert_eq!(xr_tag, expected_xr, "XR mismatch for flag {wanted}");
        assert_eq!(xg_tag, expected_xg, "XG mismatch for flag {wanted}");
        assert_eq!(xm_tag.len(), record.seq_len(), "XM length mismatch for flag {wanted}");

        // Every flag should have at least one annotation (Z or z)
        let has_annotation = xm_tag.chars().any(|c| c != '.');
        assert!(has_annotation, "XM has no annotations for flag {wanted}: {xm_tag}");

        // === Standard ===
        let mut standard_record = record.clone();
        rewrite_record(&calls, &mut standard_record, BamMode::Standard, ref_lookup)?;

        // Cross-check: both modes find the same number of methylated positions.
        // MM/ML are absent when no methylation evidence is found (absent = 0 positions).
        let xm_methylated = decode_xm_to_positions(&xm_tag);
        let mm_methylated_count = if has_tag(&standard_record, b"MM") {
            let mm_tag = string_tag(&standard_record, b"MM")?;
            let fwd_seq = seq_for_mm_tag(&standard_record);
            let (_, positions) = decode_mm_to_positions(&mm_tag, &fwd_seq)?;
            positions.len()
        } else {
            0
        };
        assert_eq!(
            xm_methylated.len(),
            mm_methylated_count,
            "Methylated position count mismatch for flag {wanted}: \
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

/// A rerun on rastair's own output gets the next free `@PG` ID, chained to
/// the last one, rather than one whose `PP` names itself.
#[test]
fn a_rerun_gets_its_own_pg_id() -> Result<()> {
    let header = |pgs: &str| BamHeader::from_sam_text(&format!("@SQ\tSN:c1\tLN:10\n{pgs}"));
    let pg_ids = |header: &BamHeader| -> Vec<String> {
        header
            .header_text()
            .lines()
            .filter(|line| line.starts_with("@PG"))
            .map(|line| {
                line.split('\t')
                    .filter(|f| !f.starts_with("CL:") && !f.starts_with("VN:"))
                    .collect::<Vec<_>>()
                    .join(" ")
            })
            .collect()
    };

    let mut fresh = header("")?;
    fresh.add_pg(rastair_pg())?;
    assert_snapshot!(pg_ids(&fresh).join("\n"), @"@PG ID:rastair.rewrite_bam PN:rastair DS:Rewrote BAM with methylation information");

    let mut rerun = header(
        "@PG\tID:rastair.rewrite_bam\tPN:rastair\n@PG\tID:rastair.rewrite_bam.1\tPN:rastair\tPP:rastair.rewrite_bam\n",
    )?;
    rerun.add_pg(rastair_pg())?;
    let last = pg_ids(&rerun).pop().wrap_err("no @PG")?;
    assert_snapshot!(last, @"@PG ID:rastair.rewrite_bam.2 PN:rastair DS:Rewrote BAM with methylation information PP:rastair.rewrite_bam.1");
    Ok(())
}

mod regions {
    use super::*;
    use crate::utils::CliRegionInput;
    use hegel::{TestCase, generators as gs};

    #[hegel::composite]
    fn spans(tc: &TestCase, max_len: u64) -> (u32, u64, u64) {
        let tid = tc.draw(gs::integers::<u32>().max_value(1));
        let start = tc.draw(gs::integers::<u64>().max_value(120));
        let len = tc.draw(gs::integers::<u64>().min_value(1).max_value(max_len));
        (tid, start, start + len)
    }

    /// `--region`s in any order, overlapping, adjacent or on two contigs:
    /// querying the merged regions in turn, skipping what the region before
    /// returned, yields each record overlapping one of them once, in
    /// (contig, start) order.
    #[hegel::test]
    fn regions_query_each_overlapping_record_once(tc: TestCase) {
        let n_regions = tc.draw(gs::integers::<usize>().max_value(6));
        let regions: Vec<(u32, u64, u64)> = (0..n_regions).map(|_| tc.draw(spans(40))).collect();
        let n_records = tc.draw(gs::integers::<usize>().max_value(40));
        let mut records: Vec<(u32, u64, u64)> =
            (0..n_records).map(|_| tc.draw(spans(30))).collect();
        records.sort_unstable();

        let header = BamHeader::synthetic(2, 200);
        let tid = |n: u32| header.tid(&format!("c{n}")).wrap_err("tid");
        let check = || -> Result<()> {
            let text: Vec<String> = regions
                .iter()
                .map(|(tid, start, end)| format!("c{tid}:{}-{end}", start + 1))
                .collect();
            let input: Option<CliRegionInput> =
                (!text.is_empty()).then(|| text.join(" ").parse()).transpose()?;
            let merged = Regions::new(input.as_ref(), &header)?;
            let span = |&(t, start, end): &(u32, u64, u64)| -> Result<Span> {
                Ok(Span { tid: tid(t)?, start, end })
            };
            let overlapping: Vec<Span> = records
                .iter()
                .filter(|(t, start, end)| {
                    regions.iter().any(|(rt, rs, re)| rt == t && start < re && rs < end)
                })
                .map(span)
                .collect::<Result<_>>()?;

            let mut queried = Vec::new();
            let mut previous = None;
            for &region in &merged.0 {
                for record in &records {
                    let record = span(record)?;
                    if record.overlaps(region) && !Regions::queried_before(previous, Some(record)) {
                        queried.push(record);
                    }
                }
                previous = Some(region);
            }
            ensure!(queried == overlapping, "queried {queried:?}, overlapping {overlapping:?}");
            Ok(())
        };
        let result = check();
        assert!(result.is_ok(), "{regions:?} {records:?}: {result:?}");
    }
}

/// A decoded CRAM container is several BAM batches, so a CRAM source runs one
/// batch per worker ahead of the writer, a BAM source the default.
#[test]
fn a_cram_runs_one_batch_per_worker_ahead() -> Result<()> {
    let window = |input: &str| -> Result<Option<usize>> {
        let params = ReaderParams {
            bam_file: ClioPath::new(input)?,
            fasta_file: ClioPath::new("tests/data/test.fasta.gz")?,
            regions: None,
        };
        Ok(Input::open(&params)?.source.window(7))
    };
    assert_eq!(window("tests/data/test.cram")?, Some(7));
    assert_eq!(window("tests/data/test.bam")?, None);
    Ok(())
}
