#![cfg(feature = "external-tool-tests")]

mod utils;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use utils::*;

struct TestBams {
    _temp_dir: TempDir,
    legacy: PathBuf,
    standard: PathBuf,
    calls_bed: PathBuf,
}

fn tool_is_available(name: &str) -> bool {
    std::process::Command::new(name)
        .arg("--version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

macro_rules! require_tool {
    ($name:expr) => {
        if !tool_is_available($name) {
            eprintln!("Skipping test: {} not found in PATH", $name);
            return Ok(());
        }
    };
}

fn setup_test_bams() -> Result<TestBams> {
    let temp_dir = TempDir::new()?;
    let calls_bed = temp_dir.path().join("calls.bed.gz");
    let legacy_bam = temp_dir.path().join("legacy.bam");
    let standard_bam = temp_dir.path().join("standard.bam");

    let region = "chr19:6103075-6103300";

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--no-ml",
            &format!("--region={region}"),
            "--cpgs-only",
            "-o",
        ])
        .arg(&calls_bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to generate calls")?;

    Command::new("tabix")
        .args(["-p", "bed"])
        .arg(&calls_bed)
        .output()
        .wrap_err("Failed to index calls file")?;

    rastair()
        .args([
            "bam",
            "legacy",
            "--fasta-file=tests/data/test.fasta.gz",
            &format!("--region={region}"),
            "tests/data/test.bam",
        ])
        .arg(&calls_bed)
        .args(["-o"])
        .arg(&legacy_bam)
        .silent()
        .succeeds()
        .wrap_err("Failed to rewrite legacy BAM")?;

    rastair()
        .args([
            "bam",
            "standard",
            "--fasta-file=tests/data/test.fasta.gz",
            &format!("--region={region}"),
            "tests/data/test.bam",
        ])
        .arg(&calls_bed)
        .args(["-o"])
        .arg(&standard_bam)
        .silent()
        .succeeds()
        .wrap_err("Failed to rewrite standard BAM")?;

    Ok(TestBams { _temp_dir: temp_dir, legacy: legacy_bam, standard: standard_bam, calls_bed })
}

fn samtools_count_with_tag(bam: &Path, tag: &str) -> Result<u32> {
    let output = Command::new("samtools")
        .args(["view", "-c", "-d", tag])
        .arg(bam)
        .output()
        .wrap_err("Failed to run samtools")?;
    ensure!(output.status.success(), "samtools view failed");
    let count: u32 = String::from_utf8_lossy(&output.stdout)
        .trim()
        .parse()
        .wrap_err("Failed to parse samtools count")?;
    Ok(count)
}

fn samtools_count_total(bam: &Path) -> Result<u32> {
    let output = Command::new("samtools")
        .args(["view", "-c"])
        .arg(bam)
        .output()
        .wrap_err("Failed to run samtools")?;
    ensure!(output.status.success(), "samtools view failed");
    let count: u32 = String::from_utf8_lossy(&output.stdout)
        .trim()
        .parse()
        .wrap_err("Failed to parse samtools count")?;
    Ok(count)
}

struct BedCall {
    pos: u32,
    strand: char,
    mod_count: u32,
    unmod_count: u32,
}

fn parse_calls_bed(path: &Path) -> Result<Vec<BedCall>> {
    let output = Command::new("gunzip")
        .args(["-c"])
        .arg(path)
        .output()
        .wrap_err("Failed to decompress calls BED")?;
    ensure!(output.status.success(), "gunzip failed");
    let text = String::from_utf8_lossy(&output.stdout);
    let mut calls = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        ensure!(fields.len() >= 8, "BED line has too few fields: {line}");
        let pos: u32 = fields[1].parse().wrap_err("parse start")?;
        let strand = fields[5].chars().next().unwrap_or('?');
        let unmod_count: u32 = fields[6].parse().wrap_err("parse unmod")?;
        let mod_count: u32 = fields[7].parse().wrap_err("parse mod")?;
        calls.push(BedCall { pos, strand, mod_count, unmod_count });
    }
    Ok(calls)
}

/// Collect per-position (ref_pos, strand) → (methylated, unmethylated) counts
/// from XM tags in a legacy BAM.
fn count_xm_per_position(bam_path: &Path) -> Result<HashMap<(u32, char), (u32, u32)>> {
    use rust_htslib::bam::{self, Read as BamRead};

    let mut counts: HashMap<(u32, char), (u32, u32)> = HashMap::new();
    let mut reader = bam::Reader::from_path(bam_path)?;
    let mut record = bam::Record::new();

    while let Some(result) = reader.read(&mut record) {
        result?;

        let xm = match record.aux(b"XM") {
            Ok(rust_htslib::bam::record::Aux::String(s)) => s.to_string(),
            _ => continue,
        };
        let is_ob = match record.aux(b"XG") {
            Ok(rust_htslib::bam::record::Aux::String(s)) => s == "GA",
            _ => continue,
        };
        let strand_char = if is_ob { '-' } else { '+' };

        let xm_bytes = xm.as_bytes();
        let cigar = record.cigar();
        let mut read_pos: usize = 0;
        let mut ref_pos = record.pos() as u32;

        for op in cigar.iter() {
            use rust_htslib::bam::record::Cigar::*;
            let len = op.len() as u32;
            match op {
                Match(_) | Equal(_) | Diff(_) => {
                    for _ in 0..len {
                        if let Some(&ch) = xm_bytes.get(read_pos) {
                            match ch {
                                b'Z' | b'X' | b'H' => {
                                    counts.entry((ref_pos, strand_char)).or_default().0 += 1;
                                }
                                b'z' | b'x' | b'h' => {
                                    counts.entry((ref_pos, strand_char)).or_default().1 += 1;
                                }
                                _ => {}
                            }
                        }
                        read_pos += 1;
                        ref_pos += 1;
                    }
                }
                Ins(_) | SoftClip(_) => {
                    read_pos += op.len() as usize;
                }
                Del(_) | RefSkip(_) => {
                    ref_pos += len;
                }
                HardClip(_) | Pad(_) => {}
            }
        }
    }
    Ok(counts)
}

/// Use `modkit extract calls` to collect per-position methylation counts
/// from MM/ML tags. Returns (ref_pos, strand) → (methylated, unmethylated).
fn count_modkit_per_position(
    bam_path: &Path,
    temp_dir: &Path,
) -> Result<HashMap<(u32, char), (u32, u32)>> {
    let extract_out = temp_dir.join("modkit_extract.tsv");
    let output = Command::new("modkit")
        .args(["extract", "calls", "--no-filtering"])
        .arg(bam_path)
        .arg(&extract_out)
        .output()
        .wrap_err("Failed to run modkit extract")?;

    ensure!(
        output.status.success(),
        "modkit extract failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let extract_text = std::fs::read_to_string(&extract_out)?;
    let mut counts: HashMap<(u32, char), (u32, u32)> = HashMap::new();

    // modkit extract calls columns (tab-separated):
    // 0:read_id 1:forward_read_position 2:ref_position 3:chrom
    // 4:mod_strand ... 13:call_code ... 21:flag
    for line in extract_text.lines() {
        if line.starts_with("read_id") {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 22 {
            continue;
        }
        let ref_pos: u32 = match fields[2].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let mod_strand = fields[4].chars().next().unwrap_or('?');
        let call_code = fields[13];

        let entry = counts.entry((ref_pos, mod_strand)).or_default();
        match call_code {
            "m" => entry.0 += 1,
            "-" => entry.1 += 1,
            _ => {}
        }
    }
    Ok(counts)
}

// ---------------------------------------------------------------------------
// Tier 1: Basic tool validation
// ---------------------------------------------------------------------------

#[test]
fn samtools_tag_presence() -> Result<()> {
    require_tool!("samtools");
    let bams = setup_test_bams()?;

    let total_standard = samtools_count_total(&bams.standard)?;
    let mm_count = samtools_count_with_tag(&bams.standard, "MM")?;
    let ml_count = samtools_count_with_tag(&bams.standard, "ML")?;

    // MM/ML are only written for reads with methylation evidence — not all reads.
    // Both tags must be present together and on the same reads.
    ensure!(mm_count > 0, "No standard-mode reads have MM tags");
    ensure!(mm_count < total_standard, "Expected some reads without MM (no methylation evidence)");
    assert_eq!(mm_count, ml_count, "MM and ML tag counts must match ({mm_count} vs {ml_count})");

    let total_legacy = samtools_count_total(&bams.legacy)?;
    let xm_count = samtools_count_with_tag(&bams.legacy, "XM")?;
    assert_eq!(
        xm_count, total_legacy,
        "Not all legacy-mode reads have XM tags ({xm_count}/{total_legacy})"
    );

    let xr_count = samtools_count_with_tag(&bams.legacy, "XR")?;
    assert_eq!(
        xr_count, total_legacy,
        "Not all legacy-mode reads have XR tags ({xr_count}/{total_legacy})"
    );

    let xg_count = samtools_count_with_tag(&bams.legacy, "XG")?;
    assert_eq!(
        xg_count, total_legacy,
        "Not all legacy-mode reads have XG tags ({xg_count}/{total_legacy})"
    );

    assert_eq!(
        total_standard, total_legacy,
        "Record count mismatch between standard ({total_standard}) and legacy ({total_legacy})"
    );

    ensure!(total_standard > 0, "No records in output BAMs");

    Ok(())
}

#[test]
fn modkit_summary_validates_standard_bam() -> Result<()> {
    require_tool!("modkit");
    let bams = setup_test_bams()?;

    let output = Command::new("modkit")
        .args(["summary", "--no-sampling"])
        .arg(&bams.standard)
        .output()
        .wrap_err("Failed to run modkit")?;

    ensure!(
        output.status.success(),
        "modkit summary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8_lossy(&output.stdout);

    ensure!(!stdout.trim().is_empty(), "modkit summary produced no output");

    let has_modifications = stdout.lines().any(|line| {
        line.contains('m')
            && line.split_whitespace().any(|tok| tok.parse::<u64>().is_ok_and(|n| n > 0))
    });

    ensure!(has_modifications, "modkit summary reports no 5mC modifications. Output:\n{stdout}");

    Ok(())
}

#[test]
fn bismark_extractor_reads_legacy_bam() -> Result<()> {
    require_tool!("bismark_methylation_extractor");
    let bams = setup_test_bams()?;
    let output_dir = bams._temp_dir.path().join("bismark_out");
    std::fs::create_dir_all(&output_dir)?;

    let output = Command::new("bismark_methylation_extractor")
        .args(["--comprehensive", "--single-end", "--output"])
        .arg(&output_dir)
        .arg(&bams.legacy)
        .output()
        .wrap_err("Failed to run bismark_methylation_extractor")?;

    ensure!(
        output.status.success(),
        "bismark_methylation_extractor failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let cpg_files: Vec<_> = std::fs::read_dir(&output_dir)?
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().contains("CpG"))
        .collect();

    ensure!(!cpg_files.is_empty(), "No CpG output files from bismark");

    for file in &cpg_files {
        let size = file.metadata()?.len();
        ensure!(size > 0, "CpG output file is empty: {:?}", file.path());
    }

    Ok(())
}

#[test]
fn xm_annotations_only_at_reference_cpgs() -> Result<()> {
    use rust_htslib::bam::{self, Read as BamRead};
    use rust_htslib::faidx;

    let bams = setup_test_bams()?;
    let fasta =
        faidx::Reader::from_path("tests/data/test.fasta.gz").wrap_err("Failed to open ref")?;

    let mut reader = bam::Reader::from_path(&bams.legacy).wrap_err("open legacy BAM")?;
    let mut record = bam::Record::new();
    let mut records_checked = 0u32;

    while let Some(result) = reader.read(&mut record) {
        result?;
        records_checked += 1;

        let xm = match record.aux(b"XM") {
            Ok(rust_htslib::bam::record::Aux::String(s)) => s.to_string(),
            _ => bail!("Missing/bad XM tag at record {records_checked}"),
        };

        let xm_bytes = xm.as_bytes();
        let flag = record.flags();

        let is_ob = match record.aux(b"XG") {
            Ok(rust_htslib::bam::record::Aux::String(s)) => s == "GA",
            _ => bail!("Missing/bad XG tag at record {records_checked}"),
        };

        let cigar = record.cigar();
        let mut read_pos: usize = 0;
        let mut ref_pos = record.pos() as usize;

        for op in cigar.iter() {
            use rust_htslib::bam::record::Cigar::*;
            let len = op.len() as usize;
            match op {
                Match(_) | Equal(_) | Diff(_) => {
                    for _ in 0..len {
                        if read_pos < xm_bytes.len() {
                            let ch = xm_bytes.get(read_pos).copied().unwrap_or(b'.');
                            if ch == b'z' || ch == b'Z' {
                                let is_cpg = if is_ob {
                                    ref_pos > 0
                                        && ref_base_at(&fasta, "chr19", ref_pos.saturating_sub(1))
                                            .is_some_and(|b| b == b'C')
                                        && ref_base_at(&fasta, "chr19", ref_pos)
                                            .is_some_and(|b| b == b'G')
                                } else {
                                    ref_base_at(&fasta, "chr19", ref_pos).is_some_and(|b| b == b'C')
                                        && ref_base_at(&fasta, "chr19", ref_pos + 1)
                                            .is_some_and(|b| b == b'G')
                                };

                                ensure!(
                                    is_cpg,
                                    "z/Z annotation at read_pos={read_pos} ref_pos={ref_pos} \
                                     is not a reference CpG (flag={flag}, record #{records_checked})"
                                );
                            }
                        }
                        read_pos += 1;
                        ref_pos += 1;
                    }
                }
                Ins(_) | SoftClip(_) => {
                    read_pos += len;
                }
                Del(_) | RefSkip(_) => {
                    ref_pos += len;
                }
                HardClip(_) | Pad(_) => {}
            }
        }
    }

    ensure!(records_checked > 0, "No records were checked");
    eprintln!("Verified {records_checked} records: all z/Z annotations are at reference CpGs");

    Ok(())
}

// ---------------------------------------------------------------------------
// Tier 2: Cross-validation tests
// ---------------------------------------------------------------------------

/// For each CpG position in the calls BED, count Z/z reads in the XM tag
/// and verify methylation fractions match. The rewritten BAM may contain a
/// subset of reads (those starting within the region), so we compare
/// fractions rather than exact counts.
#[test]
fn xm_counts_match_calls_bed() -> Result<()> {
    let bams = setup_test_bams()?;
    let bed_calls = parse_calls_bed(&bams.calls_bed)?;
    ensure!(!bed_calls.is_empty(), "No calls in BED file");

    let xm_counts = count_xm_per_position(&bams.legacy)?;

    let mut positions_checked = 0u32;
    let mut mismatches = Vec::new();

    for call in &bed_calls {
        let key = (call.pos, call.strand);
        let (obs_mod, obs_unmod) = xm_counts.get(&key).copied().unwrap_or((0, 0));
        let obs_total = obs_mod + obs_unmod;
        let bed_total = call.mod_count + call.unmod_count;

        // Skip positions with no XM coverage (reads may start outside region)
        if obs_total == 0 || bed_total == 0 {
            continue;
        }

        let bed_frac = call.mod_count as f64 / bed_total as f64;
        let obs_frac = obs_mod as f64 / obs_total as f64;
        let diff = (bed_frac - obs_frac).abs();

        // Fractions should match exactly since both come from the same reads.
        // Use small tolerance for rounding in very low coverage positions.
        if diff > 0.01 {
            mismatches.push(format!(
                "pos={} strand={}: BED beta={bed_frac:.3} ({}/{}), \
                 XM beta={obs_frac:.3} ({obs_mod}/{obs_total})",
                call.pos, call.strand, call.mod_count, bed_total,
            ));
        }
        positions_checked += 1;
    }

    ensure!(
        mismatches.is_empty(),
        "XM methylation fraction mismatches vs BED:\n{}",
        mismatches.join("\n")
    );

    ensure!(positions_checked > 0, "No CpG positions checked");
    eprintln!("Verified {positions_checked} CpG positions: XM fractions match BED calls");

    Ok(())
}

/// Use `modkit extract calls` to get per-read modification calls from the
/// standard BAM, then compare per-position fractions against the calls BED.
///
/// MM/ML tags only encode modifications at C bases in the stored SEQ, so for
/// paired reads overlapping a CpG, only one mate contributes (the other has G
/// at that position). This means modkit sees roughly half the reads that BED
/// does. We compare fractions, not exact counts, and require a minimum
/// coverage in modkit to avoid noise from very low coverage positions.
#[test]
fn modkit_extract_matches_calls_bed() -> Result<()> {
    require_tool!("modkit");
    let bams = setup_test_bams()?;
    let bed_calls = parse_calls_bed(&bams.calls_bed)?;

    let modkit_counts = count_modkit_per_position(&bams.standard, bams._temp_dir.path())?;

    let mut positions_checked = 0u32;
    let mut mismatches = Vec::new();

    for call in &bed_calls {
        let key = (call.pos, call.strand);
        let (mk_mod, mk_unmod) = modkit_counts.get(&key).copied().unwrap_or((0, 0));
        let mk_total = mk_mod + mk_unmod;
        let bed_total = call.mod_count + call.unmod_count;

        // Require minimum coverage in modkit to avoid noise from paired-read
        // asymmetry at very low coverage positions
        if mk_total < 5 || bed_total == 0 {
            continue;
        }

        let bed_frac = call.mod_count as f64 / bed_total as f64;
        let mk_frac = mk_mod as f64 / mk_total as f64;
        let diff = (bed_frac - mk_frac).abs();

        if diff > 0.05 {
            mismatches.push(format!(
                "pos={} strand={}: BED beta={bed_frac:.3} ({}/{}), \
                 modkit beta={mk_frac:.3} ({mk_mod}/{mk_total})",
                call.pos, call.strand, call.mod_count, bed_total,
            ));
        }
        positions_checked += 1;
    }

    ensure!(
        mismatches.is_empty(),
        "Methylation fraction mismatches between modkit and BED:\n{}",
        mismatches.join("\n")
    );

    ensure!(positions_checked > 0, "No CpG positions checked against modkit");
    eprintln!(
        "Verified {positions_checked} CpG positions: modkit extract fractions match BED calls"
    );

    Ok(())
}

/// Parse bismark's CpG extraction output and verify per-position methylation
/// fractions match the calls BED.
#[test]
fn bismark_counts_match_calls_bed() -> Result<()> {
    require_tool!("bismark_methylation_extractor");
    let bams = setup_test_bams()?;
    let bed_calls = parse_calls_bed(&bams.calls_bed)?;

    let output_dir = bams._temp_dir.path().join("bismark_xval");
    std::fs::create_dir_all(&output_dir)?;

    let output = Command::new("bismark_methylation_extractor")
        .args(["--comprehensive", "--single-end", "--output"])
        .arg(&output_dir)
        .arg(&bams.legacy)
        .output()
        .wrap_err("Failed to run bismark_methylation_extractor")?;

    ensure!(
        output.status.success(),
        "bismark_methylation_extractor failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Bismark CpG output: read_id \t methylation_state(+/-) \t chromosome \t position(1-based) \t context
    let cpg_file = std::fs::read_dir(&output_dir)?
        .filter_map(|e| e.ok())
        .find(|e| e.file_name().to_string_lossy().contains("CpG"))
        .ok_or_else(|| eyre!("No CpG output file from bismark"))?;

    let cpg_text = std::fs::read_to_string(cpg_file.path())?;

    // Aggregate: (0-based position) → (methylated, unmethylated)
    let mut bismark_counts: HashMap<u32, (u32, u32)> = HashMap::new();

    for line in cpg_text.lines() {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            continue;
        }
        let meth_state = fields[1];
        let pos_1based: u32 = match fields[3].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let pos_0based = pos_1based - 1;

        let entry = bismark_counts.entry(pos_0based).or_default();
        match meth_state {
            "+" => entry.0 += 1,
            "-" => entry.1 += 1,
            _ => {}
        }
    }

    let mut positions_checked = 0u32;
    let mut mismatches = Vec::new();

    // Bismark reports the C position in reference coords. For + strand (OT),
    // this is the C of the CpG. For - strand (OB), this is the reference G
    // position (where the bottom-strand C lives). Both match the BED position
    // directly: BED + strand = C position, BED - strand = G position.
    for call in &bed_calls {
        let bismark_pos = call.pos;

        let (bk_mod, bk_unmod) = bismark_counts.get(&bismark_pos).copied().unwrap_or((0, 0));
        let bk_total = bk_mod + bk_unmod;
        let bed_total = call.mod_count + call.unmod_count;

        if bk_total == 0 || bed_total == 0 {
            continue;
        }

        let bed_frac = call.mod_count as f64 / bed_total as f64;
        let bk_frac = bk_mod as f64 / bk_total as f64;
        let diff = (bed_frac - bk_frac).abs();

        if diff > 0.01 {
            mismatches.push(format!(
                "pos={} strand={}: BED beta={bed_frac:.3} ({}/{}), \
                 bismark beta={bk_frac:.3} ({bk_mod}/{bk_total})",
                call.pos, call.strand, call.mod_count, bed_total,
            ));
        }
        positions_checked += 1;
    }

    ensure!(
        mismatches.is_empty(),
        "Methylation fraction mismatches between bismark and BED:\n{}",
        mismatches.join("\n")
    );

    ensure!(positions_checked > 0, "No CpG positions checked against bismark");
    eprintln!("Verified {positions_checked} CpG positions: bismark extraction matches BED calls");

    Ok(())
}

/// For each CpG position, verify that methylation fractions from XM tags
/// (legacy) match those from MM/ML tags (standard) via modkit.
///
/// MM/ML tags only encode modifications at C bases in the stored SEQ, so for
/// paired reads overlapping a CpG, only one mate contributes. XM annotates
/// both mates. This means exact counts will differ (typically 2:1), but
/// fractions should agree at positions with adequate coverage.
#[test]
fn legacy_standard_per_position_agreement() -> Result<()> {
    require_tool!("modkit");
    let bams = setup_test_bams()?;

    let xm_counts = count_xm_per_position(&bams.legacy)?;
    let mm_counts = count_modkit_per_position(&bams.standard, bams._temp_dir.path())?;

    let mut positions_checked = 0u32;
    let mut mismatches = Vec::new();

    for (key, (xm_mod, xm_unmod)) in &xm_counts {
        let xm_total = xm_mod + xm_unmod;
        let (mm_mod, mm_unmod) = mm_counts.get(key).copied().unwrap_or((0, 0));
        let mm_total = mm_mod + mm_unmod;

        // Require minimum coverage to avoid noise from paired-read asymmetry
        if xm_total < 5 || mm_total < 5 {
            continue;
        }

        let xm_frac = *xm_mod as f64 / xm_total as f64;
        let mm_frac = mm_mod as f64 / mm_total as f64;
        let diff = (xm_frac - mm_frac).abs();

        if diff > 0.05 {
            mismatches.push(format!(
                "pos={} strand={}: XM beta={xm_frac:.3} ({xm_mod}/{xm_total}), \
                 modkit beta={mm_frac:.3} ({mm_mod}/{mm_total})",
                key.0, key.1,
            ));
        }
        positions_checked += 1;
    }

    ensure!(
        mismatches.is_empty(),
        "Legacy/standard per-position methylation fraction mismatches:\n{}",
        mismatches.join("\n")
    );

    ensure!(positions_checked > 0, "No positions checked");
    eprintln!(
        "Verified {positions_checked} CpG positions: \
         legacy XM and standard MM/ML fractions agree"
    );

    Ok(())
}

/// Verify that modbedtools bam2mod does not crash on our standard modBAM output.
///
/// Empty ML:B:C arrays (written for reads with no methylated positions) cause modbedtools
/// to segfault. This test ensures we no longer produce empty ML arrays.
#[test]
fn modbedtools_bam2mod_does_not_crash() -> Result<()> {
    require_tool!("modbedtools");
    let bams = setup_test_bams()?;
    let temp_dir = bams._temp_dir.path().join("modbedtools_out");
    std::fs::create_dir_all(&temp_dir)?;

    let output_prefix = temp_dir.join("out");

    let output = Command::new("modbedtools")
        .args(["bam2mod", "-o"])
        .arg(&output_prefix)
        .arg(&bams.standard)
        .output()
        .wrap_err("Failed to run modbedtools bam2mod")?;

    ensure!(
        output.status.success(),
        "modbedtools bam2mod failed (exit {}): {}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );

    // At least one output modbed file should be non-empty
    let has_output = ["C", "G"].iter().any(|base| {
        let path = temp_dir.join(format!("out.{base}.modbed"));
        path.exists() && path.metadata().is_ok_and(|m| m.len() > 0)
    });
    ensure!(has_output, "modbedtools produced no modbed output");

    Ok(())
}

fn ref_base_at(fasta: &rust_htslib::faidx::Reader, chrom: &str, pos: usize) -> Option<u8> {
    let seq = fasta.fetch_seq_string(chrom, pos, pos).ok()?;
    seq.bytes().next().map(|b| b.to_ascii_uppercase())
}
