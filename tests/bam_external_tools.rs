#![cfg(feature = "external-tool-tests")]

mod utils;
use std::path::{Path, PathBuf};
use utils::*;

struct TestBams {
    _temp_dir: TempDir,
    legacy: PathBuf,
    standard: PathBuf,
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

    let region = "chr19:6103075-6103200";

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

    Ok(TestBams { _temp_dir: temp_dir, legacy: legacy_bam, standard: standard_bam })
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

// ---------------------------------------------------------------------------
// Tier 1 tests
// ---------------------------------------------------------------------------

#[test]
fn samtools_tag_presence() -> Result<()> {
    require_tool!("samtools");
    let bams = setup_test_bams()?;

    let total_standard = samtools_count_total(&bams.standard)?;
    let mm_count = samtools_count_with_tag(&bams.standard, "MM")?;
    assert_eq!(
        mm_count, total_standard,
        "Not all standard-mode reads have MM tags ({mm_count}/{total_standard})"
    );

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

    // modkit summary prints a table; look for a row with modification code "m"
    // (5mC) and a nonzero count. The exact format varies by version, so just
    // check that 'm' appears somewhere in the output alongside a number > 0.
    let has_modifications = stdout.lines().any(|line| {
        // Lines with modification data contain the code letter and numeric counts
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

    // Use --single-end because our regional BAM doesn't contain all mates,
    // which makes bismark's paired-end mode fail. Single-end mode still
    // fully validates the XM/XR/XG tags on each read independently.
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

// ---------------------------------------------------------------------------
// Internal validation (no external tool needed, but feature-gated alongside
// the external tests since it validates the same properties)
// ---------------------------------------------------------------------------

#[test]
fn xm_annotations_only_at_reference_cpgs() -> Result<()> {
    use rust_htslib::bam::{self, Read as BamRead};
    use rust_htslib::faidx;

    let bams = setup_test_bams()?;
    let fasta = faidx::Reader::from_path("tests/data/test.fasta.gz")
        .wrap_err("Failed to open reference FASTA")?;

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

        // Determine methylation strand from XG tag (CT = OT, GA = OB)
        let is_ob = match record.aux(b"XG") {
            Ok(rust_htslib::bam::record::Aux::String(s)) => s == "GA",
            _ => bail!("Missing/bad XG tag at record {records_checked}"),
        };

        // Walk through aligned pairs to check z/Z annotations
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
                                // This position is annotated as CpG -- verify reference
                                let is_cpg = if is_ob {
                                    // OB strand: G in CpG context means ref[pos-1] == 'C'
                                    ref_pos > 0
                                        && ref_base_at(&fasta, "chr19", ref_pos.saturating_sub(1))
                                            .is_some_and(|b| b == b'C')
                                        && ref_base_at(&fasta, "chr19", ref_pos)
                                            .is_some_and(|b| b == b'G')
                                } else {
                                    // OT strand: C in CpG context means ref[pos+1] == 'G'
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

fn ref_base_at(fasta: &rust_htslib::faidx::Reader, chrom: &str, pos: usize) -> Option<u8> {
    let seq = fasta.fetch_seq_string(chrom, pos, pos).ok()?;
    seq.bytes().next().map(|b| b.to_ascii_uppercase())
}
