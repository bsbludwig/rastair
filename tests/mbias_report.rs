#![cfg(feature = "external-tool-tests")]

mod utils;

use std::path::Path;
use std::process::{Command, Stdio};
use utils::*;

fn tool_is_available(name: &str) -> bool {
    Command::new(name)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// The QC report only renders if R plus the packages the M-bias path loads are
/// installed. Without them the test self-skips rather than failing.
fn r_report_toolchain_available() -> bool {
    if !tool_is_available("tabix") || !tool_is_available("bgzip") {
        return false;
    }
    Command::new("Rscript")
        .args(["-e", "suppressMessages({library(rmarkdown); library(argparser); library(data.table); library(ggplot2)})"])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

macro_rules! require_r_report_toolchain {
    () => {
        if !r_report_toolchain_available() {
            eprintln!("Skipping test: Rscript (with rmarkdown/data.table/ggplot2), tabix or bgzip not available");
            return Ok(());
        }
    };
}

/// Per-read BED header expected by the QC report (mirrors
/// `PerRead::HEADER` in `src/bed/per_read/format.rs`).
const PER_READ_HEADER: &str = "#chr\tstart\tend\tread_id\tmapq\torientation\tinsert_size\tread_length\tflag\tnum_cpg\tnum_mod\tmod_cpgs\tunmod_cpgs\tsnp_cpgs\tmod_denovos\tunmod_denovos";

/// A per-read BED with one healthy contig (`chr_ok`, several First/OT reads
/// spanning >=3 distinct read positions) and one sparse contig (`chr_sparse`,
/// a single read with one CpG). `flag=99` is paired+proper+mate-reverse+first,
/// i.e. First/OT, and passes the default include=3 / exclude=3852 filters.
/// Coordinate-sorted so `tabix -p bed` can index it.
fn fixture_bed() -> String {
    let rows = [
        // chr_ok: positions {10,20,30,40,50,60} in the First/OT group
        "chr_ok\t100\t200\tok1\t60\t+\t150\t100\t99\t3\t3\t10,30,50\t\t\t\t",
        "chr_ok\t200\t300\tok2\t60\t+\t150\t100\t99\t3\t0\t\t10,30,50\t\t\t",
        "chr_ok\t300\t400\tok3\t60\t+\t150\t100\t99\t3\t1\t20\t40,60\t\t\t",
        // chr_sparse: a single read with one CpG -> one position -> skipped
        "chr_sparse\t100\t200\tsp1\t60\t+\t150\t100\t99\t1\t1\t10\t\t\t\t",
    ];
    let mut bed = String::from(PER_READ_HEADER);
    bed.push('\n');
    for row in rows {
        bed.push_str(row);
        bed.push('\n');
    }
    bed
}

fn bgzip(bed_path: &Path) -> Result<std::path::PathBuf> {
    let gz_path = bed_path.with_extension("bed.gz");
    let output = Command::new("bgzip").arg("-c").arg(bed_path).output().wrap_err("run bgzip")?;
    ensure!(output.status.success(), "bgzip failed with status: {}", output.status);
    std::fs::write(&gz_path, output.stdout).wrap_err("write bgzipped BED")?;
    Ok(gz_path)
}

/// A sparse contig must not abort the whole report: it is skipped (no plot, no
/// cutoffs file) while healthy contigs are still processed.
#[test]
fn mbias_report_skips_sparse_contig_instead_of_aborting() -> Result<()> {
    require_r_report_toolchain!();

    let temp_dir = TempDir::new()?;
    let bed_path = temp_dir.path().join("reads.bed");
    std::fs::write(&bed_path, fixture_bed())?;
    let bed_gz = bgzip(&bed_path)?;

    let out_dir = temp_dir.path().join("out");
    std::fs::create_dir(&out_dir)?;

    rastair()
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .args(["mbias", "--no-vbias", "--no-gc"])
        .arg("--bed")
        .arg(&bed_gz)
        .arg("--output-prefix")
        .arg(&out_dir)
        .succeeds()
        .wrap_err("mbias render should succeed despite the sparse contig")?;

    ensure!(out_dir.join("qc_report.html").exists(), "qc_report.html was not produced");
    ensure!(
        out_dir.join("chr_ok_cutoffs.txt").exists(),
        "healthy contig chr_ok should still get a cutoffs file"
    );
    ensure!(
        !out_dir.join("chr_sparse_cutoffs.txt").exists(),
        "sparse contig chr_sparse should be skipped (no cutoffs file)"
    );

    Ok(())
}

/// When the run is explicitly scoped to a sparse contig, that is an error
/// rather than a silent skip: the user asked for exactly this data.
#[test]
fn mbias_report_errors_when_explicit_region_is_sparse() -> Result<()> {
    require_r_report_toolchain!();

    let temp_dir = TempDir::new()?;
    let bed_path = temp_dir.path().join("reads.bed");
    std::fs::write(&bed_path, fixture_bed())?;
    let bed_gz = bgzip(&bed_path)?;

    let out_dir = temp_dir.path().join("out");
    std::fs::create_dir(&out_dir)?;

    let output = rastair()
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .args(["mbias", "--no-vbias", "--no-gc", "--region", "chr_sparse"])
        .arg("--bed")
        .arg(&bed_gz)
        .arg("--output-prefix")
        .arg(&out_dir)
        .output()
        .wrap_err("run mbias scoped to the sparse contig")?;

    ensure!(
        !output.status.success(),
        "mbias should fail when the explicitly requested contig is too sparse"
    );
    ensure!(
        !out_dir.join("chr_sparse_cutoffs.txt").exists(),
        "no cutoffs file should be written for the failed sparse contig"
    );

    Ok(())
}
