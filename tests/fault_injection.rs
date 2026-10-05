//! How the pipelines (`call`, `per-read`, `bam`) behave when a thread crashes,
//! an error occurs or the process is killed, using the `RASTAIR_INJECT_*`
//! variables.
//!
//! Fault injection is compiled out of release builds, so these tests only run
//! in debug builds.
#![cfg(unix)]
#![cfg_attr(
    not(debug_assertions),
    allow(dead_code, unused_imports, reason = "the tests are ignored in release builds")
)]

mod utils;
use utils::{faults::*, *};

// ---------------------------------------------------------------------------
// call

#[test]
fn successful_run_moves_outputs_into_place() -> Result<()> {
    let dir = TempDir::new()?;
    let (vcf, bed) = (dir.path().join("out.vcf.gz"), dir.path().join("out.bed.gz"));
    ensure_success(&call(&vcf, &[NO_ML, "--bed", &bed.to_string_lossy()], &[])?)?;

    for file in ["out.vcf.gz", "out.vcf.gz.csi", "out.bed.gz", "out.bed.gz.tbi"] {
        ensure!(dir.path().join(file).exists(), "{file} should exist");
    }
    let leftovers: Vec<_> = std::fs::read_dir(dir.path())?
        .map(|entry| Ok(entry?.file_name().to_string_lossy().into_owned()))
        .filter(|name: &Result<String>| name.as_ref().map_or(true, |name| name.contains("partial")))
        .collect::<Result<_>>()?;
    ensure!(leftovers.is_empty(), "partial files left behind: {leftovers:?}");
    ensure!(!records(&vcf)?.is_empty() && !records(&bed)?.is_empty());
    Ok(())
}

#[test]
#[cfg_attr(not(debug_assertions), ignore = "fault injection only exists in debug builds")]
fn worker_panic_stops_the_run_and_leaves_a_valid_prefix() -> Result<()> {
    let dir = TempDir::new()?;
    let (full, crashed) = (dir.path().join("full.vcf"), dir.path().join("crashed.vcf"));
    ensure_success(&call(&full, &[NO_ML], &[])?)?;

    let output = call(&crashed, &[NO_ML], &[(INJECT_PANIC, "worker@6")])?;
    let stderr = output.stderr();
    ensure!(!output.status.success(), "should fail after a worker panic");
    ensure!(stderr.contains("A worker panicked (see the crash report above)"), "{stderr}");
    ensure_crash_reported(&stderr, "src/call.rs")?;
    ensure!(stderr.contains("Output is incomplete"), "{stderr}");
    ensure!(stderr.contains("VCF output left behind under its partial name"), "{stderr}");
    // URL-encoded: the version, and the span of the segment that panicked
    let link = issue_link(&stderr)?;
    ensure!(link.contains(&format!("%7C{}%7C", env!("CARGO_PKG_VERSION"))), "{link}");
    ensure!(link.contains("process_region_wrapper%0A+++++++++++with+region%3Dchr19"), "{link}");

    let crashed = only_partial_output(&crashed)?;
    ensure_truncated_prefix(&records(&full)?, &records(&crashed)?)
}

#[test]
#[cfg_attr(not(debug_assertions), ignore = "fault injection only exists in debug builds")]
fn writer_panic_stops_the_run() -> Result<()> {
    let dir = TempDir::new()?;
    let out = dir.path().join("out.vcf");
    let output = call(&out, &[NO_ML], &[(INJECT_PANIC, "writer@3")])?;
    let stderr = output.stderr();
    ensure!(!output.status.success(), "should fail after a writer panic");
    ensure!(
        stderr.contains("Writing the output panicked (see the crash report above)"),
        "{stderr}"
    );
    ensure_crash_reported(&stderr, "src/call/writer.rs")?;
    only_partial_output(&out)?;
    Ok(())
}

#[test]
#[cfg_attr(not(debug_assertions), ignore = "fault injection only exists in debug builds")]
fn writer_error_stops_the_run() -> Result<()> {
    let dir = TempDir::new()?;
    let out = dir.path().join("out.vcf");
    let output = call(&out, &[NO_ML], &[(INJECT_ERROR, "writer@3")])?;
    let stderr = output.stderr();
    ensure!(!output.status.success(), "should fail after a writer error");
    ensure!(stderr.contains("Failed to write the output"), "{stderr}");
    ensure!(stderr.contains("Injected error at writer"), "{stderr}");
    only_partial_output(&out)?;
    Ok(())
}

#[test]
#[cfg_attr(not(debug_assertions), ignore = "fault injection only exists in debug builds")]
fn failed_segment_is_skipped_and_fails_the_run() -> Result<()> {
    let dir = TempDir::new()?;
    let (full, lossy) = (dir.path().join("full.vcf"), dir.path().join("lossy.vcf"));
    ensure_success(&call(&full, &[NO_ML], &[])?)?;

    let output = call(&lossy, &[NO_ML], &[(INJECT_ERROR, "worker@6")])?;
    let stderr = output.stderr();
    ensure!(output.status.code() == Some(1), "should fail after a segment error: {stderr}");
    ensure!(stderr.contains("Failed to process segment, its output is missing"), "{stderr}");
    ensure!(stderr.contains("Injected error at worker"), "{stderr}");
    ensure!(stderr.contains("1 of 20 segments failed to process"), "{stderr}");
    ensure!(issue_link(&stderr).is_err(), "not a bug: {stderr}");

    let lossy = only_partial_output(&lossy)?;
    ensure_missing_one_stretch(&records(&full)?, &records(&lossy)?)
}

#[test]
#[cfg_attr(not(debug_assertions), ignore = "fault injection only exists in debug builds")]
fn a_killed_run_leaves_only_partial_output() -> Result<()> {
    let dir = TempDir::new()?;
    for (name, signal) in [("TERM", libc::SIGTERM), ("INT", libc::SIGINT)] {
        let out = dir.path().join(format!("{name}.vcf.gz"));
        let spec = format!("{name}:worker@6");
        let output = call(&out, &[NO_ML], &[(INJECT_SIGNAL, &spec)])?;
        ensure_killed_by(&output, signal)?;
        only_partial_output(&out)?;
    }
    Ok(())
}

#[test]
#[cfg_attr(not(debug_assertions), ignore = "fault injection only exists in debug builds")]
fn gpu_thread_panic_falls_back_to_the_cpu() -> Result<()> {
    let dir = TempDir::new()?;
    let (cpu, crashed) = (dir.path().join("cpu.vcf"), dir.path().join("crashed.vcf"));

    // More workers than the GPU thread keeps up with, so jobs are queued when it
    // crashes: those used to wait forever.
    let threads = ["--threads", "4"];
    let fault = &[(INJECT_PANIC, "gpu-dispatch@2")];
    let output = call(&crashed, &[&threads[..], &["--gpu"]].concat(), fault)?;
    let stderr = output.stderr();
    if stderr.contains("Failed to initialise GPU context") {
        eprintln!("Skipping: no GPU available");
        return Ok(());
    }
    ensure!(output.status.success(), "should recover from a GPU thread panic: {stderr}");
    ensure!(stderr.contains("A thread panicked and recovered"), "{stderr}");
    ensure!(stderr.contains("scoring the rest of the run on the CPU"), "{stderr}");
    ensure!(!stderr.contains("The application panicked"), "no crash report: {stderr}");

    // GPU and CPU scores may differ in the last digits, so compare record count
    ensure_success(&call(&cpu, &threads, &[])?)?;
    ensure!(records(&cpu)?.len() == records(&crashed)?.len(), "output should be complete");
    Ok(())
}

#[test]
#[cfg_attr(not(debug_assertions), ignore = "fault injection only exists in debug builds")]
fn panic_report_shows_where_the_panic_happened() -> Result<()> {
    let dir = TempDir::new()?;
    let mut cmd = rastair();
    cmd.args(CALL_TEST_BAM).args([REGION, SMALL_SEGMENTS, NO_ML]);
    cmd.arg("--vcf").arg(dir.path().join("out.vcf"));
    with_faults(&mut cmd, &[(INJECT_PANIC, "worker@2")]).env("RUST_BACKTRACE", "1");
    let output = run_with_timeout(&mut cmd)?;
    let stderr = output.stderr();
    ensure!(!output.status.success(), "should fail after a worker panic");

    // Captured where the panic happened, not where it was caught
    let backtrace =
        stderr.split_once("Backtrace:").ok_or_else(|| eyre!("no panic backtrace in: {stderr}"))?.1;
    ensure!(backtrace.contains("rastair::call::process_segment"), "{stderr}");
    // Where the panic was caught is the same few lines every time: not shown
    ensure!(!stderr.contains("BACKTRACE"), "{stderr}");
    // Shown, but not in the issue link: too long for GitHub with it
    let link = issue_link(&stderr)?;
    ensure!(!link.contains("Backtrace") && link.len() < 8000, "{stderr}");
    Ok(())
}

// ---------------------------------------------------------------------------
// per-read

#[test]
#[cfg_attr(not(debug_assertions), ignore = "fault injection only exists in debug builds")]
fn per_read_worker_panic_leaves_a_valid_prefix() -> Result<()> {
    let dir = TempDir::new()?;
    let (full, crashed) = (dir.path().join("full.bed.gz"), dir.path().join("crashed.bed.gz"));
    ensure_success(&per_read(&full, &["--threads=2"], &[])?)?;
    ensure!(dir.path().join("full.bed.gz.tbi").exists(), "index moved into place");

    let output = per_read(&crashed, &["--threads=2"], &[(INJECT_PANIC, "worker@6")])?;
    let stderr = output.stderr();
    ensure!(!output.status.success(), "should fail after a worker panic");
    ensure_crash_reported(&stderr, "src/call_reads.rs")?;

    let crashed = only_partial_output(&crashed)?;
    ensure_truncated_prefix(&records(&full)?, &records(&crashed)?)
}

#[test]
#[cfg_attr(not(debug_assertions), ignore = "fault injection only exists in debug builds")]
fn per_read_writer_error_stops_the_run() -> Result<()> {
    let dir = TempDir::new()?;
    let out = dir.path().join("out.bed");
    let output = per_read(&out, &["--threads=2"], &[(INJECT_ERROR, "writer@3")])?;
    let stderr = output.stderr();
    ensure!(!output.status.success(), "should fail after a writer error");
    ensure!(stderr.contains("Injected error at writer"), "{stderr}");
    only_partial_output(&out)?;
    Ok(())
}

#[test]
#[cfg_attr(not(debug_assertions), ignore = "fault injection only exists in debug builds")]
fn per_read_failed_segment_is_skipped_and_fails_the_run() -> Result<()> {
    let dir = TempDir::new()?;
    let (full, lossy) = (dir.path().join("full.bed"), dir.path().join("lossy.bed"));
    ensure_success(&per_read(&full, &["--threads=2"], &[])?)?;

    let output = per_read(&lossy, &["--threads=2"], &[(INJECT_ERROR, "worker@6")])?;
    let stderr = output.stderr();
    ensure!(output.status.code() == Some(1), "should fail after a segment error: {stderr}");
    ensure!(stderr.contains("1 of 20 segments failed to process"), "{stderr}");

    let lossy = only_partial_output(&lossy)?;
    ensure_missing_one_stretch(&records(&full)?, &records(&lossy)?)
}

// ---------------------------------------------------------------------------
// bam

#[test]
#[cfg_attr(not(debug_assertions), ignore = "fault injection only exists in debug builds")]
fn bam_worker_panic_leaves_a_valid_prefix() -> Result<()> {
    let dir = TempDir::new()?;
    let calls = bam_calls(dir.path())?;
    let (full, crashed) = (dir.path().join("full.bam"), dir.path().join("crashed.bam"));
    ensure_success(&bam(&calls, &full, &["--threads=2"], &[])?)?;

    let output = bam(&calls, &crashed, &["--threads=2"], &[(INJECT_PANIC, "worker@6")])?;
    let stderr = output.stderr();
    ensure!(!output.status.success(), "should fail after a worker panic");
    ensure_crash_reported(&stderr, "src/bam.rs")?;

    let crashed = only_partial_output(&crashed)?;
    ensure_truncated_prefix(&bam_records(&full)?, &bam_records(&crashed)?)
}

/// The pre-filled GitHub issue link in a crash report.
fn issue_link(stderr: &str) -> Result<&str> {
    stderr
        .split_once("Consider reporting this error using this URL: ")
        .and_then(|(_, rest)| rest.split_whitespace().next())
        .ok_or_else(|| eyre!("no issue link"))
}
