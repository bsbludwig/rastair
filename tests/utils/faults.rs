//! Running the pipelines with injected faults (`RASTAIR_INJECT_*`), and
//! checking what they leave behind.

use super::*;
use std::{
    io::Read as _,
    os::unix::process::ExitStatusExt as _,
    path::PathBuf,
    process::{Output, Stdio},
    time::{Duration, Instant},
};

pub const INJECT_PANIC: &str = "RASTAIR_INJECT_PANIC";
pub const INJECT_ERROR: &str = "RASTAIR_INJECT_ERROR";
pub const INJECT_SIGNAL: &str = "RASTAIR_INJECT_SIGNAL";
/// Its first three segments have no reads, so faults are injected at the sixth
/// segment (`worker@6`) to have records before and after it
pub const REGION: &str = "--region=chr19:6100000-6120000";
/// Many small segments, so a fault lands in the middle of the run
pub const SMALL_SEGMENTS: &str = "--segment-max-length=1000";
/// What [`REGION`] is split into with [`SMALL_SEGMENTS`], in every pipeline
pub const SEGMENTS: usize = 20;
/// A crash must never turn into a hang; fail the test instead of blocking CI
pub const TIMEOUT: Duration = Duration::from_secs(120);
/// Written as the last block of every properly closed BGZF file
pub const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Faults to inject, as `(variable, spec)`
pub type Faults<'a> = &'a [(&'a str, &'a str)];

pub fn with_faults<'c>(cmd: &'c mut Command, faults: Faults<'_>) -> &'c mut Command {
    for var in [INJECT_PANIC, INJECT_ERROR, INJECT_SIGNAL] {
        cmd.env_remove(var);
    }
    for (var, spec) in faults {
        cmd.env(var, spec);
    }
    cmd
}

/// `call` writing VCF to `out`
pub fn call(out: &Path, extra: &[&str], faults: Faults<'_>) -> Result<Output> {
    call_to("--vcf", out, extra, faults)
}

/// `call` writing to `out` with the given output flag (`--vcf`, `--bed`)
pub fn call_to(flag: &str, out: &Path, extra: &[&str], faults: Faults<'_>) -> Result<Output> {
    let mut cmd = rastair();
    cmd.args(CALL_TEST_BAM).args([REGION, SMALL_SEGMENTS, "--all"]).args(extra);
    cmd.arg(flag).arg(out);
    run_with_timeout(with_faults(&mut cmd, faults))
}

pub fn per_read(out: &Path, extra: &[&str], faults: Faults<'_>) -> Result<Output> {
    let mut cmd = rastair();
    cmd.args(["per-read", "--fasta-file=tests/data/test.fasta.gz", "tests/data/test.bam"])
        .args([REGION, SMALL_SEGMENTS])
        .args(extra);
    cmd.arg("--bed").arg(out);
    run_with_timeout(with_faults(&mut cmd, faults))
}

pub fn bam(calls: &Path, out: &Path, extra: &[&str], faults: Faults<'_>) -> Result<Output> {
    let mut cmd = rastair();
    cmd.args(["bam", "legacy", "--fasta-file=tests/data/test.fasta.gz", "tests/data/test.bam"])
        .args([REGION, SMALL_SEGMENTS])
        .args(extra)
        .arg(calls);
    cmd.arg("-o").arg(out);
    run_with_timeout(with_faults(&mut cmd, faults))
}

/// Like `succeeds`, but shows what went wrong
pub fn ensure_success(output: &Output) -> Result<()> {
    ensure!(output.status.success(), "should succeed, got {}: {}", output.status, output.stderr());
    Ok(())
}

pub fn run_with_timeout(cmd: &mut Command) -> Result<Output> {
    let child = cmd
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .wrap_err("Failed to start rastair")?;
    let pid = libc::pid_t::try_from(child.id())?;
    // The output is read while waiting: a child whose output fills the pipe
    // would otherwise block, and look like it hangs
    let (sender, receiver) = std::sync::mpsc::channel();
    std::thread::spawn(move || sender.send(child.wait_with_output()));
    match receiver.recv_timeout(TIMEOUT) {
        Ok(output) => output.wrap_err("Failed to collect rastair output"),
        Err(_) => {
            // SAFETY: kill(2) has no memory-safety preconditions, and the pid
            // is still the child's: the thread above has not reaped it.
            unsafe { libc::kill(pid, libc::SIGKILL) };
            bail!("rastair did not finish within {TIMEOUT:?}, it probably hangs");
        }
    }
}

pub fn partial(path: &Path) -> PathBuf {
    let mut partial = path.as_os_str().to_owned();
    partial.push(".partial");
    PathBuf::from(partial)
}

/// A failed run leaves its output only at `<path>.partial`.
pub fn only_partial_output(path: &Path) -> Result<PathBuf> {
    ensure!(!path.exists(), "{} should not exist after a failed run", path.display());
    let partial = partial(path);
    ensure!(partial.exists(), "{} should be left behind", partial.display());
    Ok(partial)
}

/// Data lines of a plain or BGZF-compressed text file.
pub fn records(path: &Path) -> Result<Vec<String>> {
    let text = if path.to_string_lossy().contains(".gz") {
        let bytes = std::fs::read(path)?;
        ensure!(bytes.ends_with(&BGZF_EOF), "{} should be closed properly", path.display());
        let mut text = String::new();
        rust_htslib::bgzf::Reader::from_path(path)
            .wrap_err("Failed to open BGZF file")?
            .read_to_string(&mut text)
            .wrap_err("Failed to read BGZF file")?;
        text
    } else {
        std::fs::read_to_string(path).wrap_err("Failed to read file")?
    };
    Ok(vcf_content_lines(&text).map(str::to_string).collect())
}

pub fn bam_records(path: &Path) -> Result<Vec<String>> {
    use rust_htslib::bam::Read as _;
    let bytes = std::fs::read(path)?;
    ensure!(bytes.ends_with(&BGZF_EOF), "{} should be closed properly", path.display());
    let mut reader = rust_htslib::bam::Reader::from_path(path).wrap_err("Failed to open BAM")?;
    reader
        .records()
        .map(|record| {
            let record = record.wrap_err("Failed to read BAM record")?;
            Ok(format!(
                "{} {} {}",
                String::from_utf8_lossy(record.qname()),
                record.pos(),
                record.flags()
            ))
        })
        .collect()
}

/// The records of a stopped run are the first ones of a complete run.
pub fn ensure_truncated_prefix(full: &[String], stopped: &[String]) -> Result<()> {
    ensure!(!stopped.is_empty(), "segments before the stop should be written");
    ensure!(stopped.len() < full.len(), "output should be truncated");
    ensure!(full.starts_with(stopped), "output should be a prefix of the full output");
    Ok(())
}

/// The records of a run that lost one segment are the complete run's minus a
/// contiguous stretch in the middle.
pub fn ensure_missing_one_stretch(full: &[String], lossy: &[String]) -> Result<()> {
    ensure!(lossy.len() < full.len(), "the failed segment's records should be missing");
    let head = full.iter().zip(lossy).take_while(|(a, b)| a == b).count();
    let tail = full.iter().rev().zip(lossy.iter().rev()).take_while(|(a, b)| a == b).count();
    ensure!(head > 0 && tail > 0, "segments before and after the failed one should be written");
    ensure!(head + tail >= lossy.len(), "only one stretch should be missing");
    Ok(())
}

pub fn ensure_killed_by(output: &Output, signal: i32) -> Result<()> {
    ensure!(
        output.status.signal() == Some(signal),
        "should be terminated by signal {signal}, got {:?}: {}",
        output.status,
        output.stderr()
    );
    Ok(())
}

pub fn ensure_reported_once(stderr: &str) -> Result<()> {
    ensure!(!stderr.contains("The application panicked"), "reported once, as an error: {stderr}");
    Ok(())
}

/// Calls for `rastair bam`, which needs them compressed and indexed
pub fn bam_calls(dir: &Path) -> Result<PathBuf> {
    let calls = dir.join("calls.bed.gz");
    let mut cmd = rastair();
    cmd.args(CALL_TEST_BAM).args([REGION, NO_ML, "--cpgs-only", "--bed"]).arg(&calls);
    ensure_success(&run_with_timeout(with_faults(&mut cmd, &[]))?)?;
    ensure!(dir.join("calls.bed.gz.tbi").exists(), "calls should be indexed");
    Ok(calls)
}
