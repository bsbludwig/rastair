//! Property test: whatever fault hits a pipeline, wherever and whenever, the
//! run ends (no hang), its exit status says what happened, and what it leaves
//! behind is trustworthy:
//!
//! - a successful run leaves exactly the complete output, at its final path;
//! - a failed run never leaves a file at the final path, and its `.partial`
//!   file is a properly closed prefix of the complete output;
//! - a killed run never leaves a file at the final path either;
//! - a segment that failed on its own is the only thing missing;
//! - a crash is reported once.
//!
//! Fault injection is compiled out of release builds, so this only runs in
//! debug builds. Set `PROPTEST_CASES` to run more cases.
#![cfg(unix)]
#![cfg_attr(
    not(debug_assertions),
    allow(dead_code, unused_imports, reason = "the test is ignored in release builds")
)]

mod utils;
use proptest::{
    prelude::*,
    test_runner::{Config, TestCaseError, TestRunner},
};
use std::{collections::HashMap, path::PathBuf, process::Output};
use utils::{faults::*, *};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Pipeline {
    CallVcf,
    CallVcfGz,
    CallBedGz,
    PerReadBed,
    PerReadBedGz,
    Bam,
}

impl Pipeline {
    const ALL: [Self; 6] = [
        Self::CallVcf,
        Self::CallVcfGz,
        Self::CallBedGz,
        Self::PerReadBed,
        Self::PerReadBedGz,
        Self::Bam,
    ];

    fn file_name(self) -> &'static str {
        match self {
            Self::CallVcf => "out.vcf",
            Self::CallVcfGz => "out.vcf.gz",
            Self::CallBedGz | Self::PerReadBedGz => "out.bed.gz",
            Self::PerReadBed => "out.bed",
            Self::Bam => "out.bam",
        }
    }

    fn run(self, calls: &Path, out: &Path, threads: usize, faults: Faults<'_>) -> Result<Output> {
        let threads = format!("--threads={threads}");
        let extra = [threads.as_str()];
        match self {
            Self::CallVcf | Self::CallVcfGz => call(out, &[&extra[..], &[NO_ML]].concat(), faults),
            Self::CallBedGz => call_to("--bed", out, &[&extra[..], &[NO_ML]].concat(), faults),
            Self::PerReadBed | Self::PerReadBedGz => per_read(out, &extra, faults),
            Self::Bam => bam(calls, out, &extra, faults),
        }
    }

    /// Checks that compressed files are closed properly, too
    fn records(self, path: &Path) -> Result<Vec<String>> {
        match self {
            Self::Bam => bam_records(path),
            _ => records(path),
        }
    }

    /// `bam` fails the run on a segment error instead of skipping the segment
    fn skips_failed_segments(self) -> bool {
        self != Self::Bam
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Panic,
    Error,
    Signal(&'static str),
}

/// Both are reached once per segment
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Point {
    Worker,
    Writer,
}

impl Point {
    fn name(self) -> &'static str {
        match self {
            Self::Worker => "worker",
            Self::Writer => "writer",
        }
    }
}

#[derive(Debug, Clone)]
struct Scenario {
    pipeline: Pipeline,
    kind: Kind,
    point: Point,
    /// Also past the number of segments, where the fault never happens
    nth: usize,
    threads: usize,
}

impl Scenario {
    fn fault(&self) -> (&'static str, String) {
        let at = format!("{}@{}", self.point.name(), self.nth);
        match self.kind {
            Kind::Panic => (INJECT_PANIC, at),
            Kind::Error => (INJECT_ERROR, at),
            Kind::Signal(signal) => (INJECT_SIGNAL, format!("{signal}:{at}")),
        }
    }

    fn expected(&self) -> Expected {
        if self.nth > SEGMENTS {
            return Expected::Complete;
        }
        match (self.kind, self.point) {
            (Kind::Error, Point::Worker) if self.pipeline.skips_failed_segments() => {
                Expected::MissingSegment
            }
            (Kind::Panic | Kind::Error, _) => Expected::Prefix,
            (Kind::Signal(signal), _) => Expected::Killed(signal_number(signal)),
        }
    }
}

fn signal_number(name: &str) -> i32 {
    match name {
        "INT" => libc::SIGINT,
        _ => libc::SIGTERM,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Expected {
    /// Successful run with all output
    Complete,
    /// Failed, the output is a properly closed start of the complete one, also
    /// when the writer panicked
    Prefix,
    /// Failed, the output lacks one segment's records
    MissingSegment,
    /// Killed before it could close its outputs
    Killed(i32),
}

fn scenarios() -> impl Strategy<Value = Scenario> {
    let pipeline = proptest::sample::select(Pipeline::ALL.to_vec());
    let kind = prop_oneof![
        Just(Kind::Panic),
        Just(Kind::Error),
        Just(Kind::Signal("INT")),
        Just(Kind::Signal("TERM")),
    ];
    let point = prop_oneof![Just(Point::Worker), Just(Point::Writer)];
    (pipeline, kind, point, 1..=SEGMENTS + 4, 1_usize..=4).prop_map(
        |(pipeline, kind, point, nth, threads)| Scenario { pipeline, kind, point, nth, threads },
    )
}

/// Complete outputs of each pipeline, computed once
struct References {
    dir: TempDir,
    calls: PathBuf,
    records: HashMap<Pipeline, Vec<String>>,
}

impl References {
    fn new() -> Result<Self> {
        let dir = TempDir::new()?;
        let calls = bam_calls(dir.path())?;
        let mut records = HashMap::new();
        for pipeline in Pipeline::ALL {
            let out = dir.path().join(format!("{pipeline:?}-{}", pipeline.file_name()));
            ensure_success(&pipeline.run(&calls, &out, 2, &[])?)?;
            let reference = pipeline.records(&out)?;
            ensure!(!reference.is_empty(), "{pipeline:?} should produce records");
            records.insert(pipeline, reference);
        }
        Ok(Self { dir, calls, records })
    }
}

fn check(references: &References, scenario: &Scenario) -> Result<()> {
    let dir = TempDir::new_in(references.dir.path())?;
    let out = dir.path().join(scenario.pipeline.file_name());
    let (var, spec) = scenario.fault();
    let output =
        scenario.pipeline.run(&references.calls, &out, scenario.threads, &[(var, &spec)])?;
    let stderr = output.stderr();
    let full = references.records.get(&scenario.pipeline).ok_or_else(|| eyre!("no reference"))?;

    ensure!(stderr.matches("The application panicked").count() <= 1, "reported twice: {stderr}");
    let expected = scenario.expected();
    match expected {
        Expected::Complete => {
            ensure_success(&output)?;
            ensure!(!partial(&out).exists(), "no partial output should be left");
            ensure!(&scenario.pipeline.records(&out)? == full, "output should be complete");
        }
        Expected::Killed(signal) => {
            ensure_killed_by(&output, signal)?;
            ensure!(!out.exists(), "a killed run never moves its output into place");
        }
        Expected::Prefix => {
            ensure_failed(&output)?;
            let written = scenario.pipeline.records(&only_partial_output(&out)?)?;
            ensure!(full.starts_with(&written), "output should be a prefix of the full one");
        }
        Expected::MissingSegment => {
            ensure_failed(&output)?;
            let written = scenario.pipeline.records(&only_partial_output(&out)?)?;
            ensure_one_stretch_missing(full, &written)?;
        }
    }
    Ok(())
}

fn ensure_failed(output: &Output) -> Result<()> {
    let status = output.status;
    ensure!(status.code() == Some(1), "expected a failure, got {status:?}: {}", output.stderr());
    Ok(())
}

/// `written` is `full` without one contiguous stretch (which may be empty: the
/// failed segment may have had no records)
fn ensure_one_stretch_missing(full: &[String], written: &[String]) -> Result<()> {
    ensure!(written.len() <= full.len(), "more records than in the complete output");
    let head = full.iter().zip(written).take_while(|(a, b)| a == b).count();
    let tail = full.iter().rev().zip(written.iter().rev()).take_while(|(a, b)| a == b).count();
    ensure!(
        head + tail.min(written.len() - head) == written.len(),
        "more than one stretch of records is missing"
    );
    Ok(())
}

#[test]
#[cfg_attr(not(debug_assertions), ignore = "fault injection only exists in debug builds")]
fn pipelines_are_crash_safe() -> Result<()> {
    let references = References::new()?;
    let cases = std::env::var("PROPTEST_CASES").ok().and_then(|n| n.parse().ok()).unwrap_or(48);
    let mut runner = TestRunner::new(Config {
        cases,
        // Each case runs Rastair, so shrinking is slow; the cases are small anyway
        max_shrink_iters: 32,
        source_file: Some(file!()),
        ..Config::default()
    });
    runner
        .run(&scenarios(), |scenario| {
            check(&references, &scenario)
                .map_err(|error| TestCaseError::fail(format!("{scenario:?}: {error:?}")))
        })
        .map_err(|error| eyre!("{error}"))
}
