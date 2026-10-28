//! Generated command lines: hegel draws flag combinations for `call`, and for
//! `convert`, `per-read` and `bam` on a call's output, over a small scenario.
//! Every run either succeeds with valid, sorted, header-consistent output, or
//! fails with a clean error (exit 1, or 2 for a usage error): never a panic,
//! a signal or a hang.
//!
//! Each case runs Rastair, so the case counts are pinned; set
//! `RASTAIR_PINNED_CASES` for a deeper run.
#![cfg(unix)]

#[path = "../src/pinned.rs"]
mod pinned;
mod utils;

use hegel::{TestCase, generators as gs};
use std::{
    collections::HashMap,
    path::PathBuf,
    process::Output,
    sync::atomic::{AtomicUsize, Ordering},
    time::Duration,
};
use utils::{
    faults::{Bounded, run_bounded},
    planted::{CONTIG_LEN, Carrier, Event, FIRST_EVENT, IndelFlanks, Plan},
    scenario::{CONTIG, Hap, READ_LEN, Strand},
    *,
};

/// A CpG past the patch, by its C, 0-based.
const LONE_CPG: usize = CONTIG_LEN - 400;

/// One event of each kind, methylated unevenly, and a CpG without a BED beta,
/// written once.
struct Prepared {
    dir: TempDir,
    fasta: PathBuf,
    bam: PathBuf,
    /// A default `call` of the scenario, for `convert`.
    calls: PathBuf,
    /// The same calls as an indexed BED, for `per-read` and `bam`.
    bed_calls: PathBuf,
}

impl Prepared {
    fn new() -> Result<Self> {
        let plan = Plan {
            seed: 3,
            events: vec![
                (FIRST_EVENT, Event::Snv { alt: 1, carrier: Carrier::One }),
                (FIRST_EVENT + 20, Event::Cpg),
                (FIRST_EVENT + 40, Event::DenovoC { carrier: Carrier::Two }),
                (
                    FIRST_EVENT + 60,
                    Event::Insertion { bases: b"TT".to_vec(), carrier: Carrier::One },
                ),
                (FIRST_EVENT + 80, Event::Deletion { len: 2, carrier: Carrier::Both }),
                (FIRST_EVENT + 100, Event::Snv { alt: 0, carrier: Carrier::Both }),
                (FIRST_EVENT + 120, Event::DenovoG { carrier: Carrier::Both }),
            ],
            step: 25,
            copies: 1,
            methylation: [0.25, 0.75],
            indel_flanks: IndelFlanks::Random,
        };
        let dir = TempDir::new()?;
        // Past the patch, a CpG whose only reads are top-strand ones ending
        // on its C: none pairs the C with its G, so its BED beta is empty.
        let scenario = plan.scenario().reference_at(LONE_CPG, b"CG").pairs(
            Hap::One,
            Strand::Top,
            LONE_CPG + 1 - READ_LEN,
            3,
        );
        let fasta = scenario.write_fasta(dir.path())?;
        let bam = dir.path().join("reads.bam");
        scenario.write_bam(&bam)?;
        let calls = dir.path().join("calls.vcf.gz");
        let bed_calls = dir.path().join("calls.bed.gz");
        rastair()
            .arg("call")
            .arg("--fasta-file")
            .arg(&fasta)
            .arg(&bam)
            .arg(scenario.region())
            .arg("--vcf")
            .arg(&calls)
            .arg("--bed")
            .arg(&bed_calls)
            .silent()
            .succeeds()?;
        let lone = read_text(&bed_calls)?
            .lines()
            .find(|line| line.split('\t').nth(1) == Some(LONE_CPG.to_string().as_str()))
            .map(str::to_owned);
        ensure!(
            lone.as_deref().and_then(|line| line.split('\t').nth(4)) == Some(""),
            "the lone CpG has a beta: {lone:?}"
        );
        Ok(Self { dir, fasta, bam, calls, bed_calls })
    }
}

/// A `-l` list of up to four pieces in any order, overlapping or not, some
/// reaching past the contig's end; now and then one is backwards or on a
/// contig the BAM does not have.
#[hegel::composite]
fn regions(tc: &TestCase) -> String {
    let count = tc.draw(gs::integers::<usize>().min_value(1).max_value(4));
    let pieces: Vec<String> = (0..count)
        .map(|_| {
            let contigs = [CONTIG, CONTIG, CONTIG, CONTIG, CONTIG, CONTIG, CONTIG, "elsewhere"];
            let contig = tc.draw(gs::sampled_from(contigs.as_slice()));
            let position = || gs::integers::<usize>().max_value(CONTIG_LEN + 100);
            match tc.draw(gs::integers::<u8>().max_value(7)) {
                0 => contig.to_owned(),
                1 => format!("{contig}:{}", tc.draw(position())),
                2 => format!("{contig}:{}-{}", tc.draw(position()), tc.draw(position())),
                _ => {
                    let (a, b) = (tc.draw(position()), tc.draw(position()));
                    format!("{contig}:{}-{}", a.min(b), a.max(b))
                }
            }
        })
        .collect();
    pieces.join(" ")
}

/// A third of the time, so that a command line holds a few flags at once and
/// shrinks to none.
fn included(tc: &TestCase) -> bool {
    tc.draw(gs::integers::<u8>().max_value(2)) == 2
}

/// A `--segment-max-length`: mostly one the commands accept, the shortest
/// first, and a fifth of the time one they refuse.
fn segment_length(tc: &TestCase) -> u64 {
    if tc.draw(gs::integers::<u8>().max_value(4)) == 4 {
        tc.draw(gs::sampled_from([0, 1, 999, 1_000_001].as_slice()))
    } else {
        tc.draw(gs::integers::<u64>().min_value(1_000).max_value(3_000))
    }
}

/// An out-of-range `--segment-max-length` is a usage error, whatever else
/// the command line holds.
fn ensure_segment_length_refused(args: &[String], output: &Output) -> Result<()> {
    let refused =
        args.iter().filter_map(|arg| arg.strip_prefix("--segment-max-length=")).any(|length| {
            length.parse().is_ok_and(|length: u64| !(1_000..=1_000_000).contains(&length))
        });
    if refused {
        ensure!(
            output.status.code() == Some(2),
            "out-of-range segment length ran: {}",
            output.status
        );
    }
    Ok(())
}

/// Where `call` writes, and what the command line says about it.
#[derive(Debug, Clone)]
struct Outputs {
    vcf: Option<&'static str>,
    bed: Option<&'static str>,
}

#[hegel::composite]
fn call_command_lines(tc: &TestCase) -> (Vec<String>, Outputs) {
    let mut args: Vec<String> = Vec::new();
    let mut flag = |tc: &TestCase, arg: String| {
        if included(tc) {
            args.push(arg);
        }
    };
    let int = |tc: &TestCase, max: u64| tc.draw(gs::integers::<u64>().max_value(max));
    // Each list ends with one value or so the command rejects.
    flag(tc, format!("--region={}", tc.draw(regions())));
    flag(tc, format!("--segment-max-length={}", segment_length(tc)));
    flag(tc, format!("--segment-overlap={}", int(tc, 600)));
    flag(tc, format!("--segment-max-bytes={}", int(tc, 100_000)));
    flag(tc, "--guess-read-orientation".to_owned());
    let models = ["novaseq6000", "miseq", "hiseqxten", "0.005", "0.02", "0", "1.5", "nonsense"];
    flag(tc, format!("--error-model={}", tc.draw(gs::sampled_from(models.as_slice()))));
    flag(tc, format!("--vcf-threads={}", 1 + int(tc, 3)));
    flag(tc, format!("--threads={}", 1 + int(tc, 3)));
    flag(tc, "--unpaired".to_owned());
    flag(tc, "--keep-overlapping-reads".to_owned());
    let coverages = [0, 1, 2, 5, 1_000, 70_000];
    flag(tc, format!("--max-coverage={}", tc.draw(gs::sampled_from(coverages.as_slice()))));
    flag(tc, format!("--min-mapq={}", int(tc, 255)));
    flag(tc, format!("--min-baseq={}", int(tc, 255)));
    for strand in ["--nOT", "--nOB"] {
        let values = tc.draw(gs::sampled_from([4, 4, 4, 4, 4, 3].as_slice()));
        let values: Vec<String> = (0..values).map(|_| int(tc, 120).to_string()).collect();
        flag(tc, format!("{strand}={}", values.join(",")));
    }
    flag(tc, format!("--include-flags={}", int(tc, 4_095)));
    flag(tc, format!("--exclude-flags={}", int(tc, 4_095)));
    flag(tc, "--rescue-soft-clip-cpg".to_owned());
    let m_coverages = [1, 5, 1_000, 0];
    flag(tc, format!("--m-max-coverage={}", tc.draw(gs::sampled_from(m_coverages.as_slice()))));
    let ml = ["--ml", "--ml=0.4", "--ml=0.6", "--ml=0.9", "--ml=1", "--ml=0", "--ml=1.5"];
    flag(tc, tc.draw(gs::sampled_from(ml.as_slice())).to_owned());
    let engine = ["--snv-engine", "--snv-engine=true", "--snv-engine=false"];
    flag(tc, tc.draw(gs::sampled_from(engine.as_slice())).to_owned());
    flag(tc, "--cpgs-only".to_owned());
    flag(tc, format!("--min-indel-ao={}", int(tc, 5)));
    flag(tc, format!("--phase-min-switch-quality={}", int(tc, 99)));
    flag(tc, "--all".to_owned());
    let tags = ["NM=0", "XX=a", "RG=x", "MQ=60", "NM=x", "NM:i:0"];
    flag(tc, format!("--require-tags={}", tc.draw(gs::sampled_from(tags.as_slice()))));
    let info = ["AD", "DP", "CPG", "M5mC_Strands", "AF,NOI", "CPGnovo,MQ", "nonsense"];
    flag(tc, format!("--vcf-info-fields={}", tc.draw(gs::sampled_from(info.as_slice()))));
    let format = ["GT", "PL,DP", "M5mC,DPM5mC,ADM5mC", "ML", "GQ,PL", "PS"];
    flag(tc, format!("--vcf-format-fields={}", tc.draw(gs::sampled_from(format.as_slice()))));
    flag(tc, "--vcf-all-fields".to_owned());
    let vcf_names = ["out.vcf", "out.vcf.gz", "out.bcf", "-"];
    let vcf = tc.draw(gs::booleans()).then(|| tc.draw(gs::sampled_from(vcf_names.as_slice())));
    let bed_names = ["out.bed", "out.bed.gz"];
    let bed = tc.draw(gs::booleans()).then(|| tc.draw(gs::sampled_from(bed_names.as_slice())));
    let bed_formats = ["bed", "bed-gz"];
    if bed.is_some() {
        flag(tc, format!("--bed-format={}", tc.draw(gs::sampled_from(bed_formats.as_slice()))));
    }
    (args, Outputs { vcf, bed })
}

/// Runs that exited 0, so that a property whose every run fails at the
/// argument parser is noticed.
static SUCCEEDED: AtomicUsize = AtomicUsize::new(0);

/// Exit 0, or a clean error: 1 from Rastair or 2 from the argument parser,
/// with no panic reported.
fn ensure_clean(output: &Output) -> Result<bool> {
    let stderr = output.stderr();
    ensure!(!stderr.contains("panicked"), "a panic: {stderr}");
    match output.status.code() {
        Some(0) => {
            SUCCEEDED.fetch_add(1, Ordering::Relaxed);
            Ok(true)
        }
        Some(1 | 2) => Ok(false),
        _ => bail!("no clean exit: {}: {stderr}", output.status),
    }
}

/// Far above what the scenario needs (about 150 MB in a debug build).
const MAX_RSS_BYTES: u64 = 2 << 30;

/// Runs a generated command, which may ask for anything, bounded in time and
/// memory: a runaway fails the case before it takes the machine down.
fn run(cmd: &mut Command) -> Result<Output> {
    match run_bounded(cmd, Duration::from_secs(60), MAX_RSS_BYTES)? {
        Bounded::Finished(output) => Ok(output),
        runaway => bail!("{runaway:?}"),
    }
}

/// The regions a `-l` list asks for on the scenario's contig, 1-based and
/// inclusive; `None` when a piece is not one this test can read.
fn requested(regions: &str) -> Option<Vec<(usize, usize)>> {
    regions
        .split_whitespace()
        .map(|piece| match piece.split_once(':') {
            None => (piece == CONTIG).then_some((1, CONTIG_LEN)),
            Some((contig, span)) if contig == CONTIG => {
                let (start, end) = match span.split_once('-') {
                    Some((start, end)) => (start.parse().ok()?, end.parse().ok()?),
                    None => (span.parse().ok()?, CONTIG_LEN),
                };
                Some((start, end))
            }
            Some(_) => None,
        })
        .collect()
}

/// The header IDs a VCF declares, by kind (`INFO`, `FORMAT`, `FILTER`,
/// `contig`), in order.
fn declared(header: &str) -> HashMap<String, Vec<String>> {
    let mut ids: HashMap<String, Vec<String>> = HashMap::new();
    for line in header.lines().filter_map(|line| line.strip_prefix("##")) {
        let Some((kind, rest)) = line.split_once("=<ID=") else { continue };
        let id = rest.split([',', '>']).next().unwrap_or_default();
        ids.entry(kind.to_owned()).or_default().push(id.to_owned());
    }
    ids
}

/// VCF text: every record names declared IDs, and the records are sorted by
/// contig in header order, then position, within `regions` when given.
fn check_vcf_text(text: &str, regions: Option<&[(usize, usize)]>) -> Result<()> {
    ensure!(text.starts_with("##fileformat=VCF"), "not a VCF: {text:.200}");
    let ids = declared(text);
    let list = |kind: &str| ids.get(kind).map(Vec::as_slice).unwrap_or_default();
    let mut last: Option<(usize, usize)> = None;
    for line in vcf_content_lines(text) {
        let columns: Vec<&str> = line.split('\t').collect();
        ensure!(columns.len() == 10, "not ten columns: {line}");
        let column = |index: usize| columns.get(index).copied().unwrap_or_default();
        let contig = list("contig")
            .iter()
            .position(|id| id == column(0))
            .ok_or_else(|| eyre!("undeclared contig: {line}"))?;
        let pos: usize = column(1).parse()?;
        ensure!(last <= Some((contig, pos)), "out of order: {line}");
        last = Some((contig, pos));
        if let Some(regions) = regions {
            ensure!(
                regions.iter().any(|&(start, end)| (start..=end).contains(&pos)),
                "outside {regions:?}: {line}"
            );
        }
        for filter in column(6).split(';').filter(|filter| !matches!(*filter, "PASS" | ".")) {
            ensure!(list("FILTER").iter().any(|id| id == filter), "undeclared FILTER: {line}");
        }
        for key in column(7).split(';').filter(|field| *field != ".") {
            let key = key.split('=').next().unwrap_or_default();
            ensure!(list("INFO").iter().any(|id| id == key), "undeclared INFO {key}: {line}");
        }
        let keys: Vec<&str> = column(8).split(':').collect();
        for key in &keys {
            ensure!(list("FORMAT").iter().any(|id| id == key), "undeclared FORMAT {key}: {line}");
        }
        ensure!(column(9).split(':').count() <= keys.len(), "more values than keys: {line}");
    }
    Ok(())
}

/// A BCF reads back, and its records are sorted.
fn check_bcf(path: &Path, regions: Option<&[(usize, usize)]>) -> Result<()> {
    use rust_htslib::bcf::{Read as _, Reader};
    let mut reader = Reader::from_path(path).wrap_err("open BCF")?;
    let mut last = None;
    for record in reader.records() {
        let record = record.wrap_err("read BCF record")?;
        let here = (record.rid(), record.pos());
        ensure!(last <= Some(here), "out of order at {here:?}");
        last = Some(here);
        let pos = usize::try_from(record.pos())? + 1;
        if let Some(regions) = regions {
            ensure!(
                regions.iter().any(|&(start, end)| (start..=end).contains(&pos)),
                "{pos} is outside {regions:?}"
            );
        }
    }
    Ok(())
}

/// BED rows: tab-separated, `start < end`, sorted by contig then start.
fn check_bed_text(text: &str) -> Result<()> {
    let mut last: Option<(String, u64)> = None;
    for line in text.lines().filter(|line| !line.starts_with('#') && !line.is_empty()) {
        let columns: Vec<&str> = line.split('\t').collect();
        let (Some(contig), Some(start), Some(end)) =
            (columns.first(), columns.get(1), columns.get(2))
        else {
            bail!("short BED row: {line}");
        };
        let (start, end): (u64, u64) = (start.parse()?, end.parse()?);
        ensure!(start < end, "empty interval: {line}");
        let here = ((*contig).to_owned(), start);
        if let Some(last) = &last {
            ensure!(last.0 != here.0 || last.1 <= here.1, "out of order: {line}");
        }
        last = Some(here);
    }
    Ok(())
}

/// Plain or BGZF-compressed text, whatever the name says: a format flag
/// wins over the extension.
fn read_text(path: &Path) -> Result<String> {
    use std::io::Read as _;
    let mut text = String::new();
    if std::fs::read(path)?.starts_with(&[0x1f, 0x8b]) {
        rust_htslib::bgzf::Reader::from_path(path)
            .wrap_err("open BGZF")?
            .read_to_string(&mut text)?;
    } else {
        text = std::fs::read_to_string(path)?;
    }
    Ok(text)
}

fn check_output(path: &Path, regions: Option<&[(usize, usize)]>) -> Result<()> {
    let name = path.to_string_lossy();
    if name.ends_with(".bcf") {
        check_bcf(path, regions)
    } else if name.contains(".vcf") {
        check_vcf_text(&read_text(path)?, regions)
    } else if name.contains(".bed") {
        check_bed_text(&read_text(path)?)
    } else if name.ends_with(".bam") {
        check_bam(path)
    } else {
        bail!("unknown output {name}")
    }
}

/// A BAM reads back, sorted by position.
fn check_bam(path: &Path) -> Result<()> {
    use rust_htslib::bam::{Read as _, Reader};
    let mut reader = Reader::from_path(path).wrap_err("open BAM")?;
    let mut last = None;
    for record in reader.records() {
        let record = record.wrap_err("read BAM record")?;
        let here = (record.tid(), record.pos());
        ensure!(last <= Some(here), "out of order at {here:?}");
        last = Some(here);
    }
    Ok(())
}

fn check_call(prepared: &Prepared, args: &[String], outputs: &Outputs) -> Result<()> {
    let dir = TempDir::new_in(prepared.dir.path())?;
    let mut cmd = rastair();
    cmd.arg("call").arg("--fasta-file").arg(&prepared.fasta).arg(&prepared.bam).args(args);
    let path = |name: &str| if name == "-" { PathBuf::from(name) } else { dir.path().join(name) };
    if let Some(vcf) = outputs.vcf {
        cmd.arg("--vcf").arg(path(vcf));
    }
    if let Some(bed) = outputs.bed {
        cmd.arg("--bed").arg(path(bed));
    }
    let output = run(&mut cmd)?;
    ensure_segment_length_refused(args, &output)?;
    if !ensure_clean(&output)? {
        return Ok(());
    }
    let regions = args
        .iter()
        .find_map(|arg| arg.strip_prefix("--region="))
        .map_or(Some(vec![(1, CONTIG_LEN)]), requested);
    let regions = regions.as_deref();
    // With no file named, the output goes to stdout: VCF, or BED under `-c`.
    let stdout = output.stdout();
    if outputs.vcf == Some("-") || stdout.starts_with("##fileformat") {
        check_vcf_text(&stdout, regions)?;
    } else if outputs.vcf.is_none() && outputs.bed.is_none() {
        check_bed_text(&stdout)?;
    }
    for entry in std::fs::read_dir(dir.path())? {
        let path = entry?.path();
        if path.extension().is_some_and(|extension| extension == "csi" || extension == "tbi") {
            continue;
        }
        check_output(&path, regions).wrap_err_with(|| path.display().to_string())?;
    }
    Ok(())
}

#[hegel::standalone_function(test_cases = pinned::cases(16))]
fn call_runs_cleanly(tc: TestCase, prepared: &Prepared) {
    let (args, outputs) = tc.draw(call_command_lines());
    let result = check_call(prepared, &args, &outputs);
    assert!(result.is_ok(), "{args:?} {outputs:?}: {result:?}");
}

/// Any `call` command line runs to valid output or a clean error.
#[test]
fn generated_call_command_lines_run_cleanly() -> Result<()> {
    call_runs_cleanly(&Prepared::new()?);
    let succeeded = SUCCEEDED.load(Ordering::Relaxed);
    eprintln!("{succeeded} runs succeeded");
    ensure!(succeeded > 0, "every run failed");
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Downstream {
    Convert,
    PerRead,
    BamStandard,
    BamLegacy,
}

impl Downstream {
    /// The command reading the prepared calls and writing `out`.
    fn command(self, prepared: &Prepared, out: &Path) -> Command {
        let mut cmd = rastair();
        match self {
            Self::Convert => {
                cmd.args(["convert", "--input"]).arg(&prepared.calls).arg("--output").arg(out);
            }
            Self::PerRead => {
                cmd.args(["per-read", "--fasta-file"])
                    .arg(&prepared.fasta)
                    .arg(&prepared.bam)
                    .arg("--calls")
                    .arg(&prepared.bed_calls)
                    .arg("--bed")
                    .arg(out);
            }
            Self::BamStandard | Self::BamLegacy => {
                let mode = if self == Self::BamStandard { "standard" } else { "legacy" };
                cmd.args(["bam", mode, "--fasta-file"])
                    .arg(&prepared.fasta)
                    .arg(&prepared.bam)
                    .arg(&prepared.bed_calls)
                    .arg("--output")
                    .arg(out);
            }
        }
        cmd
    }
}

#[hegel::composite]
fn downstream_command_lines(tc: &TestCase) -> (Downstream, Vec<String>, &'static str) {
    let commands =
        [Downstream::Convert, Downstream::PerRead, Downstream::BamStandard, Downstream::BamLegacy];
    let command = tc.draw(gs::sampled_from(commands.as_slice()));
    let mut args: Vec<String> = Vec::new();
    let mut flag = |tc: &TestCase, arg: String| {
        if included(tc) {
            args.push(arg);
        }
    };
    let int = |tc: &TestCase, max: u64| tc.draw(gs::integers::<u64>().max_value(max));
    flag(tc, format!("--threads={}", 1 + int(tc, 3)));
    let output = match command {
        Downstream::Convert => {
            let formats = ["vcf", "bcf", "vcf-compressed", "bed", "bed-gz"];
            flag(tc, format!("--output-format={}", tc.draw(gs::sampled_from(formats.as_slice()))));
            flag(tc, "--bed-include-empty".to_owned());
            let names = ["out.vcf", "out.vcf.gz", "out.bcf", "out.bed", "out.bed.gz"];
            tc.draw(gs::sampled_from(names.as_slice()))
        }
        Downstream::PerRead => {
            flag(tc, format!("--region={}", tc.draw(regions())));
            flag(tc, format!("--segment-max-length={}", segment_length(tc)));
            flag(tc, format!("--segment-overlap={}", int(tc, 600)));
            flag(tc, "--guess-read-orientation".to_owned());
            flag(tc, format!("--include-flags={}", int(tc, 4_095)));
            flag(tc, format!("--exclude-flags={}", int(tc, 4_095)));
            flag(tc, format!("--max-read-length={}", int(tc, 300)));
            flag(tc, format!("--min-mapq={}", int(tc, 255)));
            flag(tc, "--exclude-ambiguous".to_owned());
            flag(tc, "--unpaired".to_owned());
            flag(tc, "--count-clipped".to_owned());
            flag(tc, "--all-reads".to_owned());
            tc.draw(gs::sampled_from(["out.bed", "out.bed.gz"].as_slice()))
        }
        Downstream::BamStandard | Downstream::BamLegacy => {
            flag(tc, format!("--region={}", tc.draw(regions())));
            flag(tc, format!("--segment-max-length={}", segment_length(tc)));
            "out.bam"
        }
    };
    (command, args, output)
}

fn check_downstream(
    prepared: &Prepared,
    command: Downstream,
    args: &[String],
    output: &str,
) -> Result<()> {
    let dir = TempDir::new_in(prepared.dir.path())?;
    let out = dir.path().join(output);
    let mut cmd = command.command(prepared, &out);
    cmd.args(args);
    let result = run(&mut cmd)?;
    ensure_segment_length_refused(args, &result)?;
    let stderr = result.stderr();
    ensure!(!stderr.contains("WARN rastair::bed::reader"), "a calls row was skipped: {stderr}");
    if !ensure_clean(&result)? {
        return Ok(());
    }
    // `convert` writes the format it is told to, whatever the name says.
    let format = args.iter().find_map(|arg| arg.strip_prefix("--output-format="));
    let checked = match format {
        Some("vcf") => dir.path().join("as.vcf"),
        Some("vcf-compressed") => dir.path().join("as.vcf.gz"),
        Some("bcf") => dir.path().join("as.bcf"),
        Some("bed") => dir.path().join("as.bed"),
        Some(_) => dir.path().join("as.bed.gz"),
        None => out.clone(),
    };
    if checked != out {
        std::fs::rename(&out, &checked)?;
    }
    check_output(&checked, None)
}

#[hegel::standalone_function(test_cases = pinned::cases(12))]
fn downstream_runs_cleanly(tc: TestCase, prepared: &Prepared) {
    let (command, args, output) = tc.draw(downstream_command_lines());
    let result = check_downstream(prepared, command, &args, output);
    assert!(result.is_ok(), "{command:?} {args:?} {output}: {result:?}");
}

/// Any `convert`, `per-read` or `bam` command line on a call's output runs
/// to valid output or a clean error.
#[test]
fn generated_downstream_command_lines_run_cleanly() -> Result<()> {
    downstream_runs_cleanly(&Prepared::new()?);
    let succeeded = SUCCEEDED.load(Ordering::Relaxed);
    eprintln!("{succeeded} runs succeeded");
    ensure!(succeeded > 0, "every run failed");
    Ok(())
}

/// `--segment-max-length=0` is a usage error in every command that takes it,
/// refused before any work. Each run is cut off at 1 GiB of memory or 15 s.
#[test]
fn a_zero_segment_length_is_no_runaway() -> Result<()> {
    let prepared = Prepared::new()?;
    let out = |name: &str| prepared.dir.path().join(name);
    let mut call = rastair();
    call.arg("call").arg("--fasta-file").arg(&prepared.fasta).arg(&prepared.bam);
    call.arg("--vcf").arg(out("zero.vcf"));
    let commands = [
        ("call", call),
        ("per-read", Downstream::PerRead.command(&prepared, &out("zero.bed"))),
        ("bam standard", Downstream::BamStandard.command(&prepared, &out("standard.bam"))),
        ("bam legacy", Downstream::BamLegacy.command(&prepared, &out("legacy.bam"))),
    ];
    for (name, mut cmd) in commands {
        cmd.arg("--segment-max-length=0");
        match run_bounded(&mut cmd, Duration::from_secs(15), 1 << 30)? {
            Bounded::Finished(output) => {
                ensure_clean(&output).wrap_err(name)?;
                ensure!(output.status.code() == Some(2), "{name}: {}", output.status);
                ensure!(
                    output.stderr().contains("--segment-max-length"),
                    "{name}: {}",
                    output.stderr()
                );
            }
            runaway => bail!("{name}: {runaway:?}"),
        }
    }
    Ok(())
}
