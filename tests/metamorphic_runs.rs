//! Metamorphic whole runs: a transformation of the input or the command line
//! that must not change the answer, applied to a generated scenario, and the
//! two outputs compared.
//!
//! Every run writes every field (`--vcf-all-fields`), so any difference shows.
//! Each case runs Rastair at least twice, so the case counts are pinned; set
//! `RASTAIR_PINNED_CASES` for a deeper run.

#[path = "../src/pinned.rs"]
mod pinned;
mod utils;

use hegel::{TestCase, generators as gs};
use std::path::PathBuf;
use utils::{
    planted::{
        self, CONTIG_LEN, Carrier, Event, FIRST_EVENT, IndelFlanks, Kinds, ORIGIN, PATCH, Plan,
    },
    scenario::{self, CONTIG, Hap, Scenario},
    *,
};

/// A scenario written to disk once, called as often as a property needs.
struct Prepared {
    dir: TempDir,
    fasta: PathBuf,
    bam: PathBuf,
}

impl Prepared {
    fn new(scenario: &Scenario) -> Result<Self> {
        let dir = TempDir::new()?;
        let fasta = scenario.write_fasta(dir.path())?;
        let bam = dir.path().join("reads.bam");
        scenario.write_bam(&bam)?;
        Ok(Self { dir, fasta, bam })
    }

    /// The records of a `call` writing `--vcf <out>` (any of `.vcf`,
    /// `.vcf.gz` and `.bcf`), as VCF lines.
    fn call(&self, out: &str, args: &[&str]) -> Result<Vec<String>> {
        let out = self.dir.path().join(out);
        let output = rastair()
            .arg("call")
            .arg("--fasta-file")
            .arg(&self.fasta)
            .arg(&self.bam)
            .args(["--vcf-all-fields", "--vcf"])
            .arg(&out)
            .args(args)
            .output()?;
        ensure!(output.status.success(), "call {args:?} failed: {}", output.stderr());
        vcf_records(&out)
    }
}

/// The data lines of a VCF, `.vcf.gz` or BCF, as htslib spells them.
fn vcf_records(path: &Path) -> Result<Vec<String>> {
    use rust_htslib::bcf::{Read as _, Reader};
    let mut reader = Reader::from_path(path).wrap_err("open VCF")?;
    reader.records().map(|record| Ok(record?.to_vcf_string()?.trim_end().to_owned())).collect()
}

fn whole() -> String {
    format!("--region={CONTIG}:1-{CONTIG_LEN}")
}

/// The whole contig, and `--all` when drawn, so that the records of
/// rejected alts and uncovered CpGs are compared as well.
fn base_args(all: bool) -> Vec<String> {
    let mut args = vec![whole()];
    if all {
        args.push("--all".to_owned());
    }
    args
}

/// Owned arguments as the `&str`s `Prepared::call` takes.
fn strs(args: &[String]) -> Vec<&str> {
    args.iter().map(String::as_str).collect()
}

fn kinds(indels: bool) -> Kinds {
    Kinds { indels, max_events: 6, indel_flanks: IndelFlanks::Random }
}

/// The worker count changes how segments are scheduled, never what is called.
#[hegel::test(test_cases = pinned::cases(4))]
fn thread_count_does_not_change_the_calls(tc: TestCase) {
    let plan = tc.draw(planted::plans(kinds(true)));
    let threads = tc.draw(gs::integers::<usize>().min_value(2).max_value(6));
    let all = tc.draw(gs::booleans());
    let check = || -> Result<()> {
        let prepared = Prepared::new(&plan.scenario())?;
        // The shortest segments, so there is more than one to schedule.
        let segments = ["--segment-max-length=1000", "--segment-overlap=200"];
        let mut args = base_args(all);
        args.extend(segments.map(str::to_owned));
        let run = |threads: usize| {
            let mut args = args.clone();
            args.push(format!("--threads={threads}"));
            prepared.call(&format!("{threads}.vcf"), &strs(&args))
        };
        let one = run(1)?;
        ensure!(run(threads)? == one, "--threads={threads} differs from one thread");
        Ok(())
    };
    let result = check();
    assert!(result.is_ok(), "{result:?}");
}

/// Segment cuts with enough overlap for every read of a column do not change
/// the calls, but for phase: a block is solved within its segment and cut to
/// the segment's core.
#[hegel::test(test_cases = pinned::cases(4))]
fn segment_cuts_do_not_change_the_calls(tc: TestCase) {
    let plan = tc.draw(planted::plans(kinds(true)));
    // A read pair spans INSERT + READ_LEN = 350 bases. The shortest legal
    // segment is 1,000 bases; up to 1,600 the first cut falls among the events.
    let length = tc.draw(gs::integers::<usize>().min_value(1_000).max_value(1_600));
    let overlap = tc.draw(gs::integers::<usize>().min_value(350).max_value(450));
    let all = tc.draw(gs::booleans());
    let check = || -> Result<()> {
        let prepared = Prepared::new(&plan.scenario())?;
        let mut args = base_args(all);
        let one: Vec<String> =
            prepared.call("one.vcf", &strs(&args))?.iter().map(|line| unphased(line)).collect();
        args.extend([
            format!("--segment-max-length={length}"),
            format!("--segment-overlap={overlap}"),
        ]);
        let cut: Vec<String> =
            prepared.call("cut.vcf", &strs(&args))?.iter().map(|line| unphased(line)).collect();
        ensure!(cut == one, "cut:\n{cut:#?}\none segment:\n{one:#?}");
        Ok(())
    };
    let result = check();
    assert!(result.is_ok(), "{result:?}");
}

/// Regions given as overlapping, repeated, adjacent pieces in any order are
/// merged into the region they cover (H1).
#[hegel::test(test_cases = pinned::cases(4))]
fn overlapping_regions_are_the_region_they_cover(tc: TestCase) {
    let plan = tc.draw(planted::plans(kinds(true)));
    let cuts = tc.draw(
        gs::vecs(gs::integers::<usize>().min_value(1).max_value(CONTIG_LEN))
            .min_size(1)
            .max_size(4),
    );
    let all = tc.draw(gs::booleans());
    // Pieces from consecutive cuts, each reaching a little past the next cut
    // or stopping right at it, then rotated so they come out of order.
    let mut bounds: Vec<usize> = cuts.clone();
    bounds.extend([1, CONTIG_LEN]);
    bounds.sort_unstable();
    bounds.dedup();
    let overlap = tc.draw(gs::integers::<usize>().max_value(300));
    let mut pieces: Vec<String> = bounds
        .windows(2)
        .filter_map(|pair| match pair {
            [start, end] => Some(format!("{CONTIG}:{start}-{}", (end + overlap).min(CONTIG_LEN))),
            _ => None,
        })
        .collect();
    let rotate = tc.draw(gs::integers::<usize>().max_value(pieces.len().saturating_sub(1)));
    pieces.rotate_left(rotate);
    if tc.draw(gs::booleans()) {
        pieces.extend(pieces.first().cloned());
    }
    let check = || -> Result<()> {
        let prepared = Prepared::new(&plan.scenario())?;
        let mut args = base_args(all);
        let whole_run = prepared.call("whole.vcf.gz", &strs(&args))?;
        args[0] = format!("--region={}", pieces.join(" "));
        ensure!(prepared.call("pieces.vcf.gz", &strs(&args))? == whole_run, "{pieces:?} differs");
        Ok(())
    };
    let result = check();
    assert!(result.is_ok(), "{result:?}");
}

/// Calling two halves of the contig on their own and concatenating the
/// records gives the records of one run over the whole contig, but for phase:
/// a half cannot phase a site with one in the other half.
#[hegel::test(test_cases = pinned::cases(4))]
fn separate_halves_concatenate_to_the_whole(tc: TestCase) {
    let plan = tc.draw(planted::plans(kinds(true)));
    let cut = tc.draw(gs::integers::<usize>().min_value(ORIGIN).max_value(ORIGIN + PATCH));
    let all = tc.draw(gs::booleans());
    let check = || -> Result<()> {
        let prepared = Prepared::new(&plan.scenario())?;
        let run = |name: &str, region: String| {
            let mut args = base_args(all);
            args[0] = region;
            prepared.call(name, &strs(&args))
        };
        let whole_run: Vec<String> =
            run("whole.vcf", whole())?.iter().map(|l| unphased(l)).collect();
        let mut halves = run("left.vcf", format!("--region={CONTIG}:1-{cut}"))?;
        halves.extend(run("right.vcf", format!("--region={CONTIG}:{}-{CONTIG_LEN}", cut + 1))?);
        let halves: Vec<String> = halves.iter().map(|l| unphased(l)).collect();
        ensure!(halves == whole_run, "cut after {cut}:\n{halves:#?}\n{whole_run:#?}");
        Ok(())
    };
    let result = check();
    assert!(result.is_ok(), "{result:?}");
}

/// The order of the reads at one position and their names carry no
/// information.
#[hegel::test(test_cases = pinned::cases(4))]
fn read_order_and_names_do_not_change_the_calls(tc: TestCase) {
    let plan = tc.draw(planted::plans(kinds(true)));
    let seed = tc.draw(gs::integers::<u64>().min_value(1));
    let all = tc.draw(gs::booleans());
    let check = || -> Result<()> {
        let scenario = plan.scenario();
        let args = base_args(all);
        let before = Prepared::new(&scenario)?.call("before.vcf", &strs(&args))?;
        let shuffled = scenario.shuffled(seed).read_names("fragment.");
        let after = Prepared::new(&shuffled)?.call("after.vcf", &strs(&args))?;
        ensure!(after == before, "shuffled with seed {seed}");
        Ok(())
    };
    let result = check();
    assert!(result.is_ok(), "{result:?}");
}

/// Reads somewhere else on the contig, with their own events, do not change
/// the records of the patch.
#[hegel::test(test_cases = pinned::cases(4))]
fn reads_elsewhere_do_not_change_the_calls(tc: TestCase) {
    let plan = tc.draw(planted::plans(kinds(true)));
    let elsewhere = tc.draw(gs::integers::<usize>().min_value(1).max_value(4));
    let all = tc.draw(gs::booleans());
    let check = || -> Result<()> {
        // A het SNV every 97 bases, 500 bases past the patch, and a CpG there
        // that only the reads of `after` cover: under `--all` `before` writes
        // it uncovered.
        let far = ORIGIN + PATCH + 500..CONTIG_LEN;
        let scenario = plan.scenario().reference_at(far.start + 50, b"CG");
        let args = base_args(all);
        let patch = |records: Vec<String>| -> Vec<String> {
            records
                .into_iter()
                .filter(|line| position(line).is_some_and(|pos| pos <= far.start))
                .collect()
        };
        let before = patch(Prepared::new(&scenario)?.call("before.vcf", &strs(&args))?);
        let mut other = scenario;
        for (index, pos) in far.clone().step_by(97).skip(4).take(elsewhere).enumerate() {
            other = other.snv(pos, if index % 2 == 0 { Hap::One } else { Hap::Two });
        }
        let other = other.tiled(far.start, far.end, 25, 1);
        let after = patch(Prepared::new(&other)?.call("after.vcf", &strs(&args))?);
        ensure!(after == before, "{after:#?}\n{before:#?}");
        Ok(())
    };
    let result = check();
    assert!(result.is_ok(), "{result:?}");
}

/// 1-based POS of a VCF line.
fn position(line: &str) -> Option<usize> {
    line.split('\t').nth(1)?.parse().ok()
}

/// The records of `.vcf`, `.vcf.gz` and `.bcf` output are the same.
#[hegel::test(test_cases = pinned::cases(4))]
fn output_formats_hold_the_same_records(tc: TestCase) {
    let plan = tc.draw(planted::plans(kinds(true)));
    let all = tc.draw(gs::booleans());
    let check = || -> Result<()> {
        let prepared = Prepared::new(&plan.scenario())?;
        let args = base_args(all);
        let vcf = prepared.call("out.vcf", &strs(&args))?;
        let text = std::fs::read_to_string(prepared.dir.path().join("out.vcf"))?;
        ensure!(vcf_content_lines(&text).collect::<Vec<_>>() == vcf, "htslib respells the VCF");
        for name in ["out.vcf.gz", "out.bcf"] {
            ensure!(prepared.call(name, &strs(&args))? == vcf, "{name} differs");
        }
        Ok(())
    };
    let result = check();
    assert!(result.is_ok(), "{result:?}");
}

/// `convert` of a call's VCF writes the BED `call --bed` writes.
#[hegel::test(test_cases = pinned::cases(4))]
fn convert_writes_the_bed_of_call(tc: TestCase) {
    let plan = tc.draw(planted::plans(kinds(true)));
    let all = tc.draw(gs::booleans());
    let check = || -> Result<()> {
        let prepared = Prepared::new(&plan.scenario())?;
        let dir = prepared.dir.path();
        let (vcf, called, converted) =
            (dir.join("out.vcf"), dir.join("called.bed"), dir.join("converted.bed"));
        rastair()
            .arg("call")
            .arg("--fasta-file")
            .arg(&prepared.fasta)
            .arg(&prepared.bam)
            .args(base_args(all))
            .arg("--vcf")
            .arg(&vcf)
            .arg("--bed")
            .arg(&called)
            .silent()
            .succeeds()?;
        rastair()
            .args(["convert", "--input"])
            .arg(&vcf)
            .arg("--output")
            .arg(&converted)
            .silent()
            .succeeds()?;
        let (called, converted) =
            (std::fs::read_to_string(called)?, std::fs::read_to_string(converted)?);
        ensure!(converted == called, "call:\n{called}\nconvert:\n{converted}");
        Ok(())
    };
    let result = check();
    assert!(result.is_ok(), "{result:?}");
}

/// The reverse complement of a scenario, every read on the other strand, is
/// called as the mirror image of the scenario: at mirrored positions, with
/// complemented alleles, the same genotypes, filters and methylation.
///
/// The forest scores and genotype likelihoods are compared by the two ignored
/// properties below, each failing on a bug of its own.
#[hegel::test(test_cases = pinned::cases(4))]
fn the_strand_mirror_is_called_mirrored(tc: TestCase) {
    let plan = tc.draw(planted::plans(kinds(false)));
    let all = tc.draw(gs::booleans());
    let result = check_mirror(&plan, all, Compared::Calls);
    assert!(result.is_ok(), "{result:?}");
}

/// The forest scores (`ML`, `QUAL`) of the mirror are the scores of the
/// scenario. Smallest failure: one reference CpG, haplotype two 25 %
/// methylated, 32x. Its `C>T` and `G>A` alts score 0.000320460 and
/// 0.000320459; mirrored back they score 0.000340878 and 0.000325447.
#[hegel::test(test_cases = pinned::cases(4))]
#[hegel::explicit_test_case(plan = Plan {
    seed: 1,
    events: vec![(FIRST_EVENT, Event::Cpg)],
    step: 25,
    copies: 1,
    methylation: [0.0, 0.25],
    indel_flanks: IndelFlanks::Random,
}, all = false)]
#[ignore = "bug: the mirror of a CpG scores differently, still after the M5 fixes"]
fn the_strand_mirror_is_scored_mirrored(tc: TestCase) {
    let plan = tc.draw(planted::plans(kinds(false)));
    let all = tc.draw(gs::booleans());
    let result = check_mirror(&plan, all, Compared::ForestScores);
    assert!(result.is_ok(), "{result:?}");
}

/// The genotype likelihoods (`PL`, `GQ`) of the mirror are those of the
/// scenario. The case: three het `T>C` de-novo CpGs and a het `G>A`, 8 bp
/// apart on haplotype one, 32x, unmethylated. An event scored beside its
/// linked neighbours has `0/0` at PL 255 and on its own at about 90, so
/// SNV clusters cut left to right (three and one) gave the mirror other PLs.
#[hegel::test(test_cases = pinned::cases(4))]
#[hegel::explicit_test_case(plan = Plan {
    seed: 1,
    events: vec![
        (FIRST_EVENT, Event::DenovoC { carrier: Carrier::One }),
        (FIRST_EVENT + 8, Event::DenovoC { carrier: Carrier::One }),
        (FIRST_EVENT + 16, Event::DenovoC { carrier: Carrier::One }),
        (FIRST_EVENT + 24, Event::Snv { alt: 0, carrier: Carrier::One }),
    ],
    step: 25,
    copies: 1,
    methylation: [0.0, 0.0],
    indel_flanks: IndelFlanks::Random,
}, all = false)]
fn the_strand_mirror_has_mirrored_likelihoods(tc: TestCase) {
    let plan = tc.draw(planted::plans(kinds(false)));
    let all = tc.draw(gs::booleans());
    let result = check_mirror(&plan, all, Compared::Likelihoods);
    assert!(result.is_ok(), "{result:?}");
}

/// What a mirror check compares beyond the position and alleles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Compared {
    /// `FILTER`, `GT` and the methylation fields.
    Calls,
    /// `QUAL` and `ML`.
    ForestScores,
    /// `GQ` and `PL`.
    Likelihoods,
}

fn check_mirror(plan: &Plan, all: bool, compared: Compared) -> Result<()> {
    let scenario = plan.scenario();
    let args = base_args(all);
    let forward = Prepared::new(&scenario)?.call("forward.vcf", &strs(&args))?;
    let mirrored = Prepared::new(&scenario.mirrored()?)?.call("mirrored.vcf", &strs(&args))?;
    let mut mirrored: Vec<MirrorView> =
        mirrored.iter().map(|line| MirrorView::mirrored(line, compared)).collect::<Result<_>>()?;
    // Stable, so records sharing a position keep their order.
    mirrored.sort_by_key(|view| view.pos);
    let forward: Vec<MirrorView> =
        forward.iter().map(|line| MirrorView::new(line, compared)).collect::<Result<_>>()?;
    ensure!(mirrored == forward, "forward:\n{forward:#?}\nmirrored back:\n{mirrored:#?}");
    Ok(())
}

/// A VCF line without phase: `PS` and `PQ` dropped, `GT` unordered.
fn unphased(line: &str) -> String {
    let mut columns: Vec<String> = line.split('\t').map(str::to_owned).collect();
    let (Some(keys), Some(values)) = (columns.get(8), columns.get(9)) else {
        return line.to_owned();
    };
    let (keys, values): (Vec<&str>, Vec<String>) = keys
        .split(':')
        .zip(values.split(':'))
        .filter(|(key, _)| !matches!(*key, "PS" | "PQ"))
        .map(|(key, value)| (key, if key == "GT" { unordered(value) } else { value.to_owned() }))
        .unzip();
    let (keys, values) = (keys.join(":"), values.join(":"));
    columns.truncate(8);
    columns.extend([keys, values]);
    columns.join("\t")
}

/// A genotype as a sorted allele list: a phase block opens with its lower
/// allele in position order, which the mirror reverses.
fn unordered(gt: &str) -> String {
    let mut alleles: Vec<&str> = gt.split(['/', '|']).collect();
    alleles.sort_unstable();
    alleles.join("/")
}

/// The strand-free columns of an SNV or CpG record that a mirror check
/// compares; the others are `None`.
#[derive(Debug, PartialEq)]
struct MirrorView {
    pos: usize,
    reference: String,
    alts: String,
    /// `FILTER`, `GT` and `M5mC`, `DPM5mC`, `ADM5mC`.
    calls: Option<[String; 5]>,
    /// `QUAL` and `ML`.
    scores: Option<[String; 2]>,
    /// `GQ` and `PL`.
    likelihoods: Option<[String; 2]>,
}

impl MirrorView {
    fn new(line: &str, compared: Compared) -> Result<Self> {
        let columns: Vec<&str> = line.split('\t').collect();
        let column = |index: usize| -> Result<String> {
            Ok((*columns.get(index).ok_or_else(|| eyre!("short line {line}"))?).to_owned())
        };
        let field = |key: &str| format_field(line, key).unwrap_or(".").to_owned();
        Ok(Self {
            pos: column(1)?.parse()?,
            reference: column(3)?,
            alts: column(4)?,
            calls: (compared == Compared::Calls)
                .then(|| {
                    Ok::<_, color_eyre::Report>([
                        column(6)?,
                        unordered(&field("GT")),
                        field("M5mC"),
                        field("DPM5mC"),
                        field("ADM5mC"),
                    ])
                })
                .transpose()?,
            scores: (compared == Compared::ForestScores)
                .then(|| Ok::<_, color_eyre::Report>([column(5)?, field("ML")]))
                .transpose()?,
            likelihoods: (compared == Compared::Likelihoods).then(|| [field("GQ"), field("PL")]),
        })
    }

    /// The record of the mirrored scenario, at the forward position.
    fn mirrored(line: &str, compared: Compared) -> Result<Self> {
        let view = Self::new(line, compared)?;
        let complement = |alleles: &str| -> String {
            alleles.bytes().map(|base| char::from(scenario::complement(base))).collect()
        };
        Ok(Self {
            pos: CONTIG_LEN + 1 - view.pos,
            reference: complement(&view.reference),
            alts: complement(&view.alts),
            ..view
        })
    }
}
