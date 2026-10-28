#![expect(non_snake_case, reason = "readable test names")]

mod utils;
use insta::assert_compact_debug_snapshot;
use std::collections::{BTreeMap, BTreeSet};
use utils::*;

#[derive(Debug, Clone, PartialEq)]
struct CpgBedCall {
    pos: u32,
    strand: char,
    /// Missing where no read shows the CpG.
    beta_est: Option<f64>,
    unmod_count: u32,
    mod_count: u32,
}

#[derive(Debug, PartialEq)]
struct CpgEvidenceComparisonSummary {
    flag_calls: usize,
    evidence_calls: usize,
    shared_calls: usize,
    flag_only: usize,
    evidence_only: usize,
    mod_diff_positions: usize,
    unmod_diff_positions: usize,
    mean_abs_beta_diff: f64,
    max_abs_beta_diff: f64,
    mean_abs_mod_diff: f64,
    max_abs_mod_diff: u32,
    mean_abs_unmod_diff: f64,
    max_abs_unmod_diff: u32,
    largest_beta_diffs: Vec<CpgEvidenceDiff>,
}

#[derive(Debug, PartialEq)]
struct CpgEvidenceDiff {
    pos: u32,
    strand: char,
    flag_beta: Option<f64>,
    evidence_beta: Option<f64>,
    beta_abs_diff: f64,
    flag_mod: u32,
    evidence_mod: u32,
    mod_abs_diff: u32,
    flag_unmod: u32,
    evidence_unmod: u32,
    unmod_abs_diff: u32,
}

fn parse_cpg_bed(path: &std::path::Path) -> Result<BTreeMap<(u32, char), CpgBedCall>> {
    let text = std::fs::read_to_string(path).wrap_err("read BED output")?;
    let mut calls = BTreeMap::new();

    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        ensure!(fields.len() >= 8, "BED line has too few fields: {line}");

        let pos: u32 = fields[1].parse().wrap_err("parse start")?;
        let beta_est = match fields[4] {
            "" => None,
            beta => Some(beta.parse::<f64>().wrap_err("parse beta_est")?),
        };
        let strand = fields[5].chars().next().ok_or_else(|| eyre!("missing strand"))?;
        let unmod_count: u32 = fields[6].parse().wrap_err("parse unmod")?;
        let mod_count: u32 = fields[7].parse().wrap_err("parse mod")?;

        calls.insert((pos, strand), CpgBedCall { pos, strand, beta_est, unmod_count, mod_count });
    }

    Ok(calls)
}

fn compare_cpg_calls(
    flag_calls: &BTreeMap<(u32, char), CpgBedCall>,
    evidence_calls: &BTreeMap<(u32, char), CpgBedCall>,
) -> CpgEvidenceComparisonSummary {
    let flag_keys: BTreeSet<_> = flag_calls.keys().copied().collect();
    let evidence_keys: BTreeSet<_> = evidence_calls.keys().copied().collect();

    let shared_keys: Vec<_> = flag_keys.intersection(&evidence_keys).copied().collect();
    let flag_only = flag_keys.difference(&evidence_keys).count();
    let evidence_only = evidence_keys.difference(&flag_keys).count();

    let mut diffs = Vec::with_capacity(shared_keys.len());
    let mut beta_sum = 0.0_f64;
    let mut mod_sum = 0.0_f64;
    let mut unmod_sum = 0.0_f64;
    let mut max_beta = 0.0_f64;
    let mut max_mod = 0_u32;
    let mut max_unmod = 0_u32;
    let mut mod_diff_positions = 0_usize;
    let mut unmod_diff_positions = 0_usize;

    for key in shared_keys {
        let flag = &flag_calls[&key];
        let evidence = &evidence_calls[&key];
        let beta_abs_diff = match (flag.beta_est, evidence.beta_est) {
            (Some(flag), Some(evidence)) => (flag - evidence).abs(),
            (flag, evidence) => f64::from(u8::from(flag.is_some() != evidence.is_some())),
        };
        let mod_abs_diff = flag.mod_count.abs_diff(evidence.mod_count);
        let unmod_abs_diff = flag.unmod_count.abs_diff(evidence.unmod_count);

        beta_sum += beta_abs_diff;
        mod_sum += f64::from(mod_abs_diff);
        unmod_sum += f64::from(unmod_abs_diff);
        max_beta = max_beta.max(beta_abs_diff);
        max_mod = max_mod.max(mod_abs_diff);
        max_unmod = max_unmod.max(unmod_abs_diff);
        mod_diff_positions += usize::from(mod_abs_diff > 0);
        unmod_diff_positions += usize::from(unmod_abs_diff > 0);

        diffs.push(CpgEvidenceDiff {
            pos: flag.pos,
            strand: flag.strand,
            flag_beta: flag.beta_est,
            evidence_beta: evidence.beta_est,
            beta_abs_diff,
            flag_mod: flag.mod_count,
            evidence_mod: evidence.mod_count,
            mod_abs_diff,
            flag_unmod: flag.unmod_count,
            evidence_unmod: evidence.unmod_count,
            unmod_abs_diff,
        });
    }

    diffs.sort_by(|a, b| {
        b.beta_abs_diff
            .total_cmp(&a.beta_abs_diff)
            .then_with(|| b.mod_abs_diff.cmp(&a.mod_abs_diff))
            .then_with(|| b.unmod_abs_diff.cmp(&a.unmod_abs_diff))
            .then_with(|| a.pos.cmp(&b.pos))
            .then_with(|| a.strand.cmp(&b.strand))
    });
    diffs.truncate(10);

    let shared_count = flag_keys.intersection(&evidence_keys).count();

    CpgEvidenceComparisonSummary {
        flag_calls: flag_calls.len(),
        evidence_calls: evidence_calls.len(),
        shared_calls: shared_count,
        flag_only,
        evidence_only,
        mod_diff_positions,
        unmod_diff_positions,
        mean_abs_beta_diff: if shared_count == 0 { 0.0 } else { beta_sum / shared_count as f64 },
        max_abs_beta_diff: max_beta,
        mean_abs_mod_diff: if shared_count == 0 { 0.0 } else { mod_sum / shared_count as f64 },
        max_abs_mod_diff: max_mod,
        mean_abs_unmod_diff: if shared_count == 0 { 0.0 } else { unmod_sum / shared_count as f64 },
        max_abs_unmod_diff: max_unmod,
        largest_beta_diffs: diffs,
    }
}

#[test]
fn simple_call_gives_you_vcf_on_stdout() -> Result<()> {
    apply_common_filters!();

    let call = rastair().args(CALL_TEST_BAM).args([CHR19_SMALL]).output()?;

    assert_snapshot!(call.stderr(), @r#"
    [TIME] INFO rastair::runtime::partial_output: Wrote VCF output file="-"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    let stdout = call.stdout();
    assert!(stdout.trim().starts_with("##fileformat=VCF"));
    assert_snapshot!(stdout);

    Ok(())
}

/// `--max-coverage 0` means "no limit", both to the seqair depth cap, which
/// loads every read, and to the per-column count. Read as a cap of zero there,
/// it would count none of them: an empty VCF at unbounded memory.
#[test]
fn max_coverage_zero_means_no_limit() -> Result<()> {
    let unlimited =
        rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, "--max-coverage=0"]).output()?;
    let huge =
        rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, "--max-coverage=4000000000"]).output()?;

    let body = |out: &str| {
        out.lines().filter(|l| !l.starts_with("##")).map(str::to_owned).collect::<Vec<_>>()
    };
    let unlimited_body = body(&unlimited.stdout());
    assert!(unlimited_body.len() > 1, "--max-coverage=0 produced no records at all");
    assert_eq!(unlimited_body, body(&huge.stdout()), "a cap of 0 must not cap anything");

    Ok(())
}

/// `--snv-engine=false` scores with the first forests alone, and the engine
/// still genotypes the called SNVs: the test BAM's one called SNV keeps its
/// genotype, `GQ` and `PL`, and only its forest score (`ML`, and with it `QUAL`)
/// changes.
#[test]
fn without_the_engine_forests_snvs_keep_their_engine_genotype() -> Result<()> {
    let snvs = |engine: bool| -> Result<Vec<Vec<String>>> {
        let out = rastair()
            .args(CALL_TEST_BAM)
            .args(["--region=chr19:6128400-6128450", &format!("--snv-engine={engine}")])
            .output()?;
        ensure!(out.status.success(), "{}", out.stderr());
        Ok(out
            .stdout()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').map(str::to_owned).collect::<Vec<_>>())
            .filter(|columns| columns.get(4).is_some_and(|alt| alt != "."))
            .collect())
    };
    let (off, on) = (snvs(false)?, snvs(true)?);
    let fields = |record: &Vec<String>| {
        let keys = record.get(8).map(|k| k.split(':').map(str::to_owned).collect::<Vec<_>>());
        let values = record.get(9).map(|v| v.split(':').map(str::to_owned).collect::<Vec<_>>());
        keys.zip(values).map(|(k, v)| k.into_iter().zip(v).collect::<BTreeMap<_, _>>())
    };
    let [off] = off.as_slice() else { bail!("one SNV without the engine forests: {off:?}") };
    let [on] = on.as_slice() else { bail!("one SNV with them: {on:?}") };
    let (off_fields, on_fields) = (fields(off).unwrap_or_default(), fields(on).unwrap_or_default());
    assert_eq!(off.get(1..5), on.get(1..5), "position and alleles");
    assert_eq!(off.get(6), on.get(6), "FILTER");
    assert_eq!(off_fields.get("GT").map(String::as_str), Some("1/1"));
    assert_eq!(off_fields.get("GT"), on_fields.get("GT"));
    assert_eq!(off_fields.get("GQ"), on_fields.get("GQ"));
    assert_eq!(off_fields.get("PL"), on_fields.get("PL"));
    assert_ne!(off_fields.get("ML"), on_fields.get("ML"), "the engine forests re-score it");
    Ok(())
}

#[test]
fn vcf_with_ml() -> Result<()> {
    apply_common_filters!();

    assert_cmd_snapshot!(rastair().args(CALL_TEST_BAM).arg(CHR19_SMALL).arg(
        "--ml=0.8", // explicitly set ML threshold
    ));

    Ok(())
}

#[test]
fn vcf_with_all_fields() -> Result<()> {
    apply_common_filters!();

    assert_cmd_snapshot!(rastair().args(CALL_TEST_BAM).arg(CHR19_SMALL).arg("--vcf-all-fields",));

    Ok(())
}

#[test]
fn call_accepts_guess_read_orientation_flag() -> Result<()> {
    apply_common_filters!();

    let mut call =
        rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, "--guess-read-orientation"]).output()?;

    call.succeeds()?;
    assert!(call.stdout().trim().starts_with("##fileformat=VCF"));

    Ok(())
}

#[test]
fn guess_read_orientation_stays_close_to_flag_strand_calls() -> Result<()> {
    const REGION: &str = "--region=chr19:6103000-6106000";

    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let flag_bed = temp_dir.path().join("flag-strands.bed");
    let evidence_bed = temp_dir.path().join("guess-read-orientation.bed");

    rastair()
        .args(CALL_TEST_BAM)
        .args([REGION, "--cpgs-only", "--bed"])
        .arg(&flag_bed)
        .succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([REGION, "--cpgs-only", "--guess-read-orientation", "--bed"])
        .arg(&evidence_bed)
        .succeeds()?;

    let flag_calls = parse_cpg_bed(&flag_bed)?;
    let evidence_calls = parse_cpg_bed(&evidence_bed)?;
    let summary = compare_cpg_calls(&flag_calls, &evidence_calls);

    //eprintln!("largest beta diffs: {:#?}", summary.largest_beta_diffs);

    let min_shared = ((summary.flag_calls as f64) * 0.8).ceil() as usize;
    assert!(
        summary.shared_calls >= min_shared,
        "expected at least 80% shared CpG calls, got {} shared out of {} flag-mode calls",
        summary.shared_calls,
        summary.flag_calls,
    );
    assert!(
        summary.mean_abs_beta_diff <= 0.1,
        "mean absolute beta difference too large: {}",
        summary.mean_abs_beta_diff,
    );
    assert!(
        summary.mean_abs_mod_diff <= 1.1,
        "mean absolute mod-count difference too large: {}",
        summary.mean_abs_mod_diff,
    );
    assert!(
        summary.mean_abs_unmod_diff <= 1.1,
        "mean absolute unmod-count difference too large: {}",
        summary.mean_abs_unmod_diff,
    );
    // assert!(
    //     summary.mod_diff_positions * 2 <= summary.shared_calls,
    //     "too many positions with mod-count differences: {} of {}",
    //     summary.mod_diff_positions,
    //     summary.shared_calls,
    // );
    // assert!(
    //     summary.unmod_diff_positions * 2 <= summary.shared_calls,
    //     "too many positions with unmod-count differences: {} of {}",
    //     summary.unmod_diff_positions,
    //     summary.shared_calls,
    // );
    assert_compact_debug_snapshot!(summary);

    Ok(())
}

/// Overlapping-mate dedup, end to end on a real BAM: linking happens inside
/// `Readers::pileup`, the column loop acts on it, and the counts land in the
/// BED. The synthetic differential tests in `from_seqair.rs` pin the *rule*;
/// this pins that the rule is reached at all — a mate link that silently never
/// formed (stale mate coordinates, a reader that stopped linking) would leave
/// both runs identical.
#[test]
fn overlapping_mate_dedup_lowers_coverage_end_to_end() -> Result<()> {
    const REGION: &str = "--region=chr19:6103000-6106000";

    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let deduped = temp_dir.path().join("deduped.bed");
    let kept = temp_dir.path().join("kept.bed");

    rastair()
        .args(CALL_TEST_BAM)
        .args([REGION, "--cpgs-only", "--bed"])
        .arg(&deduped)
        .succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([REGION, "--cpgs-only", "--keep-overlapping-reads", "--bed"])
        .arg(&kept)
        .succeeds()?;

    let deduped = parse_cpg_bed(&deduped)?;
    let kept = parse_cpg_bed(&kept)?;
    assert!(!deduped.is_empty(), "no CpG calls to compare");

    let total = |calls: &BTreeMap<(u32, char), CpgBedCall>| -> u32 {
        calls.values().map(|c| c.mod_count + c.unmod_count).sum()
    };
    let (deduped_total, kept_total) = (total(&deduped), total(&kept));
    assert!(
        deduped_total < kept_total,
        "dedup removed nothing: {deduped_total} observations with it, {kept_total} without"
    );

    // Every position must lose observations, never gain them: dedup only ever
    // drops one half of a pair.
    for (key, kept_call) in &kept {
        let Some(deduped_call) = deduped.get(key) else { continue };
        let before = kept_call.mod_count + kept_call.unmod_count;
        let after = deduped_call.mod_count + deduped_call.unmod_count;
        assert!(after <= before, "{key:?}: dedup increased coverage from {before} to {after}");
    }

    Ok(())
}

#[test]
fn rescue_soft_clip_cpg_adds_methylation_evidence() -> Result<()> {
    // A region wide enough to contain soft-clips that land on CpG partners.
    const REGION: &str = "--region=chr19:6100000-6120000";

    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let off_bed = temp_dir.path().join("rescue-off.bed");
    let on_bed = temp_dir.path().join("rescue-on.bed");

    rastair()
        .args(CALL_TEST_BAM)
        .args([REGION, "--cpgs-only", "--bed"])
        .arg(&off_bed)
        .succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([REGION, "--cpgs-only", "--rescue-soft-clip-cpg", "--bed"])
        .arg(&on_bed)
        .succeeds()?;

    let off = parse_cpg_bed(&off_bed)?;
    let on = parse_cpg_bed(&on_bed)?;

    // Rescue only adds evidence, so every position called without it is still
    // called with it (it may add new ones, never drop them).
    let off_keys: BTreeSet<_> = off.keys().copied().collect();
    let on_keys: BTreeSet<_> = on.keys().copied().collect();
    assert!(off_keys.is_subset(&on_keys), "rescue dropped a CpG call");

    let evidence = |calls: &BTreeMap<(u32, char), CpgBedCall>| -> u64 {
        calls.values().map(|c| u64::from(c.mod_count) + u64::from(c.unmod_count)).sum()
    };
    let changed = off.iter().filter(|(key, off_call)| on.get(key) != Some(off_call)).count();

    assert!(changed > 0, "rescue should change at least one CpG");
    assert!(
        evidence(&on) > evidence(&off),
        "rescue should add net methylation observations ({} -> {})",
        evidence(&off),
        evidence(&on),
    );

    Ok(())
}

#[test]
fn asking_for_cpgs_defaults_to_bed_output() -> Result<()> {
    apply_common_filters!();

    let call = rastair().args(CALL_TEST_BAM).args([CHR19_SMALL]).arg("-c").output()?;

    assert_snapshot!(call.stderr(), @r#"
    [TIME] INFO rastair::runtime::partial_output: Wrote BED output file="-"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    let stdout = call.stdout();
    assert!(stdout.trim().starts_with("#chr"));
    assert_snapshot!(stdout);

    Ok(())
}

#[test]
fn bed_with_ml() -> Result<()> {
    apply_common_filters!();

    assert_cmd_snapshot!(rastair().args(CALL_TEST_BAM).arg(CHR19_SMALL).arg("-c").arg(
        "--ml=0.8", // explicitly set ML threshold
    ));

    Ok(())
}

#[test]
fn writing_vcf_to_file() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.vcf");

    assert_cmd_snapshot!(
        rastair().args(CALL_TEST_BAM).args([CHR19_SMALL]).arg("--vcf").arg(&temp_file), @r#"
    success: true
    exit_code: 0
    ----- stdout -----

    ----- stderr -----
    [TIME] INFO rastair::runtime::partial_output: Wrote VCF output file=[PATH]"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    read_bcf(&temp_file).wrap_err("validate bcf file")?;

    Ok(())
}

#[test]
fn asking_for_all_variants_includes_non_passing_ones() -> Result<()> {
    apply_common_filters!();

    let mut call = rastair().args(CALL_TEST_BAM).args([CHR19_SMALL]).arg("--vcf").output()?;
    call.succeeds()?;
    let rejected = |l: &&str| l.split('\t').nth(6).is_some_and(|filter| filter != "PASS");
    assert!(!vcf_content_lines(&call.stdout()).any(|l| rejected(&l)));

    let mut call =
        rastair().args(CALL_TEST_BAM).args([CHR19_SMALL]).arg("--vcf").arg("--all").output()?;
    call.succeeds()?;
    let stdout = call.stdout();
    assert!(vcf_content_lines(&stdout).any(|l| rejected(&l)));

    Ok(())
}

#[test]
fn ask_for_cpgs_and_vcf() -> Result<()> {
    apply_common_filters!();

    let call = rastair().args(CALL_TEST_BAM).args([CHR19_SMALL]).args(["-c", "--vcf"]).output()?;

    assert_snapshot!(call.stderr(), @r#"
    [TIME] INFO rastair::runtime::partial_output: Wrote VCF output file="-"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    assert!(call.stdout().trim().starts_with("##fileformat=VCF"));

    Ok(())
}

#[test]
fn write_bcf_to_file_and_bed_to_stdout() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.bcf");

    let mut call = rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL])
        .args(["--bed=-", "--vcf"])
        .arg(&temp_file)
        .output()?;

    call.succeeds()?;

    assert!(temp_file.exists());

    read_bcf(&temp_file).wrap_err("validate bcf file")?;

    let bed = call.stdout();
    assert_snapshot!(bed);

    Ok(())
}

/// Compressed file output is coordinate-indexed automatically: a `.csi` is
/// written next to the `.bcf`, and it can be used for region queries. Regions
/// given out of tid order on the CLI must still produce a valid sorted index.
#[test]
fn compressed_output_gets_a_csi_index() -> Result<()> {
    use rust_htslib::bcf::{IndexedReader, Read};

    let temp_dir = TempDir::new()?;
    let bcf = temp_dir.path().join("out.bcf");

    rastair()
        .args(CALL_TEST_BAM)
        // Out of tid order: bacteriophage is tid 2, chr19 is tid 0.
        .args(["--region", "bacteriophage_lambda_CpG chr19", "--vcf"])
        .arg(&bcf)
        .output()?
        .succeeds()?;

    let csi = temp_dir.path().join("out.bcf.csi");
    assert!(csi.exists(), "expected a .csi index next to the .bcf");

    // Opening via IndexedReader + fetch only works with a valid index.
    let mut reader = IndexedReader::from_path(&bcf).wrap_err("open indexed bcf")?;
    let rid = reader.header().name2rid(b"chr19").wrap_err("chr19 in header")?;
    reader.fetch(rid, 0, None).wrap_err("fetch chr19 via index")?;
    let chr19_records = reader.records().count();
    assert!(chr19_records > 0, "index region query for chr19 returned no records");

    Ok(())
}

/// Overlapping, adjacent and repeated regions are merged per contig, so each
/// position is written once and in order, and the `.bcf` can be indexed.
#[test]
fn overlapping_regions_are_merged() -> Result<()> {
    const OVERLAPPING: &str =
        "chr19:6105000-6106000 chr19:6105500-6106500 chr19:6106501-6106600 chr19:6105000-6106000";

    let temp_dir = TempDir::new()?;
    let bcf = temp_dir.path().join("overlapping.bcf");
    rastair().args(CALL_TEST_BAM).args(["--region", OVERLAPPING, "--vcf"]).arg(&bcf).succeeds()?;
    assert!(temp_dir.path().join("overlapping.bcf.csi").exists());

    let lines = |region: &str, name: &str| -> Result<Vec<String>> {
        let vcf = temp_dir.path().join(name);
        rastair().args(CALL_TEST_BAM).args(["--region", region, "--vcf"]).arg(&vcf).succeeds()?;
        Ok(vcf_content_lines(&std::fs::read_to_string(&vcf)?).map(str::to_string).collect())
    };
    let overlapping = lines(OVERLAPPING, "overlapping.vcf")?;
    let merged = lines("chr19:6105000-6106600", "merged.vcf")?;
    assert!(!merged.is_empty());
    assert_eq!(overlapping, merged);

    Ok(())
}

/// Plain (uncompressed) VCF has no coordinate index, so no `.csi` is written.
#[test]
fn plain_vcf_output_has_no_index() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let vcf = temp_dir.path().join("out.vcf");

    rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, "--vcf"]).arg(&vcf).output()?.succeeds()?;

    assert!(vcf.exists());
    assert!(!temp_dir.path().join("out.vcf.csi").exists(), "plain VCF must not be indexed");

    Ok(())
}

#[test]
fn when_asked_for_bed_file_in_vcf_param_we_are_nice() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.bed");

    assert_cmd_snapshot!(
        rastair().args(CALL_TEST_BAM).args([CHR19_SMALL]).arg("--vcf").arg(&temp_file),
        @r#"
    success: true
    exit_code: 0
    ----- stdout -----

    ----- stderr -----
    [TIME] WARN rastair::call: VCF output file name ends with `.bed`/`.bed.gz`, did you mean to use `--bed` instead of `-o`/`--vcf`? Assuming you meant `--bed` and switching the output accordingly. file=[PATH]"
    [TIME] INFO rastair::runtime::partial_output: Wrote BED output file=[PATH]"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    Ok(())
}

#[test]
fn segmentation_does_not_change_bed_output() -> Result<()> {
    const REGION: &str = "--region=chr19:6000000-7000000";

    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file1 = temp_dir.path().join("test1.bed");
    let temp_file2 = temp_dir.path().join("test2.bed");

    rastair()
        .args(CALL_TEST_BAM)
        .args([REGION])
        .args(["--segment-max-length=1000", "--segment-overlap=100", "--threads=7", "--bed"])
        .arg(&temp_file1)
        .succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([REGION])
        // The SNV engine's clusters and reads reach past a one-base overlap.
        .args(["--segment-max-length=1001", "--segment-overlap=100", "--threads=4", "--bed"])
        .arg(&temp_file2)
        .succeeds()?;

    // Check that both files have the same hash
    fn hash(content: &[u8]) -> u64 {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        let mut hasher = DefaultHasher::new();
        content.hash(&mut hasher);
        hasher.finish()
    }
    let hash1 = hash(&std::fs::read(&temp_file1)?);
    let hash2 = hash(&std::fs::read(&temp_file2)?);

    assert_eq!(hash1, hash2, "Files should have identical content");

    Ok(())
}

#[test]
fn segmentation_overlaps_do_not_cause_duplicate_records() -> Result<()> {
    const REGION: &str = "--region=chr19";

    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.vcf");

    rastair()
        .args(CALL_TEST_BAM)
        .args([REGION])
        .args(["--segment-max-length=10000", "--segment-overlap=300", "--vcf"])
        .arg(&temp_file)
        .succeeds()
        .wrap_err("rastair call failed")?;

    let text = std::fs::read_to_string(&temp_file).wrap_err("read rastair 2 vcf")?;
    vcf_content_lines(&text)
        .filter_map(|line| line.split("\t").nth(1))
        .filter_map(|x| x.parse::<u32>().ok())
        .try_fold(BTreeSet::new(), |mut set, position| {
            let is_new = set.insert(position);
            if is_new { Ok(set) } else { Err(eyre!("Duplicate position found: {}", position)) }
        })?;

    Ok(())
}

#[test]
fn vcf_with_nOT_nOB() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let a = temp_dir.path().join("defaults.bcf");
    let b = temp_dir.path().join("with_args.bcf");

    rastair().args(CALL_TEST_BAM).args([CHR19_SMALL]).arg("--vcf").arg(&a).succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, "--nOT=12,12,12,12", "--nOB=12,12,12,12"])
        .arg("--vcf")
        .arg(&b)
        .succeeds()?;

    assert_compact_debug_snapshot!(get_depths(&a), @"Ok([18, 16, 16])");
    assert_compact_debug_snapshot!(get_depths(&b), @"Ok([11, 13, 13])");

    fn get_depths(path: &std::path::Path) -> Result<Vec<i32>> {
        use rust_htslib::bcf::Read;

        let mut bcf = read_bcf(path).wrap_err("invalid bcf file")?;
        let depths = bcf
            .records()
            .map(|r| {
                let field = r.unwrap().info(b"DP").integer().unwrap().unwrap();
                *field.first().unwrap()
            })
            .collect::<Vec<_>>();
        Ok(depths)
    }

    Ok(())
}

#[test]
fn vcf_field_configuration_via_cli() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let default_vcf = temp_dir.path().join("default.vcf");
    let custom_vcf = temp_dir.path().join("custom.vcf");

    // The region holds an SNV record; reference-only records lack some of the fields.
    let region = "--region=chr19:6128000-6129000";
    rastair().args(CALL_TEST_BAM).args([region, "--vcf"]).arg(&default_vcf).succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([region, "--vcf"])
        .arg(&custom_vcf)
        .args(["--vcf-info-fields=AF,MQ0,NS"])
        .succeeds()?;

    let default_content = std::fs::read_to_string(&default_vcf)?;
    let custom_content = std::fs::read_to_string(&custom_vcf)?;

    // Note: Headers are written for all fields regardless of config
    // Check the actual data lines instead
    let default_data_lines: Vec<&str> = vcf_content_lines(&default_content).collect();
    let custom_data_lines: Vec<&str> = vcf_content_lines(&custom_content).collect();

    // Default VCF should not contain AF, MQ0, NS in data lines (they are not default)
    assert!(!default_data_lines.iter().any(|l| l.contains("AF=")), "AF should not be in default");
    assert!(!default_data_lines.iter().any(|l| l.contains("MQ0=")), "MQ0 should not be in default");
    assert!(!default_data_lines.iter().any(|l| l.contains("NS=")), "NS should not be in default");

    // Default VCF should contain default fields like AD in data lines
    assert!(default_data_lines.iter().all(|l| l.contains("AD=")), "Should have AD in default");
    assert!(default_data_lines.iter().all(|l| l.contains("DP=")), "Should have DP in default");
    assert!(default_data_lines.iter().all(|l| l.contains("BQ=")), "Should have BQ in default");

    // Custom VCF should contain the additional fields in at least some data lines
    // (not all fields are present on all variant types)
    assert!(
        custom_data_lines.iter().any(|l| l.contains("AF=")),
        "Should have AF in some data lines"
    );
    assert!(
        custom_data_lines.iter().any(|l| l.contains("MQ0=")),
        "Should have MQ0 in some data lines"
    );
    assert!(
        custom_data_lines.iter().any(|l| l.contains("NS=")),
        "Should have NS in some data lines"
    );

    // Custom VCF should still have default fields
    assert!(custom_data_lines.iter().all(|l| l.contains("AD=")), "Should still have AD");

    Ok(())
}

#[test]
fn min_baseq_filter_reduces_variant_calls() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let default_vcf = temp_dir.path().join("default.vcf");
    let filtered_vcf = temp_dir.path().join("filtered.vcf");

    rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, "--vcf"]).arg(&default_vcf).succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, "--vcf"])
        .arg(&filtered_vcf)
        .arg("--min-baseq=30")
        .succeeds()?;

    let default_content = std::fs::read_to_string(&default_vcf)?;
    let filtered_content = std::fs::read_to_string(&filtered_vcf)?;

    let default_count = vcf_content_lines(&default_content).count();
    let filtered_count = vcf_content_lines(&filtered_content).count();

    assert!(
        filtered_count <= default_count,
        "Filtered VCF should have fewer or equal records than default (filtered: {}, default: {})",
        filtered_count,
        default_count
    );

    Ok(())
}

#[test]
fn min_mapq_filter_reduces_variant_calls() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let default_vcf = temp_dir.path().join("default.vcf");
    let filtered_vcf = temp_dir.path().join("filtered.vcf");

    rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, "--vcf"]).arg(&default_vcf).succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, "--vcf"])
        .arg(&filtered_vcf)
        .arg("--min-mapq=40")
        .succeeds()?;

    let default_content = std::fs::read_to_string(&default_vcf)?;
    let filtered_content = std::fs::read_to_string(&filtered_vcf)?;

    let default_count = vcf_content_lines(&default_content).count();
    let filtered_count = vcf_content_lines(&filtered_content).count();

    assert!(
        filtered_count <= default_count,
        "Filtered VCF should have fewer or equal records than default (filtered: {}, default: {})",
        filtered_count,
        default_count
    );

    Ok(())
}

#[test]
fn mbias_filter_affects_variant_calls() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let default_vcf = temp_dir.path().join("default.vcf");
    let filtered_vcf = temp_dir.path().join("filtered.vcf");

    rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, "--vcf"]).arg(&default_vcf).succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, "--vcf"])
        .arg(&filtered_vcf)
        .args(["--nOT=10,10,10,10", "--nOB=10,10,10,10"])
        .succeeds()?;

    let default_content = std::fs::read_to_string(&default_vcf)?;
    let filtered_content = std::fs::read_to_string(&filtered_vcf)?;

    let default_count = vcf_content_lines(&default_content).count();
    let filtered_count = vcf_content_lines(&filtered_content).count();

    assert!(
        filtered_count <= default_count,
        "Filtered VCF should have fewer or equal records than default (filtered: {}, default: {})",
        filtered_count,
        default_count
    );

    Ok(())
}

#[test]
fn error_model_accepts_platform_names() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let vcf_miseq = temp_dir.path().join("miseq.vcf");
    let vcf_novaseq = temp_dir.path().join("novaseq.vcf");

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, "--vcf"])
        .arg(&vcf_miseq)
        .arg("--error-model=miseq")
        .succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, "--vcf"])
        .arg(&vcf_novaseq)
        .arg("--error-model=novaseq6000")
        .succeeds()?;

    // Both should succeed - different error models can produce different results
    assert!(vcf_miseq.exists());
    assert!(vcf_novaseq.exists());

    Ok(())
}

#[test]
fn error_model_accepts_custom_error_rate() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let vcf_custom = temp_dir.path().join("custom.vcf");

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, "--vcf"])
        .arg(&vcf_custom)
        .arg("--error-model=0.005")
        .succeeds()?;

    assert!(vcf_custom.exists());

    Ok(())
}

#[test]
fn error_model_rejects_invalid_error_rate() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let vcf = temp_dir.path().join("test.vcf");

    // Error rate > 1.0 should fail
    let result = rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, "--vcf"])
        .arg(&vcf)
        .arg("--error-model=1.5")
        .output()?;

    assert!(!result.status.success());
    let stderr = result.stderr();
    assert!(
        stderr.contains("Error rate must be between 0.0 and 1.0"),
        "Expected error message about invalid range, got: {}",
        stderr
    );

    Ok(())
}

#[test]
fn read_group_single_filter_reduces_records() -> Result<()> {
    apply_common_filters!();

    const REGION: &str = "--region=chr19:6103000-6110000";

    let temp_dir = TempDir::new()?;
    let unfiltered = temp_dir.path().join("unfiltered.vcf");
    let filtered = temp_dir.path().join("filtered.vcf");

    rastair().args(CALL_TEST_BAM).args([REGION, "--vcf"]).arg(&unfiltered).succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([REGION, "--vcf"])
        .arg(&filtered)
        .arg("--require-tags=RG=mTet1-PyBr-16h-p1_S1_L001")
        .succeeds()?;

    let unfiltered_count = vcf_content_lines(&std::fs::read_to_string(&unfiltered)?).count();
    let filtered_count = vcf_content_lines(&std::fs::read_to_string(&filtered)?).count();

    assert_snapshot!(unfiltered_count, @"184");
    assert_snapshot!(filtered_count, @"179");

    Ok(())
}

/// Copy `tests/data/test.bam` into `output`, stamping each read with two extra Z tags:
/// * `XX:Z:keep` on every read,
/// * `YY:Z:lane<N>` derived from the suffix of the read's RG (e.g. `..._L001` → `lane1`),
///   `YY:Z:unknown` when the RG does not match the expected pattern.
///
/// Also builds a BAI index next to the output BAM.
fn write_bam_with_extra_tags(output: &std::path::Path) -> Result<()> {
    use rust_htslib::bam::{self, Read as _, Record, record::Aux};

    let mut reader = bam::Reader::from_path("tests/data/test.bam")?;
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(output, &header, bam::Format::Bam)?;

    let mut record = Record::new();
    while let Some(result) = reader.read(&mut record) {
        result?;
        let lane_value = match record.aux(b"RG") {
            Ok(Aux::String(rg)) => rg
                .rsplit_once('_')
                .and_then(|(_, suffix)| suffix.strip_prefix('L'))
                .and_then(|n| n.parse::<u32>().ok())
                .map(|n| format!("lane{n}"))
                .unwrap_or_else(|| "unknown".into()),
            _ => "unknown".into(),
        };
        record.push_aux(b"XX", Aux::String("keep"))?;
        record.push_aux(b"YY", Aux::String(&lane_value))?;
        writer.write(&record)?;
    }
    drop(writer);

    bam::index::build(output, None, bam::index::Type::Bai, 1)?;
    Ok(())
}

#[test]
fn multiple_different_tags_must_all_match() -> Result<()> {
    apply_common_filters!();

    const REGION: &str = "--region=chr19:6103000-6110000";

    let temp_dir = TempDir::new()?;
    let tagged_bam = temp_dir.path().join("tagged.bam");
    write_bam_with_extra_tags(&tagged_bam)?;

    let run_with = |filters: &[&str]| -> Result<usize> {
        let vcf = temp_dir.path().join(format!("out_{}.vcf", filters.join("_")));
        let mut cmd = rastair();
        cmd.args(["call", "--fasta-file=tests/data/test.fasta.gz"])
            .arg(&tagged_bam)
            .args([REGION, "--vcf"])
            .arg(&vcf);
        if !filters.is_empty() {
            cmd.arg("--require-tags");
            for f in filters {
                cmd.arg(f);
            }
        }
        cmd.succeeds()?;
        Ok(vcf_content_lines(&std::fs::read_to_string(&vcf)?).count())
    };

    let unfiltered = run_with(&[])?;
    // Tag every read carries → identical to unfiltered.
    let always_matches = run_with(&["XX=keep"])?;
    // RG-derived tag, picks one lane.
    let lane1_only = run_with(&["YY=lane1"])?;
    // Two different tags, both satisfied → equivalent to YY=lane1 alone since XX=keep is universal.
    let xx_and_lane1 = run_with(&["XX=keep", "YY=lane1"])?;
    // One filter unsatisfiable → entire set rejected even though XX=keep matches everything.
    let xx_and_missing = run_with(&["XX=keep", "YY=nonexistent"])?;

    assert_eq!(always_matches, unfiltered, "XX=keep should not exclude any reads");
    assert!(
        lane1_only < unfiltered,
        "YY=lane1 should be a strict subset ({lane1_only} vs {unfiltered})"
    );
    assert_eq!(xx_and_lane1, lane1_only, "combining with a universal tag must not change result");
    assert_eq!(xx_and_missing, 0, "one unsatisfiable filter must reject everything");

    assert_snapshot!(unfiltered, @"184");
    assert_snapshot!(lane1_only, @"179");

    Ok(())
}

#[test]
fn error_model_rejects_invalid_platform_name() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let vcf = temp_dir.path().join("test.vcf");

    let result = rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, "--vcf"])
        .arg(&vcf)
        .arg("--error-model=invalid_platform")
        .output()?;

    assert!(!result.status.success());
    let stderr = result.stderr();
    assert!(
        stderr.contains("Invalid error model") || stderr.contains("invalid value"),
        "Expected error message about invalid platform, got: {}",
        stderr
    );

    Ok(())
}

#[test]
fn cpgs_only_with_all_reports_uncovered_reference_cpgs_in_bed() -> Result<()> {
    apply_common_filters!();

    const REGION: &str = "--region=chr19:6105700-6105800";

    // The second reference CpG pair in the region (0-based 6105743..6105745,
    // i.e. 1-based 6105744) loses its coverage: either no read reaches it, or
    // every read that does fails the default `--min-mapq 1`. Both must look
    // the same in the output.
    for gap in [Gap::NoReads, Gap::ZeroMapq] {
        let temp_dir = TempDir::new()?;
        let gapped_bam = temp_dir.path().join("gapped.bam");
        write_bam_with_gap(&gapped_bam, "chr19", 6_105_743, 6_105_745, gap)?;

        let cpgs_only_bed = temp_dir.path().join("cpgs-only.bed");
        let all_bed = temp_dir.path().join("all.bed");

        rastair()
            .args(["call", "--fasta-file=tests/data/test.fasta.gz"])
            .arg(&gapped_bam)
            .args([REGION, "--cpgs-only", "--bed"])
            .arg(&cpgs_only_bed)
            .succeeds()?;

        rastair()
            .args(["call", "--fasta-file=tests/data/test.fasta.gz"])
            .arg(&gapped_bam)
            .args([REGION, "--cpgs-only", "--all", "--bed"])
            .arg(&all_bed)
            .succeeds()?;

        let cpgs_only_positions: BTreeSet<(u32, char)> =
            parse_cpg_bed(&cpgs_only_bed)?.keys().copied().collect();
        let all_positions: BTreeSet<(u32, char)> =
            parse_cpg_bed(&all_bed)?.keys().copied().collect();

        // The second CpG pair has no coverage, so plain `--cpgs-only` drops it …
        assert_eq!(
            cpgs_only_positions,
            BTreeSet::from([(6_105_711, '+'), (6_105_712, '-')]),
            "expected only the covered CpG pair under --cpgs-only ({gap:?})"
        );
        // … but `--all --cpgs-only` must report every reference CpG in the
        // region, whether it has coverage or not (both CpG pairs here).
        assert_eq!(
            all_positions,
            BTreeSet::from([
                (6_105_711, '+'),
                (6_105_712, '-'),
                (6_105_743, '+'),
                (6_105_744, '-'),
            ]),
            "expected all reference CpGs under --all --cpgs-only ({gap:?})"
        );
    }

    Ok(())
}

/// A region no read reaches still lists every reference CpG under `--all`, so
/// the outputs of two samples line up row for row.
#[test]
fn all_reports_cpgs_in_a_region_without_reads() -> Result<()> {
    use rust_htslib::faidx;

    apply_common_filters!();

    // 1-based, inclusive; `test.bam` has no read here.
    const REGION: &str = "--region=chr19:6000001-6001000";
    const FIRST: u32 = 6_000_000;
    const LAST: u32 = 6_000_999;

    // One base either side, so a CpG straddling an edge is seen whole; only
    // its half inside the region is reported.
    let fasta = faidx::Reader::from_path("tests/data/test.fasta.gz")?;
    let seq = fasta.fetch_seq_string("chr19", FIRST as usize - 1, LAST as usize + 1)?;
    let expected: BTreeSet<(u32, char)> = seq
        .as_bytes()
        .windows(2)
        .zip(FIRST - 1..)
        .filter(|(pair, _)| pair.eq_ignore_ascii_case(b"CG"))
        .flat_map(|(_, c)| [(c, '+'), (c + 1, '-')])
        .filter(|&(pos, _)| (FIRST..=LAST).contains(&pos))
        .collect();
    assert!(!expected.is_empty(), "the region has reference CpGs");

    for flags in [&["--all", "--cpgs-only"][..], &["--all"]] {
        let temp_dir = TempDir::new()?;
        let bed = temp_dir.path().join("all.bed");
        rastair()
            .args(["call", "--fasta-file=tests/data/test.fasta.gz", "tests/data/test.bam"])
            .args([REGION])
            .args(flags)
            .arg("--bed")
            .arg(&bed)
            .succeeds()?;

        let reported: BTreeSet<(u32, char)> = parse_cpg_bed(&bed)?.keys().copied().collect();
        assert_eq!(reported, expected, "every reference CpG under {flags:?}");
    }

    Ok(())
}

/// The CpG halves `--all --cpgs-only` reports for `seq` (from 0-based
/// position `offset`) over the 0-based, inclusive `first..=last`: a CpG
/// straddling an edge contributes only its half inside.
fn reference_cpgs(seq: &[u8], offset: u32, first: u32, last: u32) -> BTreeSet<(u32, char)> {
    seq.windows(2)
        .zip(offset..)
        .filter(|(pair, _)| pair.eq_ignore_ascii_case(b"CG"))
        .flat_map(|(_, c)| [(c, '+'), (c + 1, '-')])
        .filter(|&(pos, _)| (first..=last).contains(&pos))
        .collect()
}

fn reported_cpgs(
    fasta: &std::path::Path,
    bam: &std::path::Path,
    region: &str,
) -> Result<BTreeSet<(u32, char)>> {
    let temp_dir = TempDir::new()?;
    let bed = temp_dir.path().join("all.bed");
    rastair()
        .args(["call", "--fasta-file"])
        .arg(fasta)
        .arg(bam)
        .args([region, "--all", "--cpgs-only", "--bed"])
        .arg(&bed)
        .succeeds()?;
    Ok(parse_cpg_bed(&bed)?.keys().copied().collect())
}

/// A `--region` is 1-based and inclusive at both ends, and a CpG half on
/// either edge is reported: the C at the first base, the G at the last.
#[test]
fn region_edges_are_reported_under_all_cpgs_only() -> Result<()> {
    use rust_htslib::faidx;

    apply_common_filters!();

    let fasta_path = std::path::Path::new("tests/data/test.fasta.gz");
    let fasta = faidx::Reader::from_path(fasta_path)?;
    // 0-based 6105711 and 6105743 are the Cs of the region's two CpGs.
    let (offset, seq) = (6_105_700u32, fasta.fetch_seq("chr19", 6_105_700, 6_105_760)?.to_vec());
    for (region, first, last) in [
        // Starts on a CpG's C, ends on another's G: both CpGs whole.
        ("--region=chr19:6105712-6105745", 6_105_711, 6_105_744),
        // Starts on a CpG's G, ends on another's C: one half of each.
        ("--region=chr19:6105713-6105744", 6_105_712, 6_105_743),
    ] {
        let expected = reference_cpgs(&seq, offset, first, last);
        let reported =
            reported_cpgs(fasta_path, std::path::Path::new("tests/data/test.bam"), region)?;
        assert_eq!(reported, expected, "{region}");
    }

    // The same first region as BED: 0-based, half-open.
    let temp_dir = TempDir::new()?;
    let bed = temp_dir.path().join("region.bed");
    std::fs::write(&bed, "chr19\t6105711\t6105745\n")?;
    let region = format!("--region=@{}", bed.display());
    let reported = reported_cpgs(fasta_path, std::path::Path::new("tests/data/test.bam"), &region)?;
    assert_eq!(reported, reference_cpgs(&seq, offset, 6_105_711, 6_105_744), "{region}");

    Ok(())
}

/// A contig that starts and ends with a CpG and that no read reaches: `--all
/// --cpgs-only` over the whole contig lists every CpG, both ends included,
/// and a `--region` on either end reports its edge half.
#[test]
fn contig_ends_are_reported_under_all_cpgs_only() -> Result<()> {
    apply_common_filters!();

    const SEQ: &[u8] = b"CGATTACAGCGTTTACCGGATACATTGACGGCATACG";
    let last = u32::try_from(SEQ.len())? - 1;
    let temp_dir = TempDir::new()?;
    let (fasta, bam) = write_reference_without_reads(temp_dir.path(), "ctg", SEQ)?;

    let whole = reported_cpgs(&fasta, &bam, "--region=ctg")?;
    assert_eq!(whole, reference_cpgs(SEQ, 0, 0, last), "the whole contig");
    assert!(whole.contains(&(0, '+')) && whole.contains(&(last, '-')), "both ends are CpGs");

    // 1-based: the second base is the first CpG's G, the second-to-last the
    // last CpG's C.
    let inner = format!("--region=ctg:2-{}", SEQ.len() - 1);
    assert_eq!(
        reported_cpgs(&fasta, &bam, &inner)?,
        reference_cpgs(SEQ, 0, 1, last - 1),
        "{inner}"
    );

    Ok(())
}

/// A plain FASTA holding `seq` as `contig`, and an indexed BAM whose header
/// names it but which holds no read.
fn write_reference_without_reads(
    dir: &std::path::Path,
    contig: &str,
    seq: &[u8],
) -> Result<(std::path::PathBuf, std::path::PathBuf)> {
    use rust_htslib::bam::{self, header::HeaderRecord};

    let fasta = dir.join("ref.fa");
    let header_line = format!(">{contig}\n");
    let mut text = header_line.clone().into_bytes();
    text.extend_from_slice(seq);
    text.push(b'\n');
    std::fs::write(&fasta, text)?;
    std::fs::write(
        dir.join("ref.fa.fai"),
        format!(
            "{contig}\t{len}\t{offset}\t{len}\t{width}\n",
            len = seq.len(),
            offset = header_line.len(),
            width = seq.len() + 1
        ),
    )?;

    let bam_path = dir.join("empty.bam");
    let mut header = bam::Header::new();
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", contig).push_tag(b"LN", seq.len()));
    drop(bam::Writer::from_path(&bam_path, &header, bam::Format::Bam)?);
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1)?;
    Ok((fasta, bam_path))
}

#[test]
fn convert_bed_include_empty_gates_uncovered_cpgs_independently_of_call() -> Result<()> {
    apply_common_filters!();

    const REGION: &str = "--region=chr19:6105700-6105800";

    let temp_dir = TempDir::new()?;
    let gapped_bam = temp_dir.path().join("gapped.bam");
    write_bam_with_gap(&gapped_bam, "chr19", 6_105_743, 6_105_745, Gap::NoReads)?;

    // `call --all --cpgs-only` already reports the uncovered CpG pair with no
    // extra flag needed (that's the point of this refactor), so the VCF
    // contains both pairs.
    let vcf = temp_dir.path().join("all.vcf");
    rastair()
        .args(["call", "--fasta-file=tests/data/test.fasta.gz"])
        .arg(&gapped_bam)
        .args([REGION, "--cpgs-only", "--all", "-o"])
        .arg(&vcf)
        .succeeds()?;

    let without_empty = temp_dir.path().join("without-empty.bed");
    let with_empty = temp_dir.path().join("with-empty.bed");

    rastair()
        .args(["convert", "--input"])
        .arg(&vcf)
        .arg("--output")
        .arg(&without_empty)
        .succeeds()?;

    rastair()
        .args(["convert", "--input"])
        .arg(&vcf)
        .arg("--output")
        .arg(&with_empty)
        .arg("--bed-include-empty")
        .succeeds()?;

    let without_empty_positions: BTreeSet<(u32, char)> =
        parse_cpg_bed(&without_empty)?.keys().copied().collect();
    let with_empty_positions: BTreeSet<(u32, char)> =
        parse_cpg_bed(&with_empty)?.keys().copied().collect();

    // `convert` still has to decide independently whether to carry the
    // uncovered pair, already present in the VCF, into this BED output.
    assert_eq!(
        without_empty_positions,
        BTreeSet::from([(6_105_711, '+'), (6_105_712, '-')]),
        "convert without --bed-include-empty should drop the uncovered CpG pair"
    );
    assert_eq!(
        with_empty_positions,
        BTreeSet::from([(6_105_711, '+'), (6_105_712, '-'), (6_105_743, '+'), (6_105_744, '-')]),
        "convert --bed-include-empty should keep the uncovered CpG pair"
    );

    Ok(())
}

/// How [`write_bam_with_gap`] takes a window's coverage away.
#[derive(Debug, Clone, Copy)]
enum Gap {
    /// Drop every read that overlaps the window.
    NoReads,
    /// Keep the reads, with MAPQ 0.
    ZeroMapq,
}

/// Copy `tests/data/test.bam`, taking the coverage away from every read whose
/// alignment overlaps the 0-based half-open window `[win_start, win_end)` on
/// `chrom`.
fn write_bam_with_gap(
    output: &std::path::Path,
    chrom: &str,
    win_start: i64,
    win_end: i64,
    gap: Gap,
) -> Result<()> {
    use rust_htslib::bam::{self, Read as _, Record};

    let mut reader = bam::Reader::from_path("tests/data/test.bam")?;
    let header = bam::Header::from_template(reader.header());
    let tid = i32::try_from(
        bam::HeaderView::from_header(&header)
            .tid(chrom.as_bytes())
            .ok_or_else(|| eyre!("chrom not in header"))?,
    )
    .wrap_err("tid does not fit in i32")?;

    let mut writer = bam::Writer::from_path(output, &header, bam::Format::Bam)?;

    let mut record = Record::new();
    while let Some(result) = reader.read(&mut record) {
        result?;
        let overlaps = record.tid() == tid && record.pos() >= 0 && {
            let pos = record.pos();
            let end = record.cigar().end_pos();
            pos < win_end && end > win_start
        };
        if overlaps {
            match gap {
                Gap::NoReads => continue,
                Gap::ZeroMapq => record.set_mapq(0),
            }
        }
        writer.write(&record)?;
    }
    drop(writer);

    bam::index::build(output, None, bam::index::Type::Bai, 1)?;
    Ok(())
}

// TODO: add tests that compare default output with output when
// - mbias (nOT/nOB) are set
// - min depth is set
// - max depth is set
// - min bq is set
// - min mapq is set

/// Reads carrying one shared deletion, far enough into a contig that the
/// enclosing segment cannot start at 0.
///
/// `tests/data/test.bam` has 86 indel-carrying reads but no two of them agree
/// on an allele, so it produces no indel call at any threshold — which is why
/// nothing here exercised the indel path.
fn write_deletion_bam(
    output: &std::path::Path,
    contig: &str,
    del_pos: u64,
    del_len: usize,
    reads: usize,
) -> Result<()> {
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::{self, Record, header::HeaderRecord};
    use rust_htslib::faidx;

    const READ_LEN: usize = 80;
    // Every read spans the deletion but starts at a different offset, so the
    // anchor is reached from a different query position each time.
    let first_start = del_pos - 40;

    let fasta = faidx::Reader::from_path("tests/data/test.fasta.gz")?;
    let fetch = |from: u64, to: u64| -> Result<Vec<u8>> {
        Ok(fasta.fetch_seq(contig, from as usize, to as usize - 1)?.to_ascii_uppercase())
    };

    let mut header = bam::Header::new();
    let contig_len =
        fasta.fetch_seq_len(contig).ok_or_else(|| eyre!("contig {contig} not in the FASTA"))?;
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", contig);
    sq.push_tag(b"LN", contig_len);
    header.push_record(&sq);
    let mut writer = bam::Writer::from_path(output, &header, bam::Format::Bam)?;

    for i in 0..reads {
        let start = first_start + i as u64;
        let left = (del_pos - start) as usize;
        let right = READ_LEN - left;
        // The read skips the deleted bases, so its sequence is the reference
        // either side of them.
        let mut seq = fetch(start, del_pos + 1)?;
        let after = del_pos + 1 + del_len as u64;
        seq.extend_from_slice(&fetch(after, after + right as u64 - 1)?);
        seq.truncate(READ_LEN);

        let cigar = CigarString(
            vec![
                Cigar::Match(left as u32 + 1),
                Cigar::Del(del_len as u32),
                Cigar::Match(READ_LEN as u32 - left as u32 - 1),
            ]
            .into(),
        );
        let mut record = Record::new();
        record.set(format!("del_{i}").as_bytes(), Some(&cigar), &seq, &vec![40u8; seq.len()]);
        record.set_tid(0);
        record.set_pos(start as i64);
        record.set_mapq(60);
        record.set_flags(if i % 2 == 0 { 0 } else { 16 });
        record.set_mtid(-1);
        record.set_mpos(-1);
        writer.write(&record)?;
    }
    drop(writer);

    bam::index::build(output, None, bam::index::Type::Bai, 1)?;
    Ok(())
}

/// A deletion's REF allele must be the reference *at the deletion*, not at the
/// start of the enclosing segment; a wrong REF is real-looking sequence from
/// another locus, and `bcftools norm` refuses the file. This is the CLI test
/// that holds a multi-base REF.
///
/// The contig position is load-bearing. At a segment starting at 0 the
/// segment-relative and the genomic index coincide, which is also why the
/// `from_seqair` unit fixtures cannot see this.
#[test]
fn indel_ref_alleles_come_from_the_deletion_site() -> Result<()> {
    use rust_htslib::faidx;

    const CONTIG: &str = "bacteriophage_lambda_CpG";
    const DEL_POS: u64 = 24_137; // 0-based anchor, well past a segment boundary
    const DEL_LEN: usize = 4;

    let temp_dir = TempDir::new()?;
    let bam = temp_dir.path().join("deletion.bam");
    write_deletion_bam(&bam, CONTIG, DEL_POS, DEL_LEN, 12)?;

    let vcf = temp_dir.path().join("out.vcf");
    rastair()
        .args(["call", "--fasta-file=tests/data/test.fasta.gz"])
        .arg(&bam)
        .args([
            &format!("--region={CONTIG}:{}-{}", DEL_POS - 200, DEL_POS + 200),
            "--unpaired",
            "--vcf",
        ])
        .arg(&vcf)
        .succeeds()?;

    let fasta = faidx::Reader::from_path("tests/data/test.fasta.gz")?;
    let text = std::fs::read_to_string(&vcf)?;
    let mut multi_base_refs = 0;
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let fields: Vec<&str> = line.split('\t').collect();
        let (Some(&contig), Some(&pos), Some(&reference)) =
            (fields.first(), fields.get(1), fields.get(3))
        else {
            continue;
        };
        if reference.len() < 2 {
            continue;
        }
        multi_base_refs += 1;
        // VCF POS is 1-based; fetch_seq takes a 0-based inclusive range.
        let start = pos.parse::<usize>()? - 1;
        let expected = fasta.fetch_seq(contig, start, start + reference.len() - 1)?;
        assert_eq!(
            reference,
            String::from_utf8_lossy(&expected.to_ascii_uppercase()),
            "REF at {contig}:{pos} does not match the reference"
        );
    }

    assert!(
        multi_base_refs > 0,
        "no multi-base REF was emitted, so this test proved nothing:\n{text}"
    );
    Ok(())
}

/// A compound het, an insertion on one haplotype and a deletion on the other
/// at one anchor, is one `1/2` record whose REF spans the deletion; a deletion
/// on both haplotypes is a PASS `1/1`.
#[test]
fn the_pair_hmm_calls_a_compound_het_and_a_homozygous_indel() -> Result<()> {
    use utils::scenario::{Hap, Scenario};
    const COMPOUND: usize = 150;
    const HOMOZYGOUS: usize = 180;
    let scenario = Scenario::new(1000)
        .reference_at(COMPOUND, b"AGTC")
        .insertion(COMPOUND, b"GG", Hap::One)
        .deletion(COMPOUND, 2, Hap::Two)
        .reference_at(HOMOZYGOUS, b"CTAG")
        .deletion(HOMOZYGOUS, 2, Hap::One)
        .deletion(HOMOZYGOUS, 2, Hap::Two)
        .balanced(120, 370, 5);

    let vcf = call_scenario(&scenario, &[])?;
    let indels: Vec<_> = vcf
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let (reference, alt) = (*columns.get(3)?, *columns.get(4)?);
            (reference.len() > 1 || alt.len() > 1).then_some((
                columns.get(1)?.parse::<usize>().ok()?,
                reference,
                alt,
                *columns.get(6)?,
                format_field(line, "GT")?,
            ))
        })
        .collect();
    assert_eq!(
        indels,
        [
            (COMPOUND + 1, "AGT", "AGGGT,A", "PASS", "1/2"),
            (HOMOZYGOUS + 1, "CTA", "C", "PASS", "1/1")
        ],
        "{vcf}"
    );
    Ok(())
}

/// Two events of one cluster: two deletions 13 bp apart on both haplotypes
/// are each a PASS `1/1`, which needs a haplotype carrying both; a deletion
/// on one haplotype and an insertion 10 bp on, on the other, are each `0/1`.
#[test]
fn the_pair_hmm_calls_two_events_of_one_cluster() -> Result<()> {
    use utils::scenario::{Hap, Scenario};
    const FIRST: usize = 150;
    const SECOND: usize = FIRST + 13;
    const DELETION: usize = 400;
    const INSERTION: usize = DELETION + 10;
    let mut scenario = Scenario::new(1000)
        .reference_at(FIRST, b"CTAG")
        .reference_at(SECOND, b"AGTC")
        .reference_at(DELETION, b"AGTC")
        .reference_at(INSERTION, b"AT");
    for hap in [Hap::One, Hap::Two] {
        scenario = scenario.deletion(FIRST, 2, hap).deletion(SECOND, 2, hap);
    }
    let scenario = scenario
        .deletion(DELETION, 2, Hap::One)
        .insertion(INSERTION, b"GG", Hap::Two)
        .balanced(120, 370, 5)
        .balanced(370, 620, 5);

    let vcf = call_scenario(&scenario, &[])?;
    let indels: Vec<_> = vcf
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let (reference, alt) = (*columns.get(3)?, *columns.get(4)?);
            (reference.len() > 1 || alt.len() > 1).then_some((
                columns.get(1)?.parse::<usize>().ok()?,
                reference,
                alt,
                *columns.get(6)?,
                format_field(line, "GT")?,
            ))
        })
        .collect();
    assert_eq!(
        indels,
        [
            (FIRST + 1, "CTA", "C", "PASS", "1/1"),
            (SECOND + 1, "AGT", "A", "PASS", "1/1"),
            (DELETION + 1, "AGT", "A", "PASS", "0|1"),
            (INSERTION + 1, "A", "AGG", "PASS", "1|0"),
        ],
        "{vcf}"
    );
    Ok(())
}

/// SNV events never change an indel call: the scenario of
/// `the_pair_hmm_calls_two_events_of_one_cluster` with het SNVs inside both
/// clusters gives byte-identical indel records with `--snv-engine` on and
/// off, though with the bundled model the engine scores those SNVs against
/// the indels' haplotypes.
#[test]
fn snv_events_leave_the_indel_records_alone() -> Result<()> {
    use utils::scenario::{Hap, Scenario};
    const FIRST: usize = 150;
    const SECOND: usize = FIRST + 13;
    const DELETION: usize = 400;
    const INSERTION: usize = DELETION + 10;
    let mut scenario = Scenario::new(1000)
        .reference_at(FIRST, b"CTAG")
        .reference_at(SECOND, b"AGTC")
        .reference_at(DELETION, b"AGTC")
        .reference_at(INSERTION, b"AT");
    for hap in [Hap::One, Hap::Two] {
        scenario = scenario.deletion(FIRST, 2, hap).deletion(SECOND, 2, hap);
    }
    let scenario = scenario
        .deletion(DELETION, 2, Hap::One)
        .insertion(INSERTION, b"GG", Hap::Two)
        .snv(FIRST - 6, Hap::One)
        .snv(FIRST + 8, Hap::Two)
        .snv(SECOND + 9, Hap::One)
        .snv(DELETION + 6, Hap::Two)
        .snv(INSERTION + 7, Hap::One)
        .balanced(120, 370, 5)
        .balanced(370, 620, 5);

    let dir = TempDir::new()?;
    let fasta = scenario.write_fasta(dir.path())?;
    let bam = dir.path().join("reads.bam");
    scenario.write_bam(&bam)?;
    let call = |engine: bool| -> Result<Vec<String>> {
        let vcf = dir.path().join(format!("out-{engine}.vcf"));
        rastair()
            .arg("call")
            .arg(format!("--snv-engine={engine}"))
            .arg("--fasta-file")
            .arg(&fasta)
            .arg(&bam)
            .args([scenario.region().as_str(), "--all", "--vcf"])
            .arg(&vcf)
            .succeeds()?;
        Ok(std::fs::read_to_string(&vcf)?
            .lines()
            .filter(|line| !line.starts_with('#'))
            .filter(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                columns.get(3).is_some_and(|r| r.len() > 1)
                    || columns.get(4).is_some_and(|a| a.split(',').any(|a| a.len() > 1))
            })
            .map(str::to_owned)
            .collect())
    };
    let (off, on) = (call(false)?, call(true)?);
    assert_eq!(off.len(), 4, "{off:#?}");
    assert_eq!(on, off);
    Ok(())
}

/// Indel calling end to end on a scenario: a het insertion and a het
/// deletion are called phased, `0|1`, PASS, with `QUAL` of at least 10 and a `GQ`.
/// An insertion only two of forty-two fragments carry is scored and
/// rejected: absent from the output, and written `0/0` with the
/// `low_indel_qual` FILTER under `--all`.
#[test]
fn the_pair_hmm_calls_indels_on_a_scenario() -> Result<()> {
    use utils::scenario::{Hap, Scenario, Strand};
    const INSERTION: usize = 150;
    const DELETION: usize = 180;
    const REJECTED: usize = 600;
    let scenario = Scenario::new(1000)
        .reference_at(INSERTION, b"AT")
        .insertion(INSERTION, b"GG", Hap::Two)
        .reference_at(DELETION, b"AGTC")
        .deletion(DELETION, 2, Hap::Two)
        .balanced(120, 370, 5)
        .reference_at(REJECTED, b"CA")
        .insertion(REJECTED, b"TT", Hap::Two)
        // Two carriers in forty fragments: exactly the 5 % of the depth a
        // candidate needs, so the pair-HMM scores and rejects it.
        .pairs(Hap::One, Strand::Top, 560, 19)
        .pairs(Hap::One, Strand::Bottom, 560, 19)
        .pairs(Hap::Two, Strand::Top, 560, 1)
        .pairs(Hap::Two, Strand::Bottom, 560, 1);

    // (POS, REF, ALT, QUAL, FILTER, GT, GQ) of every indel record.
    let indels = |vcf: &str| -> Vec<(usize, String, String, f32, String, String, String)> {
        vcf.lines()
            .filter(|line| !line.starts_with('#'))
            .filter_map(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                let (reference, alt) = (*columns.get(3)?, *columns.get(4)?);
                if reference.len() == 1 && alt.len() == 1 {
                    return None;
                }
                Some((
                    columns.get(1)?.parse().ok()?,
                    reference.to_owned(),
                    alt.to_owned(),
                    columns.get(5)?.parse().ok()?,
                    (*columns.get(6)?).to_owned(),
                    format_field(line, "GT")?.to_owned(),
                    format_field(line, "GQ").unwrap_or(".").to_owned(),
                ))
            })
            .collect()
    };

    let vcf = call_scenario(&scenario, &["--vcf-all-fields"])?;
    let called = indels(&vcf);
    let sites: Vec<_> = called
        .iter()
        .map(|(pos, r, a, _, f, gt, _)| (*pos, r.as_str(), a.as_str(), f.as_str(), gt.as_str()))
        .collect();
    assert_eq!(
        sites,
        [(INSERTION + 1, "A", "AGG", "PASS", "0|1"), (DELETION + 1, "AGT", "A", "PASS", "0|1")],
        "{vcf}"
    );
    for (pos, .., qual, _, _, gq) in &called {
        assert!(*qual >= 10.0 && gq.parse::<i32>().is_ok_and(|gq| gq > 0), "{pos}: {vcf}");
    }
    // AD counts the reads the pair-HMM assigns to each allele, DP the reads it
    // scored: one per fragment, so twenty pairs over each site.
    let depths: Vec<_> = vcf
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            let info = line.split('\t').nth(7)?;
            let field =
                |key: &str| info.split(';').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='));
            let pos: usize = line.split('\t').nth(1)?.parse().ok()?;
            [INSERTION + 1, DELETION + 1].contains(&pos).then(|| (pos, field("AD"), field("DP")))
        })
        .collect();
    assert_eq!(
        depths,
        [(INSERTION + 1, Some("10,10"), Some("20")), (DELETION + 1, Some("10,10"), Some("20"))],
        "{vcf}"
    );

    let vcf = call_scenario(&scenario, &["--vcf-all-fields", "--all"])?;
    let rejected: Vec<_> =
        indels(&vcf).into_iter().filter(|record| record.0 == REJECTED + 1).collect();
    let [(_, reference, alt, qual, filter, gt, _)] = rejected.as_slice() else {
        panic!("one record at the rejected insertion: {vcf}")
    };
    assert_eq!(
        (reference.as_str(), alt.as_str(), filter.as_str(), gt.as_str()),
        ("C", "CTT", "low_indel_qual", "0/0")
    );
    assert!(*qual < 10.0, "{vcf}");
    Ok(())
}

/// CpG-only output is the CpG subset of the full output: a CpG two SNVs form
/// together (`TA>CG`), which no reference neighbour shows, is in the BED of
/// `-c` exactly as in the BED a full run writes beside its VCF.
#[test]
fn cpg_only_output_reports_a_cpg_of_two_snvs() -> Result<()> {
    use utils::scenario::{Hap, Scenario};
    const AT: usize = 300;
    let scenario = Scenario::new(600)
        .reference_at(AT, b"TA")
        .het(AT, b'C', Hap::One)
        .het(AT + 1, b'G', Hap::One)
        .methylated(Hap::One)
        .balanced(250, 400, 8);
    let dir = TempDir::new()?;
    let bed = |name: &str, extra: &[&str]| -> Result<String> {
        let path = dir.path().join(name);
        let path_str = path.to_str().ok_or_else(|| eyre!("temp path is not UTF-8"))?;
        let args: Vec<&str> = extra.iter().copied().chain(["--bed", path_str]).collect();
        call_scenario_in(dir.path(), &scenario, &args)?;
        Ok(std::fs::read_to_string(&path)?)
    };
    let full = bed("full.bed", &[])?;
    let cpgs_only = bed("cpgs.bed", &["-c"])?;
    let positions: Vec<_> = full
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.split('\t').nth(1)?.parse::<usize>().ok())
        .collect();
    assert_eq!(positions, [AT, AT + 1], "{full}");
    assert_eq!(cpgs_only, full);
    Ok(())
}

/// `convert` writes the BED `call` writes beside the VCF: at a heterozygous
/// `C>T` of a reference CpG the remaining `C` keeps its β, and under `--all`
/// the partner `G`, whose methylated reads are a rejected `A` record after its
/// main record, is still one row.
#[test]
fn convert_writes_the_bed_of_call_at_a_het_cpg_snv() -> Result<()> {
    use utils::scenario::{Hap, Scenario};
    const AT: usize = 300;
    let scenario = Scenario::new(600)
        .reference_at(AT, b"CG")
        .het(AT, b'T', Hap::One)
        .methylated(Hap::Two)
        .balanced(250, 400, 8);
    for extra in [&[][..], &["--all"]] {
        let dir = TempDir::new()?;
        let called = dir.path().join("called.bed");
        let converted = dir.path().join("converted.bed");
        let called_str = called.to_str().ok_or_else(|| eyre!("temp path is not UTF-8"))?;
        let args: Vec<&str> = extra.iter().copied().chain(["--bed", called_str]).collect();
        let vcf = call_scenario_in(dir.path(), &scenario, &args)?;
        rastair()
            .args(["convert", "--input"])
            .arg(dir.path().join("out.vcf"))
            .arg("--output")
            .arg(&converted)
            .succeeds()?;
        let called = std::fs::read_to_string(&called)?;
        assert!(called.contains("C/T"), "{called}\n{vcf}");
        assert_eq!(std::fs::read_to_string(&converted)?, called, "{extra:?}\n{vcf}");
    }
    Ok(())
}

/// CpG-only VCF output holds the CpG records of a full run, and not the
/// indel record anchored on a CpG's `G`.
#[test]
fn cpg_only_vcf_is_the_cpg_records_of_a_full_run() -> Result<()> {
    use utils::scenario::{Hap, Scenario};
    const AT: usize = 300;
    let mut scenario = Scenario::new(600).reference_at(AT, b"CGTCA");
    for hap in [Hap::One, Hap::Two] {
        scenario = scenario.deletion(AT + 1, 2, hap);
    }
    let scenario = scenario.balanced(250, 400, 8);
    let records = |vcf: &str| -> Vec<String> {
        vcf.lines().filter(|line| !line.starts_with('#')).map(String::from).collect()
    };
    let is_cpg = |line: &&String| {
        line.split('\t')
            .nth(7)
            .is_some_and(|info| info.split(';').any(|key| key == "CPG" || key == "CPGnovo"))
    };
    let full = records(&call_scenario(&scenario, &[])?);
    assert!(full.iter().any(|line| line.contains("\tGTC\tG\t")), "{full:#?}");
    let cpg_records: Vec<_> = full.iter().filter(is_cpg).cloned().collect();
    assert_eq!(cpg_records.len(), 2, "{full:#?}");
    assert_eq!(records(&call_scenario(&scenario, &["-c"])?), cpg_records);
    Ok(())
}
