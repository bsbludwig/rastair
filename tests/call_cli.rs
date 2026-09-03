#![expect(non_snake_case, reason = "readable test names")]

mod utils;
use insta::assert_compact_debug_snapshot;
use std::collections::{BTreeMap, BTreeSet};
use utils::*;

#[derive(Debug, Clone, PartialEq)]
struct CpgBedCall {
    pos: u32,
    strand: char,
    beta_est: f64,
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
    flag_beta: f64,
    evidence_beta: f64,
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
        let beta_est: f64 = fields[4].parse().wrap_err("parse beta_est")?;
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
        let beta_abs_diff = (flag.beta_est - evidence.beta_est).abs();
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

    let call = rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, NO_ML]).output()?;

    assert_snapshot!(call.stderr(), @r#"
    [TIME] INFO rastair::call::writer: Wrote VCF output file="-"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    let stdout = call.stdout();
    assert!(stdout.trim().starts_with("##fileformat=VCF"));
    assert_snapshot!(stdout);

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

    let mut call = rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--guess-read-orientation"])
        .output()?;

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
        .args([REGION, NO_ML, "--cpgs-only", "--bed"])
        .arg(&flag_bed)
        .succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([REGION, NO_ML, "--cpgs-only", "--guess-read-orientation", "--bed"])
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

#[test]
fn asking_for_cpgs_defaults_to_bed_output() -> Result<()> {
    apply_common_filters!();

    let call = rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, NO_ML]).arg("-c").output()?;

    assert_snapshot!(call.stderr(), @r#"
    [TIME] INFO rastair::call::writer: Wrote BED output file="-"
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
        rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, NO_ML]).arg("--vcf").arg(&temp_file), @r#"
    success: true
    exit_code: 0
    ----- stdout -----

    ----- stderr -----
    [TIME] INFO rastair::call::writer: Wrote VCF output file=[PATH]"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    read_bcf(&temp_file).wrap_err("validate bcf file")?;

    Ok(())
}

#[test]
fn asking_for_all_variants_includes_non_passing_ones() -> Result<()> {
    apply_common_filters!();

    let mut call =
        rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, NO_ML]).arg("--vcf").output()?;
    call.succeeds()?;
    assert!(!call.stdout().lines().filter(|l| !l.starts_with("#")).any(|l| l.contains("lowDp")));

    let mut call = rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML])
        .arg("--vcf")
        .arg("--all")
        .output()?;
    call.succeeds()?;
    let stdout = call.stdout();
    assert!(vcf_content_lines(&stdout).any(|l| l.contains("lowDp")));

    Ok(())
}

#[test]
fn ask_for_cpgs_and_vcf() -> Result<()> {
    apply_common_filters!();

    let call =
        rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, NO_ML]).args(["-c", "--vcf"]).output()?;

    assert_snapshot!(call.stderr(), @r#"
    [TIME] INFO rastair::call::writer: Wrote VCF output file="-"
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
        .args([CHR19_SMALL, NO_ML])
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

#[test]
fn when_asked_for_bed_file_in_vcf_param_we_are_nice() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.bed");

    assert_cmd_snapshot!(
        rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, NO_ML]).arg("--vcf").arg(&temp_file),
        @r#"
    success: true
    exit_code: 0
    ----- stdout -----

    ----- stderr -----
    [TIME] WARN rastair::call: VCF output file name ends with `.bed`/`.bed.gz`, did you mean to use `--bed` instead of `-o`/`--vcf`? Assuming you meant `--bed` and switching the output accordingly. file=[PATH]"
    [TIME] INFO rastair::call::writer: Wrote BED output file=[PATH]"
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
        .args([NO_ML, REGION])
        .args(["--segment-max-length=1000", "--segment-overlap=100", "--threads=7", "--bed"])
        .arg(&temp_file1)
        .succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([NO_ML, REGION])
        .args(["--segment-max-length=1001", "--segment-overlap=1", "--threads=4", "--bed"])
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
        .args([NO_ML, REGION])
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

    rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, NO_ML]).arg("--vcf").arg(&a).succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--nOT=12,12,12,12", "--nOB=12,12,12,12"])
        .arg("--vcf")
        .arg(&b)
        .succeeds()?;

    assert_compact_debug_snapshot!(get_depths(&a), @"Ok([18, 16, 16])");
    assert_compact_debug_snapshot!(get_depths(&b), @"Ok([11, 13, 13])");

    fn get_depths(path: &std::path::Path) -> Result<Vec<i32>> {
        use rastair_vcf::VcfField as _;
        use rust_htslib::bcf::Read;

        let mut bcf = read_bcf(path).wrap_err("invalid bcf file")?;
        let depths = bcf
            .records()
            .map(|r| {
                let field = r
                    .unwrap()
                    .info(rastair_vcf::standard_fields::ReadDepth::ID.as_bytes())
                    .integer()
                    .unwrap()
                    .unwrap();
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

    // Create VCF with default fields
    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
        .arg(&default_vcf)
        .succeeds()?;

    // Create VCF with additional fields enabled
    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
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
fn min_depth_filter_reduces_variant_calls() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let default_vcf = temp_dir.path().join("default.vcf");
    let filtered_vcf = temp_dir.path().join("filtered.vcf");

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
        .arg(&default_vcf)
        .succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
        .arg(&filtered_vcf)
        .arg("--v-min-depth=10")
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
fn max_coverage_filter_affects_output() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let default_vcf = temp_dir.path().join("default.vcf");
    let filtered_vcf = temp_dir.path().join("filtered.vcf");

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
        .arg(&default_vcf)
        .succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
        .arg(&filtered_vcf)
        .arg("--m-max-coverage=5")
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
fn min_baseq_filter_reduces_variant_calls() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let default_vcf = temp_dir.path().join("default.vcf");
    let filtered_vcf = temp_dir.path().join("filtered.vcf");

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
        .arg(&default_vcf)
        .succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
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

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
        .arg(&default_vcf)
        .succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
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

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
        .arg(&default_vcf)
        .succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
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
        .args([CHR19_SMALL, NO_ML, "--vcf"])
        .arg(&vcf_miseq)
        .arg("--error-model=miseq")
        .succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
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
        .args([CHR19_SMALL, NO_ML, "--vcf"])
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
        .args([CHR19_SMALL, NO_ML, "--vcf"])
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

    rastair().args(CALL_TEST_BAM).args([REGION, NO_ML, "--vcf"]).arg(&unfiltered).succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([REGION, NO_ML, "--vcf"])
        .arg(&filtered)
        .arg("--require-tags=RG=mTet1-PyBr-16h-p1_S1_L001")
        .succeeds()?;

    let unfiltered_count = vcf_content_lines(&std::fs::read_to_string(&unfiltered)?).count();
    let filtered_count = vcf_content_lines(&std::fs::read_to_string(&filtered)?).count();

    assert_snapshot!(unfiltered_count, @"192");
    assert_snapshot!(filtered_count, @"180");

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
            .args([REGION, NO_ML, "--vcf"])
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

    assert_snapshot!(unfiltered, @"192");
    assert_snapshot!(lane1_only, @"180");

    Ok(())
}

#[test]
fn error_model_rejects_invalid_platform_name() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let vcf = temp_dir.path().join("test.vcf");

    let result = rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
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

    let temp_dir = TempDir::new()?;
    let gapped_bam = temp_dir.path().join("gapped.bam");
    // Set MAPQ 0 on every read covering the second reference CpG pair in the
    // region (0-based 6105743..6105745, i.e. 1-based 6105744). The reads still
    // align there, so the pileup column exists, but with the default
    // `--min-mapq 1` they are filtered out and the position has zero coverage.
    write_bam_with_zero_mapq_overlapping(&gapped_bam, "chr19", 6_105_743, 6_105_745)?;

    let cpgs_only_bed = temp_dir.path().join("cpgs-only.bed");
    let all_bed = temp_dir.path().join("all.bed");

    rastair()
        .args(["call", "--fasta-file=tests/data/test.fasta.gz"])
        .arg(&gapped_bam)
        .args([REGION, NO_ML, "--cpgs-only", "--bed"])
        .arg(&cpgs_only_bed)
        .succeeds()?;

    rastair()
        .args(["call", "--fasta-file=tests/data/test.fasta.gz"])
        .arg(&gapped_bam)
        .args([REGION, NO_ML, "--cpgs-only", "--all", "--bed"])
        .arg(&all_bed)
        .succeeds()?;

    let cpgs_only_positions: BTreeSet<(u32, char)> =
        parse_cpg_bed(&cpgs_only_bed)?.keys().copied().collect();
    let all_positions: BTreeSet<(u32, char)> = parse_cpg_bed(&all_bed)?.keys().copied().collect();

    // The second CpG pair has no coverage, so plain `--cpgs-only` drops it …
    assert_eq!(
        cpgs_only_positions,
        BTreeSet::from([(6_105_711, '+'), (6_105_712, '-')]),
        "expected only the covered CpG pair under --cpgs-only"
    );
    // … but `--all --cpgs-only` must report every reference CpG in the region,
    // whether it has coverage or not (both CpG pairs here).
    assert_eq!(
        all_positions,
        BTreeSet::from([(6_105_711, '+'), (6_105_712, '-'), (6_105_743, '+'), (6_105_744, '-'),]),
        "expected all reference CpGs under --all --cpgs-only"
    );

    Ok(())
}

#[test]
fn convert_bed_include_empty_gates_uncovered_cpgs_independently_of_call() -> Result<()> {
    apply_common_filters!();

    const REGION: &str = "--region=chr19:6105700-6105800";

    let temp_dir = TempDir::new()?;
    let gapped_bam = temp_dir.path().join("gapped.bam");
    write_bam_with_zero_mapq_overlapping(&gapped_bam, "chr19", 6_105_743, 6_105_745)?;

    // `call --all --cpgs-only` already reports the uncovered CpG pair with no
    // extra flag needed (that's the point of this refactor), so the mpk file
    // contains both pairs.
    let mpk = temp_dir.path().join("all.mpk.lz4");
    rastair()
        .args(["call", "--fasta-file=tests/data/test.fasta.gz"])
        .arg(&gapped_bam)
        .args([REGION, NO_ML, "--cpgs-only", "--all", "-o"])
        .arg(&mpk)
        .succeeds()?;

    let without_empty = temp_dir.path().join("without-empty.bed");
    let with_empty = temp_dir.path().join("with-empty.bed");

    rastair()
        .args(["convert", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(&without_empty)
        .succeeds()?;

    rastair()
        .args(["convert", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(&with_empty)
        .arg("--bed-include-empty")
        .succeeds()?;

    let without_empty_positions: BTreeSet<(u32, char)> =
        parse_cpg_bed(&without_empty)?.keys().copied().collect();
    let with_empty_positions: BTreeSet<(u32, char)> =
        parse_cpg_bed(&with_empty)?.keys().copied().collect();

    // `convert` still has to decide independently whether to carry the
    // uncovered pair, already present in the mpk, into this BED output.
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

/// Copy `tests/data/test.bam`, setting MAPQ 0 on all reads whose alignment
/// overlaps the 0-based half-open window `[win_start, win_end)` on `chrom`.
fn write_bam_with_zero_mapq_overlapping(
    output: &std::path::Path,
    chrom: &str,
    win_start: i64,
    win_end: i64,
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
            record.set_mapq(0);
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
