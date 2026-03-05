#![expect(non_snake_case, reason = "readable test names")]

mod utils;
use insta::assert_compact_debug_snapshot;
use utils::*;

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
        .arg("--read-groups=mTet1-PyBr-16h-p1_S1_L001")
        .succeeds()?;

    let unfiltered_count = vcf_content_lines(&std::fs::read_to_string(&unfiltered)?).count();
    let filtered_count = vcf_content_lines(&std::fs::read_to_string(&filtered)?).count();

    assert!(
        filtered_count < unfiltered_count,
        "single-group filter should produce fewer records (filtered: {filtered_count}, unfiltered: {unfiltered_count})",
    );

    Ok(())
}

#[test]
fn read_group_multiple_space_separated_filters_give_intermediate_records() -> Result<()> {
    apply_common_filters!();

    const REGION: &str = "--region=chr19:6103000-6110000";

    let temp_dir = TempDir::new()?;
    let one_group = temp_dir.path().join("one_group.vcf");
    let two_groups = temp_dir.path().join("two_groups.vcf");
    let unfiltered = temp_dir.path().join("unfiltered.vcf");

    rastair().args(CALL_TEST_BAM).args([REGION, NO_ML, "--vcf"]).arg(&unfiltered).succeeds()?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([REGION, NO_ML, "--vcf"])
        .arg(&one_group)
        .arg("--read-groups=mTet1-PyBr-16h-p1_S1_L001")
        .succeeds()?;

    // Space-separated: pass both groups as a single --read-groups value list
    rastair()
        .args(CALL_TEST_BAM)
        .args([REGION, NO_ML, "--vcf"])
        .arg(&two_groups)
        .arg("--read-groups")
        .arg("mTet1-PyBr-16h-p1_S1_L001")
        .arg("mTet1-PyBr-16h-p1_S1_L002")
        .succeeds()?;

    let unfiltered_count = vcf_content_lines(&std::fs::read_to_string(&unfiltered)?).count();
    let one_count = vcf_content_lines(&std::fs::read_to_string(&one_group)?).count();
    let two_count = vcf_content_lines(&std::fs::read_to_string(&two_groups)?).count();

    assert!(
        one_count <= two_count,
        "two groups should yield at least as many records as one (one: {one_count}, two: {two_count})",
    );
    assert!(
        two_count <= unfiltered_count,
        "two groups should yield no more records than unfiltered (two: {two_count}, unfiltered: {unfiltered_count})",
    );

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

// TODO: add tests that compare default output with output when
// - mbias (nOT/nOB) are set
// - min depth is set
// - max depth is set
// - min bq is set
// - min mapq is set
