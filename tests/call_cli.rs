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
#[ignore = "TODO: Fix after changing vcf output"]
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

// TODO: add tests that compare default output with output when
// - mbias (nOT/nOB) are set
// - min depth is set
// - max depth is set
// - min bq is set
// - min mapq is set
