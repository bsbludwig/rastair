mod utils;
use utils::*;

#[test]
fn missing_bam() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=test_data/test.fasta.gz",
        "/path/to/nonexistent/file.bam",
        "--vcf"
    ]), @r#"
    success: false
    exit_code: 2
    ----- stdout -----

    ----- stderr -----
    error: Invalid value for --fasta-file <FASTA_FILE>: Invalid path "test_data/test.fasta.gz": No such file or directory (os error 2)

    Usage: rastair call [OPTIONS] --fasta-file <FASTA_FILE> <BAM_FILE>

    For more information, try '--help'.
    "#);

    Ok(())
}

#[test]
fn missing_fasta() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=tests/data/test_which_doesnt_exist.fasta",
        "tests/data/test.bam"
    ]), @r#"
    success: false
    exit_code: 2
    ----- stdout -----

    ----- stderr -----
    error: Invalid value for --fasta-file <FASTA_FILE>: Invalid path "tests/data/test_which_doesnt_exist.fasta": No such file or directory (os error 2)

    Usage: rastair call [OPTIONS] --fasta-file <FASTA_FILE> <BAM_FILE>

    For more information, try '--help'.
    "#);

    Ok(())
}

#[test]
fn validates_region_arg() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "--region=chr19:6105700-xxx",
        "--vcf",
    ]), @"
    success: false
    exit_code: 2
    ----- stdout -----

    ----- stderr -----
    error: invalid value 'chr19:6105700-xxx' for '--region <REGIONS>': Malformed region string

    For more information, try '--help'.
    ");

    Ok(())
}

#[test]
fn different_paths_for_bed_and_vcf() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "--region=chr19:6105700",
        "--vcf=foo",
        "--bed=foo",
    ]), @"
    success: false
    exit_code: 1
    ----- stdout -----

    ----- stderr -----
    Error: 
       0: Unclear output choice
       1: Can't write both VCF and BED output to the same file. Please specify different output files.
    ");

    Ok(())
}

#[test]
fn invalid_bam_file() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let fake_bam = temp_dir.path().join("test.bam");
    std::fs::write(&fake_bam, b"This is not a BAM file")?;

    let mut cmd = rastair();
    let mut cmd = cmd
        .args(["call", "--region=chr19:6105700-6105800", "--fasta-file=tests/data/test.fasta.gz"])
        .arg(&fake_bam);

    #[cfg(not(feature = "experimental-seqair"))]
    assert_cmd_snapshot!(cmd, @r#"
    success: false
    exit_code: 1
    ----- stdout -----

    ----- stderr -----
    Error: 
       0: Failed to read BAM/FASTA files
       1: unable to open SAM/BAM/CRAM file at [PATH]

    Suggestion: Ensure the BAM/CRAM file is sorted and indexed with `samtools sort "[PATH]"`, respectively.
    Note: If you have a .bai/.crai file, ensure it is in the same directory as the BAM/CRAM file.
    "#);

    #[cfg(feature = "experimental-seqair")]
    assert_cmd_snapshot!(cmd, @"
    success: false
    exit_code: 1
    ----- stdout -----

    ----- stderr -----
    Error: 
       0: Failed to read BAM/FASTA files
       1: Failed to open alignment file [PATH]
       2: unrecognized file format for [PATH]), bgzf-compressed SAM (.sam.gz), CRAM (.cram).
    ");

    Ok(())
}

#[test]
fn vcf_field_configuration_invalid_field_id() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.vcf");

    // Try to enable an invalid INFO field
    let result = rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
        .arg(&temp_file)
        .arg("--vcf-info-fields=INVALID_FIELD")
        .output()?;

    assert!(!result.status.success(), "Should fail with invalid field ID");
    assert_snapshot!(result.stderr());

    // Try to enable an invalid FORMAT field
    let result = rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--vcf"])
        .arg(&temp_file)
        .arg("--vcf-format-fields=INVALID_FIELD")
        .output()?;

    assert!(!result.status.success(), "Should fail with invalid field ID");
    assert_snapshot!(result.stderr());

    Ok(())
}
