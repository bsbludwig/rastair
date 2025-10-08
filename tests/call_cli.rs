mod utils;
use utils::*;

#[test]
fn simple_call_gives_you_vcf_on_stdout() -> Result<()> {
    apply_common_filters!();

    let output = rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--thresholds", // disable ML for faster test
            "--region=chr19:6105700-6105800",
        ])
        .output()?;

    let stderr = String::from_utf8(output.stderr).wrap_err("utf8 decode")?;
    assert_snapshot!(stderr, @r#"
    [TIME] INFO rastair::call: Wrote VCF output file="-"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    let stdout = String::from_utf8(output.stdout).wrap_err("utf8 decode")?;
    assert!(stdout.trim().starts_with("##fileformat=VCF"));

    Ok(())
}

#[test]
fn asking_for_cpgs_defaults_to_bed_output() -> Result<()> {
    apply_common_filters!();

    let output = rastair()
        .args([
            "call",
            "-c",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--thresholds", // disable ML for faster test
            "--region=chr19:6105700-6105800",
        ])
        .output()?;

    let stderr = String::from_utf8(output.stderr).wrap_err("utf8 decode")?;
    assert_snapshot!(stderr, @r#"
    [TIME] INFO rastair::call: Wrote BED output file="-"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    let stdout = String::from_utf8(output.stdout).wrap_err("utf8 decode")?;
    assert!(stdout.trim().starts_with("#chr"));

    Ok(())
}

#[test]
fn writing_vcf_to_file() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.vcf");

    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "--thresholds", // disable ML for faster test
        "--region=chr19:6105700-6105800",
        "--vcf",
    ]).arg(&temp_file), @r#"
    success: true
    exit_code: 0
    ----- stdout -----

    ----- stderr -----
    [TIME] INFO rastair::call: Wrote VCF output file=[PATH]"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    Ok(())
}

#[test]
fn ask_for_cpgs_and_vcf() -> Result<()> {
    apply_common_filters!();

    let output = rastair()
        .args([
            "call",
            "-c",
            "--vcf",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--thresholds", // disable ML for faster test
            "--region=chr19:6105700-6105800",
        ])
        .output()?;

    let stderr = String::from_utf8(output.stderr).wrap_err("utf8 decode")?;
    assert_snapshot!(stderr, @r#"
    [TIME] INFO rastair::call: Wrote VCF output file="-"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    let stdout = String::from_utf8(output.stdout).wrap_err("utf8 decode")?;
    assert!(stdout.trim().starts_with("##fileformat=VCF"));

    Ok(())
}

#[test]
fn write_bcf_to_file_and_bed_to_stdout() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.bcf");

    let call = rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--thresholds", // disable ML for faster test
            "--region=chr19:6105700-6105750",
            "--bed=-",
            "--vcf",
        ])
        .arg(&temp_file)
        .output()?;

    assert!(call.status.success(), "rastair call failed");
    assert!(temp_file.exists());

    let bed = str::from_utf8(&call.stdout)?;
    assert_snapshot!(bed);

    Ok(())
}

#[test]
fn when_asked_for_bed_file_in_vcf_param_we_are_nice() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.bed");

    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "--thresholds", // disable ML for faster test
        "--region=chr19:6105700-6105800",
        "--vcf",
    ]).arg(&temp_file), @r#"
    success: true
    exit_code: 0
    ----- stdout -----

    ----- stderr -----
    [TIME] WARN rastair::call: VCF output file name ends with `.bed`/`.bed.gz`, did you mean to use `--bed` instead of `-o`/`--vcf`? Assuming you meant `--bed` and switching the output accordingly. file=[PATH]"
    [TIME] INFO rastair::call: Wrote BED output file=[PATH]"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    Ok(())
}

#[test]
fn includes_all_cpgs_when_methylation_calling() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(
        "does not include unmethylated cpgs without alts",
        rastair().args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--thresholds", // disable ML for faster test
            "--region=chr19:6117965-6118004",
            "--skip-methylation-calling",
            "--vcf"
        ])
    );

    assert_cmd_snapshot!(
        "includes all cpgs",
        rastair().args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--thresholds", // disable ML for faster test
            "--region=chr19:6117965-6118004",
            "--vcf"
        ])
    );

    Ok(())
}

#[test]
fn segmentation_overlaps_do_not_cause_duplicate_records() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.vcf");

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--thresholds", // disable ML for faster test
            "--region=chr19",
            "--segment-max-length=10000",
            "--segment-overlap=300",
            "--vcf",
        ])
        .arg(&temp_file)
        .succeeds()
        .wrap_err("rastair call failed")?;

    let text = std::fs::read_to_string(&temp_file).wrap_err("read rastair 2 vcf")?;
    text.lines()
        .filter(|line| !line.starts_with("#"))
        .filter_map(|line| line.split("\t").nth(1))
        .filter_map(|x| x.parse::<u32>().ok())
        .try_fold(BTreeSet::new(), |mut set, position| {
            let is_new = set.insert(position);
            if is_new { Ok(set) } else { Err(eyre!("Duplicate position found: {}", position)) }
        })?;

    Ok(())
}
