mod utils;
use utils::*;

#[test]
fn random_call() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.vcf");

    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "--region=chr19:6105700-6105800",
        "--vcf",
    ]).arg(&temp_file), @r#"
    success: true
    exit_code: 0
    ----- stdout -----

    ----- stderr -----
    [TIME] INFO rastair2::call: Wrote VCF output file=[PATH]"
    [TIME] INFO rastair2: Call finished [DURATION]
    "#);

    Ok(())
}

#[test]
fn pipe_to_stdout() -> Result<()> {
    apply_common_filters!();

    // Please manually verify that only the VCF output is printed to stdout and
    // logging goes to stderr.
    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "--region=chr19:6105700-6105750",
        "--vcf",
    ]));

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
fn includes_all_cpgs_when_methylation_calling() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(
        "does not include unmethylated cpgs without alts",
        rastair().args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
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
            "--region=chr19:6117965-6118004",
            "--vcf"
        ])
    );

    Ok(())
}

#[test]
fn includes_only_cpgs_when_methylation_calling() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(
        "includes only cpgs",
        rastair().args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--region=chr19:6117965-6118004",
            "--cpgs-only",
            "--vcf",
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
            "--region=chr19",
            "--segment-max-length=10000",
            "--segment-overlap=300",
            "--vcf",
        ])
        .arg(&temp_file)
        .succeeds()
        .wrap_err("rastair2 call failed")?;

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
