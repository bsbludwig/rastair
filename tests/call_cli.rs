mod utils;
use utils::*;

#[test]
fn random_call() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.vcf");

    assert_cmd_snapshot!(rastair().args([
        "call",
        "-r",
        "tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "-l",
        "chr19:6105700-6105800",
        "--calling=thresholds",
        "-o",
    ]).arg(&temp_file), @r#"
    success: true
    exit_code: 0
    ----- stdout -----

    ----- stderr -----
    [TIME] INFO rastair2::call: Wrote output file=[PATH]"
    [TIME] INFO rastair2: Call finished [DURATION]
    "#);

    let result = std::fs::read_to_string(&temp_file).wrap_err("read temp file")?;
    assert_snapshot!(result);

    Ok(())
}

#[test]
fn pipe_to_stdout() -> Result<()> {
    apply_common_filters!();

    // Please manually verify that only the VCF output is printed to stdout and
    // logging goes to stderr.
    assert_cmd_snapshot!(rastair().args([
        "call",
        "-r",
        "tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "-l",
        "chr19:6105700-6105750",
    ]));

    Ok(())
}

#[test]
fn includes_all_cpgs_when_methylation_calling() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(
        "does not include unmethylated cpgs without alts",
        rastair().args([
            "call",
            "-r",
            "tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "-l",
            "chr19:6117965-6118004",
            "-o",
            "-"
        ])
    );

    assert_cmd_snapshot!(
        "includes all cpgs",
        rastair().args([
            "call",
            "-r",
            "tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "-l",
            "chr19:6117965-6118004",
            "-o",
            "-",
            "--calling=thresholds"
        ])
    );

    Ok(())
}
