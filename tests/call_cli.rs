use color_eyre::{Result, eyre::Context as _};
use insta::assert_snapshot;
use insta_cmd::{assert_cmd_snapshot, get_cargo_bin};
use std::process::Command;
use tempfile::TempDir;

fn rastair() -> Command {
    let mut cmd = Command::new(get_cargo_bin("rastair2"));
    cmd.env("NO_COLOR", "1");
    cmd
}

macro_rules! apply_common_filters {
    {} => {
        let mut settings = insta::Settings::clone_current();
        settings.add_filter(r"\w{4}-[0-9T\-:.]+Z\s", "[TIME]");
        settings.add_filter(r"duration=[\w.]+", "[DURATION]");
        settings.add_filter(r": close time.*", " [CLOSE]");
        settings.add_filter(r"Wrote output to /.*/test.vcf", "Wrote output to [PATH]");
        let _bound = settings.bind_to_scope();
    }
}

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
        "-o",
    ]).arg(&temp_file), @r"
    success: true
    exit_code: 0
    ----- stdout -----

    ----- stderr -----
    [TIME] INFO rastair2::call: Processed 1 segments
    [TIME] INFO rastair2::call: Wrote output to [PATH]
    [TIME] INFO rastair2: Call finished [DURATION]
    ");

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
        "-o",
        "-"
    ]));

    Ok(())
}

#[test]
fn missing_bam() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file", "test_data/test.fasta.gz",
        "/path/to/nonexistent/file.bam"
    ]), @r#"
    success: false
    exit_code: 2
    ----- stdout -----

    ----- stderr -----
    error: Invalid value for --fasta-file <FASTA_FILE>: Invalid path "test_data/test.fasta.gz": No such file or directory (os error 2)

    Usage: rastair2 call [OPTIONS] --fasta-file <FASTA_FILE> --vcf-output <VCF_OUTPUT> <BAM_FILE>

    For more information, try '--help'.
    "#);

    Ok(())
}

#[test]
fn missing_fasta() -> Result<(), Box<dyn std::error::Error>> {
    apply_common_filters!();
    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file", "tests/data/test_which_doesnt_exist.fasta",
        "tests/data/test.bam"
    ]), @r#"
    success: false
    exit_code: 2
    ----- stdout -----

    ----- stderr -----
    error: Invalid value for --fasta-file <FASTA_FILE>: Invalid path "tests/data/test_which_doesnt_exist.fasta": No such file or directory (os error 2)

    Usage: rastair2 call [OPTIONS] --fasta-file <FASTA_FILE> --vcf-output <VCF_OUTPUT> <BAM_FILE>

    For more information, try '--help'.
    "#);

    Ok(())
}

#[test]
fn validates_region_arg() -> Result<(), Box<dyn std::error::Error>> {
    apply_common_filters!();
    assert_cmd_snapshot!(rastair().args([
        "call",
        "-r",
        "tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "-l",
        "chr19:6105700-xxx",
    ]), @r"
    success: false
    exit_code: 2
    ----- stdout -----

    ----- stderr -----
    error: invalid value 'chr19:6105700-xxx' for '--region <REGION>': Invalid region string:
    chr19:6105700-xxx
                 ^


    For more information, try '--help'.
    ");

    Ok(())
}
