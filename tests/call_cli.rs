use assert_cmd::Command;
use color_eyre::Result;
use predicates::prelude::*;

#[test]
fn random_call() -> Result<()> {
    let mut cmd = Command::cargo_bin("rastair2")?;
    cmd.args([
        "call",
        "-r",
        "tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "-l",
        "chr19:6105700-6105800",
        "-o",
        "-",
    ]);
    cmd.assert().success();

    Ok(())
}

#[test]
fn missing_bam() -> Result<()> {
    let mut cmd = Command::cargo_bin("rastair2")?;

    cmd.arg("call");
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
    cmd.arg("/path/to/nonexistent/file.bam");
    cmd.assert().failure().stderr(predicate::str::contains("No such file or directory"));

    Ok(())
}

#[test]
fn missing_fasta() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair2")?;

    cmd.arg("call");
    cmd.args(["--fasta-file", "tests/data/test_which_doesnt_exist.fasta"]);
    cmd.arg("tests/data/test.bam");
    cmd.assert().failure().stderr(predicate::str::contains("No such file or directory"));

    Ok(())
}

#[test]
fn validates_region_arg() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair2")?;

    cmd.args([
        "call",
        "-r",
        "tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "-l",
        "chr19:6105700-xxx",
    ]);
    cmd.assert().failure();

    insta::assert_snapshot!(str(cmd.output()?.stderr), @r"
    error: invalid value 'chr19:6105700-xxx' for '--region <REGION>': Invalid region string:
    chr19:6105700-xxx
                 ^


    For more information, try '--help'.
    ");

    Ok(())
}

#[allow(clippy::unwrap_used)]
fn str(vec: Vec<u8>) -> String {
    String::from_utf8(vec).unwrap()
}
