use assert_cmd::Command;
use predicates::prelude::*;

#[test]
fn random_call() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("rastair2")?;
    cmd.args([
        "call",
        "-r",
        "tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "-l",
        "chr19:6105700-6105800",
    ]);
    cmd.assert().success();

    Ok(())
}

#[test]
fn missing_bam() -> anyhow::Result<()> {
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
