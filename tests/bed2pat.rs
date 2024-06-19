use assert_cmd::prelude::*;
// Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs
use tempfile::NamedTempFile;

fn stage_read_bed(region: &str) -> Result<NamedTempFile, Box<dyn std::error::Error>>
{
    // Run rastair per-read on some region
    let mut cmd = Command::cargo_bin("rastair")?;
    let file = NamedTempFile::new()?;
    let write_handle = file.reopen()?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
    if region.len() > 0
    {
        cmd.args(["-l", region]);
    }
    cmd.arg("test_data/test.bam")
       .stdout(write_handle)
       .status()
       .expect("Failed to create per-read file needed for tests");

    Ok(file)
}

#[test]
fn missing_bed() -> Result<(), Box<dyn std::error::Error>> {
        // stage input file
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("bed2pat");
    cmd.args(["test_data/test.fasta.gz"]);
    cmd.arg("/path/to/nonexistent/file.bed");
    cmd.assert()
       .failure()
       .stderr(predicate::str::contains("No such file or directory"));

    Ok(())
}

#[test]
fn missing_fasta() -> Result<(), Box<dyn std::error::Error>> {
    let file = stage_read_bed("bacteriophage_lambda_CpG:1-1000")?;

    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("bed2pat");
    cmd.args(["test_data/test_which_doesnt_exist.fasta"]);
    cmd.arg(file.path());
    cmd.assert()
       .failure()
       .stderr(predicate::str::contains("No such file or directory"));

    Ok(())
}

#[test]
fn can_create_pat() -> Result<(), Box<dyn std::error::Error>> {
    let file = stage_read_bed("bacteriophage_lambda_CpG:1-1000")?;
    println!("Size of {}: {}", file.path().to_str().unwrap_or_default(), file.as_ref().metadata().unwrap().len());

    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("bed2pat");
    cmd.args(["test_data/test.fasta.gz"]);
    cmd.arg(file.path());
    cmd.assert()
       .success();

    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    let line_count = output_str
        .lines()
        .count();
    assert!(line_count > 0);
    Ok(())
}