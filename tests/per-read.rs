use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn missing_bam() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.arg("/path/to/nonexistent/file.bam");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("No such file or directory"));

    Ok(())
}

#[test]
fn missing_fasta() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test_which_doesnt_exist.fasta"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("No such file or directory"));

    Ok(())
}

#[test]
fn default_settings() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 6601); // Checked against methyldackel
    // Check header row is there
    let first_line = output_str.lines().next().unwrap();
    assert!(predicate::str::contains("#chr").eval(first_line));

    Ok(())
}

#[test]
fn report_all() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
        cmd.args(["--fasta-file", "test_data/test.fasta"]);
        cmd.arg("-A");
        cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 6659); // Checked against methyldackel
    let total: u32 = output_str.lines()
        .map(|elem| {let row_elems: Vec<&str> = elem.split("\t").collect(); row_elems[8].parse::<u32>().unwrap_or_default()})
        .sum();
    assert_eq!(total, 53493);
    Ok(())
}

#[test]
fn restrict_to_chromosome() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.args(["-l", "bacteriophage_lambda_CpG"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 5213); // Checked against methyldackel
    let total: u32 = output_str.lines()
        .map(|elem| {let row_elems: Vec<&str> = elem.split("\t").collect(); row_elems[8].parse::<u32>().unwrap_or_default()})
        .sum();
    assert_eq!(total, 35339);

    let total_mod: usize = output_str.lines()
        .map(|elem| {let row_elems: Vec<&str> = elem.split("\t").collect(); row_elems[10].split(",").collect::<Vec<&str>>().len()})
        .sum();
    assert_eq!(total_mod, 33910);
    Ok(())
}

#[test]
fn restrict_to_region() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.args(["-l", "bacteriophage_lambda_CpG:1-1000"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 47); // Checked against methyldackel

    Ok(())
}