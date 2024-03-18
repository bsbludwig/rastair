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
fn threaded() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.args(["-@", "2"]);
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
fn report_all_threaded() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
        cmd.args(["--fasta-file", "test_data/test.fasta"]);
        cmd.arg("-A");
        cmd.args(["-@", "2"]);
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
fn restrict_to_chromosome_threaded() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.args(["-l", "bacteriophage_lambda_CpG"]);
    cmd.args(["-@", "2"]);
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


#[test]
fn filter_mq_0() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.args(["-q", "5"]);
    cmd.args(["-l", "chr19"]);
    cmd.args(["-@", "2"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    let roi = output_str
        .lines()
        .filter(|l|!predicate::str::contains("mapq").eval(l))
        .fold(0, |acc, line| {
            let elems = line.split_whitespace().collect::<Vec<&str>>();
            if elems[4].parse::<usize>().unwrap_or(255) < 5
            {
                acc + 1
            }
            else {
                acc
            }
        });
    assert_eq!(roi, 0);

    cmd = Command::cargo_bin("rastair")?;
    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.args(["-q", "0"]);
    cmd.args(["-@", "2"]);
    cmd.args(["-l", "chr19"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    let roi = output_str
        .lines()
        .filter(|l|!predicate::str::contains("mapq").eval(l))
        .fold(0, |acc, line| {
            let elems = line.split_whitespace().collect::<Vec<&str>>();
            if elems[4].parse::<usize>().unwrap_or(255) < 5
            {
                acc + 1
            }
            else {
                acc
            }
        });
    assert!(roi > 0);
    Ok(())
}

#[test]
fn correct_pos_with_skips() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.args(["-l", "bacteriophage_lambda_CpG:6000-7000"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    let roi = output_str
        .lines()
        .filter(|l| predicate::str::contains("A00711:92:HMH3WDSXX:3:2662:7328:10300").eval(l))
        .filter(|l| predicate::str::contains("6048").eval(l))
        .last()
        .unwrap_or_default();
    let elems = roi.split_ascii_whitespace().collect::<Vec<&str>>();
    assert_eq!(elems[10], "11,16,29,59,67,79,95,108,122");
    Ok(())
}