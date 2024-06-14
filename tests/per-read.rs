use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn missing_bam() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
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
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 18747); // Checked against methyldackel
    // Check header row is there
    let first_line = output_str.lines().next().unwrap();
    assert!(predicate::str::contains("#chr").eval(first_line));

    Ok(())
}

#[test]
fn uncompressed_fasta() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 18747); // Checked against methyldackel
    // Check header row is there
    let first_line = output_str.lines().next().unwrap();
    assert!(predicate::str::contains("#chr").eval(first_line));

    Ok(())
}

#[test]
fn threaded() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
    cmd.args(["-@", "2"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 18747); // Checked against methyldackel
    // Check header row is there
    let first_line = output_str.lines().next().unwrap();
    assert!(predicate::str::contains("#chr").eval(first_line));

    Ok(())
}

#[test]
fn report_all() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
        cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
        cmd.arg("-A");
        cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 24743); // Checked against methyldackel
    let total: u32 = output_str.lines()
        .map(|elem| {let row_elems: Vec<&str> = elem.split("\t").collect(); row_elems[8].parse::<u32>().unwrap_or_default()})
        .sum();
    assert_eq!(total, 65393);
    Ok(())
}

#[test]
fn report_all_threaded() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
        cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
        cmd.arg("-A");
        cmd.args(["-@", "2"]);
        cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 24743); // Checked against methyldackel
    let total: u32 = output_str.lines()
        .map(|elem| {let row_elems: Vec<&str> = elem.split("\t").collect(); row_elems[8].parse::<u32>().unwrap_or_default()})
        .sum();
    assert_eq!(total, 65393);
    Ok(())
}

#[test]
fn restrict_to_chromosome() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
    cmd.args(["-l", "bacteriophage_lambda_CpG"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 9646); // Checked against methyldackel
    let total: u32 = output_str.lines()
        .map(|elem| {let row_elems: Vec<&str> = elem.split("\t").collect(); row_elems[8].parse::<u32>().unwrap_or_default()})
        .sum();
    assert_eq!(total, 35655);

    let total_mod: usize = output_str.lines()
        .map(|elem| {let row_elems: Vec<&str> = elem.split("\t").collect(); row_elems[10].split(",").collect::<Vec<&str>>().len()})
        .sum();
    assert_eq!(total_mod, 34468);
    Ok(())
}

#[test]
fn restrict_to_chromosome_threaded() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
    cmd.args(["-l", "bacteriophage_lambda_CpG"]);
    cmd.args(["-@", "2"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 9646); // Checked against methyldackel
    let total: u32 = output_str.lines()
        .map(|elem| {let row_elems: Vec<&str> = elem.split("\t").collect(); row_elems[8].parse::<u32>().unwrap_or_default()})
        .sum();
    assert_eq!(total, 35655);

    let total_mod: usize = output_str.lines()
        .map(|elem| {let row_elems: Vec<&str> = elem.split("\t").collect(); row_elems[10].split(",").collect::<Vec<&str>>().len()})
        .sum();
    assert_eq!(total_mod, 34468);
    Ok(())
}

#[test]
fn restrict_to_region() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
    cmd.args(["-l", "bacteriophage_lambda_CpG:1-1000"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 121); // Checked against methyldackel

    Ok(())
}


#[test]
fn filter_mq_0() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
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
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
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
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
    cmd.args(["-l", "bacteriophage_lambda_CpG:6000-7000"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    let roi = output_str
        .lines()
        .filter(|l| predicate::str::contains("NB502094:69:HN2H2BGX5:2:22305:19777:14996").eval(l))
        .filter(|l| predicate::str::contains("6359").eval(l))
        .last()
        .unwrap_or_default();

    let elems = roi.split("\t").collect::<Vec<&str>>();
    assert_eq!(elems[10], "22,40,54,57,64");
    Ok(())
}

#[test]
fn correct_data_in_random_region_2() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
    cmd.args(["-l", "bacteriophage_lambda_CpG:6610-6738"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    let some_row = output_str
        .lines()
        .filter(|l| predicate::str::contains("NB502094:69:HN2H2BGX5:4:11501:13313:12598").eval(l))
        .last()
        .unwrap_or_default();
    let elems = some_row.split("\t").collect::<Vec<&str>>();
    assert_eq!(elems[8], "8");
    assert_eq!(elems[9], "0");
    assert_eq!(elems[10], ""); // these will all look like unmod, with no mod
    assert_eq!(elems[11], "8,14,32,38,43,55,63,77"); // these will all look like unmod, with no mod
    Ok(())
}

#[test]
fn reports_snps_at_cpg() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
    cmd.args(["-l", "bacteriophage_lambda_CpG:4962-5006"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    let some_row = output_str
    .lines()
    .filter(|l| predicate::str::contains("NB502094:69:HN2H2BGX5:3:22406:24759:2699").eval(l))
    .last()
    .unwrap_or_default();
let elems = some_row.split("\t").collect::<Vec<&str>>();
assert_eq!(elems[8], "9");
assert_eq!(elems[9], "8");
assert_eq!(elems[10], "4,11,14,24,31,35,52,76");
assert_eq!(elems[11], "");
assert_eq!(elems[12], "49"); // one C>A SNP here
    Ok(())
}

#[test]
fn does_not_report_deleted_positions() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
    cmd.args(["-l", "bacteriophage_lambda_CpG:114-158"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
        .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    let some_row = output_str
    .lines()
    .filter(|l| predicate::str::contains("NB502094:69:HN2H2BGX5:3:22506:13203:14643").eval(l))
    .last()
    .unwrap_or_default();
let elems = some_row.split("\t").collect::<Vec<&str>>();
assert_eq!(elems[8], "2");
assert_eq!(elems[9], "2");
assert_eq!(elems[10], "26,61");
    Ok(())
}