use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn missing_bam() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("call");
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

    cmd.arg("call");
        cmd.args(["--fasta-file", "test_data/test_which_doesnt_exist.fasta"]);
        cmd.arg("test_data/test.bam");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("No such file or directory"));

    Ok(())
}

#[test]
fn refuse_0_threads() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("call");
        cmd.args(["--fasta-file", "test_data/test.fasta"]);
        cmd.args(["--threads", "0"]);
        cmd.arg("test_data/test.bam");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("invalid value"));

    Ok(())
}

#[test]
fn refuse_0_chunks() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("call");
        cmd.args(["--fasta-file", "test_data/test.fasta"]);
        cmd.args(["--chunk-size", "0"]);
        cmd.arg("test_data/test.bam");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("invalid value"));

    Ok(())
}

#[test]
fn default_settings() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("call");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
            .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 7898);
    // Check header row is there
    let first_line = output_str.lines().next().unwrap();
    assert!(predicate::str::contains("#chr").eval(first_line));
    Ok(())
}

#[test]
fn finds_right_number_of_positions_with_threads() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("call");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.args(["-l", "bacteriophage_lambda_CpG"]);
    cmd.args(["-@", "2"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
            .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 6110); // these values are empirical, ie more for detection of regressions. I checked them manually against MethylDackel and via IGV and they look right

    let mut total = 0;
    let mut total_plus = 0;
    let mut total_minus = 0;
    for line in output_str.lines()
    {
        let elems : Vec<&str> = line.split_ascii_whitespace().collect();
        if elems[1] == "start"
        {
            // skip header col
            continue;
        }
        // sum total number of mod/unmod read pos
        total = total + elems[6].parse::<i32>().unwrap();
        total = total + elems[7].parse::<i32>().unwrap();
        match elems[5]
        {
            "+" => total_plus += 1,
            "-" => total_minus += 1,
            _   => ()
        };
    }
    assert_eq!(total, 35201); // these values are empirical, ie more for detection of regressions
    assert_eq!(total_plus+1, total_minus); // there's one - strand pos that's not covered in this test dataset
    Ok(())
}

#[test]
fn finds_right_number_of_positions_with_chunks() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("call");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.args(["-s", "1000"]);
    cmd.args(["-l", "bacteriophage_lambda_CpG"]);
    cmd.args(["-@", "2"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
            .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 6110); // these values are empirical, ie more for detection of regressions. I checked them manually against MethylDackel and via IGV and they look right

    let mut total = 0;
    for line in output_str.lines()
    {
        let elems : Vec<&str> = line.split_ascii_whitespace().collect();
        if elems[1] == "start"
        {
            // skip header col
            continue;
        }
        // sum total number of mod/unmod read pos
        total = total + elems[6].parse::<i32>().unwrap();
        total = total + elems[7].parse::<i32>().unwrap();
    }
    assert_eq!(total, 35201); // these values are empirical, ie more for detection of regressions
    Ok(())
}

#[test]
fn finds_right_number_of_positions_with_trimming() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("call");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.args(["--nOT", "5,5,5,5"]);
    cmd.args(["--nOB", "5,5,5,5"]);
    cmd.args(["-l", "bacteriophage_lambda_CpG"]);
    cmd.args(["-@", "2"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
            .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 6083);
    let mut total = 0;
    for line in output_str.lines()
    {
        let elems : Vec<&str> = line.split_ascii_whitespace().collect();
        if elems[1] == "start"
        {
            // skip header col
            continue;
        }
        // sum total number of mod/unmod read pos
        total = total + elems[6].parse::<i32>().unwrap();
        total = total + elems[7].parse::<i32>().unwrap();
    }
    assert_eq!(total, 31604);
    Ok(())
}

#[test]
fn allow_restriction_to_chromosome() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("call");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.args(["--region", "bacteriophage_lambda_CpG"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
            .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 6110);

    Ok(())
}

#[test]
fn allow_restriction_to_region() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("call");
    cmd.args(["--fasta-file", "test_data/test.fasta"]);
    cmd.args(["--region", "bacteriophage_lambda_CpG:1-1000"]);
    cmd.arg("test_data/test.bam");
    cmd.assert()
            .success();
    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output_str.lines().count(), 143);

    Ok(())
}