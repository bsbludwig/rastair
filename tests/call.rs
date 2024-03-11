use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn unrecognised_command() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("something");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("unrecognized subcommand"));

    Ok(())
}

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
		cmd.arg("test_data/file.bam");
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
fn has_header() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("call");
		cmd.args(["--fasta-file", "test_data/test.fasta"]);
		cmd.arg("test_data/test.bam");
    cmd.assert()
        .success()
				.stdout(predicate::str::contains("#chr\t"));
    Ok(())
}

#[test]
fn finds_right_number_of_positions() -> Result<(), Box<dyn std::error::Error>> {
	let mut cmd = Command::cargo_bin("rastair")?;

	cmd.arg("call");
	cmd.args(["--fasta-file", "test_data/test.fasta"]);
	cmd.arg("test_data/test.bam");
	cmd.assert()
			.success();
	let output = cmd.output().unwrap();
	let output_str = String::from_utf8_lossy(&output.stdout);
	assert_eq!(output_str.lines().count(), 6486);
	Ok(())
}

#[test]
fn finds_right_number_of_positions_with_threads() -> Result<(), Box<dyn std::error::Error>> {
	let mut cmd = Command::cargo_bin("rastair")?;

	cmd.arg("call");
	cmd.args(["--fasta-file", "test_data/test.fasta"]);
	cmd.args(["-@", "2"]);
	cmd.arg("test_data/test.bam");
	cmd.assert()
			.success();
	let output = cmd.output().unwrap();
	let output_str = String::from_utf8_lossy(&output.stdout);
	assert_eq!(output_str.lines().count(), 6486);
	Ok(())
}

#[test]
fn finds_right_number_of_positions_with_chunks() -> Result<(), Box<dyn std::error::Error>> {
	let mut cmd = Command::cargo_bin("rastair")?;

	cmd.arg("call");
	cmd.args(["--fasta-file", "test_data/test.fasta"]);
	cmd.args(["-s", "200"]);
	cmd.arg("test_data/test.bam");
	cmd.assert()
			.success();
	let output = cmd.output().unwrap();
	let output_str = String::from_utf8_lossy(&output.stdout);
	assert_eq!(output_str.lines().count(), 6486); // these values are empirical, ie more for detection of regressions. I checked them manually against MethylDackel and via IGV and they look right

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
	assert_eq!(total, 48121); // these values are empirical, ie more for detection of regressions
	Ok(())
}

#[test]
fn finds_right_number_of_positions_with_trimming() -> Result<(), Box<dyn std::error::Error>> {
	let mut cmd = Command::cargo_bin("rastair")?;

	cmd.arg("call");
	cmd.args(["--fasta-file", "test_data/test.fasta"]);
	cmd.args(["--nOT", "5,5,5,5"]);
	cmd.args(["--nOB", "5,5,5,5"]);
	cmd.arg("test_data/test.bam");
	cmd.assert()
			.success();
	let output = cmd.output().unwrap();
	let output_str = String::from_utf8_lossy(&output.stdout);
	assert_eq!(output_str.lines().count(), 6470);
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
	assert_eq!(total, 45811);
	Ok(())
}