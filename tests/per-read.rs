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