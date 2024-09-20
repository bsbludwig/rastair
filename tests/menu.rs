use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn unrecognised_command() -> Result<(), Box<dyn std::error::Error>>
{
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("something");
    cmd.assert()
       .failure()
       .stderr(predicate::str::contains("unrecognized subcommand"));

    Ok(())
}

#[test]
fn help_success() -> Result<(), Box<dyn std::error::Error>>
{
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("--help");
    cmd.assert().success();

    Ok(())
}

#[test]
fn version_success() -> Result<(), Box<dyn std::error::Error>>
{
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("--version");
    cmd.assert()
       .success()
       .stdout(predicate::str::contains("rastair"));
    Ok(())
}
