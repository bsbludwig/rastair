mod utils;
use utils::*;

#[test]
fn write_mpk_then_convert_to_bcf() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let mpk = temp_dir.path().join("test.mpk.lz4");

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--thresholds",
            "-o",
        ])
        .arg(&mpk)
        .succeeds()
        .wrap_err("Failed to run rastair call")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(temp_dir.path().join("test.bcf"))
        .succeeds()
        .wrap_err("Failed to convert to bcf")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(temp_dir.path().join("test.vcf.gz"))
        .succeeds()
        .wrap_err("Failed to convert to vcf.gz")?;

    Ok(())
}

#[test]
fn write_mpk_then_convert_to_bed() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let mpk = temp_dir.path().join("test.mpk.lz4");

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--thresholds",
            "-o",
        ])
        .arg(&mpk)
        .succeeds()
        .wrap_err("Failed to run rastair call")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(temp_dir.path().join("test.bed"))
        .succeeds()
        .wrap_err("Failed to convert to bed")?;

    Ok(())
}

#[test]
fn write_bcf_then_convert_to_bed() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let mpk = temp_dir.path().join("test.bcf");

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--thresholds",
            "-o",
        ])
        .arg(&mpk)
        .succeeds()
        .wrap_err("Failed to run rastair call")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(temp_dir.path().join("test.bed"))
        .succeeds()
        .wrap_err("Failed to convert to bed")?;

    Ok(())
}

#[test]
fn can_pipe_through() -> Result<()> {
    apply_common_filters!();

    let mut cmd = Command::new(insta_cmd::get_cargo_bin("/bin/bash"));
    cmd.arg("-c");
    cmd.env("NO_COLOR", "1");
    cmd.arg("cargo run -q --release -- call --fasta-file=tests/data/test.fasta.gz tests/data/test.bam --thresholds --vcf | head -n1000 | cargo run -q --release -- convert -f bcf -F bed | head -n5");

    assert_cmd_snapshot!(cmd);

    Ok(())
}
