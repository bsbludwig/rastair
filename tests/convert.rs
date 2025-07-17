mod utils;
use utils::*;

#[test]
fn write_mpk_then_convert_to_bcf() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let mpk = temp_dir.path().join("test.mpk.lz4");

    rastair()
        .args(["call", "--fasta-file=tests/data/test.fasta.gz", "tests/data/test.bam", "-o"])
        .arg(&mpk)
        .status()?
        .is_success()
        .wrap_err("Failed to run rastair call")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(temp_dir.path().join("test.bcf"))
        .status()?
        .is_success()
        .wrap_err("Failed to convert to bcf")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(temp_dir.path().join("test.vcf.gz"))
        .status()?
        .is_success()
        .wrap_err("Failed to convert to vcf.gz")?;

    Ok(())
}

#[test]
fn write_mpk_then_convert_to_bed() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let mpk = temp_dir.path().join("test.mpk.lz4");

    rastair()
        .args(["call", "--fasta-file=tests/data/test.fasta.gz", "tests/data/test.bam", "-o"])
        .arg(&mpk)
        .status()?
        .is_success()
        .wrap_err("Failed to run rastair call")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(temp_dir.path().join("test.bed"))
        .status()?
        .is_success()
        .wrap_err("Failed to convert to bed")?;

    Ok(())
}

#[test]
fn write_bcf_then_convert_to_bed() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let mpk = temp_dir.path().join("test.bcf");

    rastair()
        .args(["call", "--fasta-file=tests/data/test.fasta.gz", "tests/data/test.bam", "-o"])
        .arg(&mpk)
        .status()?
        .is_success()
        .wrap_err("Failed to run rastair call")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(temp_dir.path().join("test.bed"))
        .status()?
        .is_success()
        .wrap_err("Failed to convert to bed")?;

    Ok(())
}
