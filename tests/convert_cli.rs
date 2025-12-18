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
            "--no-ml", // disable ML for faster test
            "--region=chr19:6105700-6105800",
            "-o",
        ])
        .arg(&mpk)
        .silent()
        .succeeds()
        .wrap_err("Failed to run rastair call")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(temp_dir.path().join("test.bcf"))
        .silent()
        .succeeds()
        .wrap_err("Failed to convert to bcf")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(temp_dir.path().join("test.vcf.gz"))
        .silent()
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
            "--no-ml", // disable ML for faster test
            "--region=chr19:6105700-6105800",
            "-o",
        ])
        .arg(&mpk)
        .silent()
        .succeeds()
        .wrap_err("Failed to run rastair call")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(temp_dir.path().join("test.bed"))
        .silent()
        .succeeds()
        .wrap_err("Failed to convert to bed")?;

    Ok(())
}

#[test]
fn write_bcf_then_convert_to_bed() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let bcf = temp_dir.path().join("test.bcf");

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--no-ml", // disable ML for faster test
            "--region=chr19:6105700-6105800",
            "-o",
        ])
        .arg(&bcf)
        .silent()
        .succeeds()
        .wrap_err("Failed to run rastair call")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&bcf)
        .arg("--output")
        .arg(temp_dir.path().join("test.bed"))
        .succeeds()
        .wrap_err("Failed to convert to bed")?;

    // use clio::ClioPath;
    // use rastair_types::Probability;
    // rastair::convert(&rastair::ConvertParams {
    //     input: ClioPath::new(&bcf)?,
    //     input_format: None,
    //     output: ClioPath::new(temp_dir.path().join("test.bed"))?,
    //     output_format: None,
    //     bed_params: Default::default(),
    //     ml_threshold: Probability::new_panicky(0.8),
    // })?;

    Ok(())
}

#[test]
fn same_bed_via_mkv_and_bcf_conversion() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let mpk = temp_dir.path().join("test.mpk.lz4");
    let vcf = temp_dir.path().join("test.vcf");
    let bed = temp_dir.path().join("test.bed");
    let bed_from_mpk = temp_dir.path().join("from_mpk.bed");
    let bed_from_vcf = temp_dir.path().join("from_vcf.bed");

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--no-ml", // disable ML for faster test
            "--region=chr19:6105700-6105800",
            "-o",
        ])
        .arg(&mpk)
        .silent()
        .succeeds()
        .wrap_err("Failed to run rastair call to mpk")?;

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--no-ml", // disable ML for faster test
            "--region=chr19:6105700-6105800",
            "-o",
        ])
        .arg(&vcf)
        .silent()
        .succeeds()
        .wrap_err("Failed to run rastair call to vcf")?;

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--no-ml", // disable ML for faster test
            "--region=chr19:6105700-6105800",
            "--bed",
        ])
        .arg(&bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to run rastair call to bed")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(&bed_from_mpk)
        .silent()
        .succeeds()
        .wrap_err("Failed to convert mpk to bed")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&vcf)
        .arg("--output")
        .arg(&bed_from_vcf)
        .silent()
        .succeeds()
        .wrap_err("Failed to convert vcf to bed")?;

    // Compare the bed files
    let bed_directly = std::fs::read_to_string(&bed)?;

    let bed_vcf = std::fs::read_to_string(&bed_from_vcf)?;
    assert_eq!(bed_directly, bed_vcf, "BED files from VCF conversions do not match");

    let bed_mpk = std::fs::read_to_string(&bed_from_mpk)?;
    assert_eq!(bed_directly, bed_mpk, "BED files from MPK conversions do not match");

    Ok(())
}

#[test]
fn can_pipe_through() -> Result<()> {
    apply_common_filters!();

    let mut cmd = Command::new(insta_cmd::get_cargo_bin("/bin/bash"));
    cmd.arg("-c");
    cmd.env("NO_COLOR", "1");
    cmd.arg("cargo run -q --release -- call --fasta-file=tests/data/test.fasta.gz tests/data/test.bam --no-ml --vcf | head -n1000 | cargo run -q --release -- convert -f bcf -F bed | head -n5");

    assert_cmd_snapshot!(cmd);

    Ok(())
}
