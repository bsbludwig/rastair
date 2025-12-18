mod utils;
use utils::*;

#[test]
fn convert_from_mpk() -> Result<()> {
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
fn convert_vcf_to_same_bed() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let vcf = temp_dir.path().join("test.vcf");
    let bed = temp_dir.path().join("test.bed");
    let bed_from_vcf = temp_dir.path().join("from_vcf.bed");

    let args = &[
        "call",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "--region=chr19:6105700-6105800",
    ];

    rastair()
        .args(args)
        .arg("-o")
        .arg(&vcf)
        .silent()
        .succeeds()
        .wrap_err("Failed to run rastair call to vcf")?;

    rastair()
        .args(args)
        .arg("-c")
        .arg("--bed")
        .arg(&bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to run rastair call to vcf")?;

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

    Ok(())
}

#[test]
fn can_pipe_through() -> Result<()> {
    apply_common_filters!();

    let mut cmd = Command::new(insta_cmd::get_cargo_bin("/bin/bash"));
    cmd.arg("-c");
    cmd.env("NO_COLOR", "1");
    cmd.arg("cargo run -q -- call --fasta-file=tests/data/test.fasta.gz tests/data/test.bam --region=chr19:6105000-6105800 --no-ml --vcf | head -n1000 | cargo run -q -- convert -f bcf -F bed | head -n5");

    assert_cmd_snapshot!(cmd);

    Ok(())
}
