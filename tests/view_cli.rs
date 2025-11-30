mod utils;
use utils::*;

#[test]
fn write_mpk_then_view_stdout() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let mpk = temp_dir.path().join("test.mpk.lz4");

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--no-ml", // disable ML for faster test
            "--region=chr19:6105700-6105750",
            "-o",
        ])
        .arg(&mpk)
        .succeeds()
        .wrap_err("Failed to run rastair call")?;

    rastair().arg("view").arg(&mpk).succeeds().wrap_err("Failed to view mpk file")?;

    Ok(())
}

#[test]
fn write_mpk_then_view_file() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let mpk = temp_dir.path().join("test.mpk.lz4");
    let json = temp_dir.path().join("test.json");

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--no-ml", // disable ML for faster test
            "--region=chr19:6105700-6105750",
            "--vcf",
        ])
        .arg(&mpk)
        .succeeds()
        .wrap_err("Failed to run rastair call")?;

    rastair()
        .arg("view")
        .arg(&mpk)
        .arg("--output")
        .arg(&json)
        .succeeds()
        .wrap_err("Failed to view mpk file")?;

    let json_content =
        std::fs::read_to_string(&json).wrap_err("Failed to read JSON output file")?;
    let lines = json_content.lines().count();

    let json_objects =
        serde_json::Deserializer::from_str(&json_content).into_iter::<serde_json::Value>().count();

    assert_eq!(lines, json_objects);

    Ok(())
}
