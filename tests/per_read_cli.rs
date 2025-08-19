mod utils;
use utils::*;

#[test]
fn simple_per_read_call() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.bed");

    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "--region=chr19:6105700-6105800",
        "--bed",
    ]).arg(&temp_file), @r#"
    success: true
    exit_code: 0
    ----- stdout -----

    ----- stderr -----
    [TIME] INFO rastair2::call: Wrote BED output file=[PATH]"
    [TIME] INFO rastair2: Call finished [DURATION]
    "#);

    Ok(())
}

#[test]
fn bed_file_is_sorted() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.bed");

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--region=chr19",
            "--bed",
        ])
        .arg(&temp_file)
        .status()?
        .is_success()
        .wrap_err("Failed to run rastair call")?;

    let bed_content =
        std::fs::read_to_string(&temp_file).wrap_err("Failed to read BED output file")?;

    let sorted_file =
        Command::new("sort").arg(&temp_file).output().wrap_err("Failed to sort BED file")?;
    let sorted_content = String::from_utf8(sorted_file.stdout)
        .wrap_err("Failed to convert sorted output to string")?;

    assert_eq!(bed_content, sorted_content, "BED file is not sorted");

    Ok(())
}
