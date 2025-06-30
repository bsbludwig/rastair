use std::collections::BTreeSet;

mod utils;
use color_eyre::eyre::eyre;
use utils::*;

#[test]
fn random_call() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.vcf");

    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "--region=chr19:6105700-6105800",
        "--calling=thresholds",
        "-o",
    ]).arg(&temp_file), @r#"
    success: true
    exit_code: 0
    ----- stdout -----

    ----- stderr -----
    [TIME] INFO rastair2::call: Wrote output file=[PATH]"
    [TIME] INFO rastair2: Call finished [DURATION]
    "#);

    Ok(())
}

#[test]
fn pipe_to_stdout() -> Result<()> {
    apply_common_filters!();

    // Please manually verify that only the VCF output is printed to stdout and
    // logging goes to stderr.
    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "--region=chr19:6105700-6105750",
    ]));

    Ok(())
}

#[test]
fn includes_all_cpgs_when_methylation_calling() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(
        "does not include unmethylated cpgs without alts",
        rastair().args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--region=chr19:6117965-6118004",
        ])
    );

    assert_cmd_snapshot!(
        "includes all cpgs",
        rastair().args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--region=chr19:6117965-6118004",
            "--calling=thresholds"
        ])
    );

    Ok(())
}

#[test]
fn includes_only_cpgs_when_methylation_calling() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(
        "includes only cpgs",
        rastair().args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--region=chr19:6117965-6118004",
            "--calling=thresholds",
            "--cpgs-only"
        ])
    );

    Ok(())
}

#[test]
fn segmentation_overlaps_do_not_cause_duplicate_records() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.vcf");

    let call = rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--region=chr19",
            "--segment-max-length=10000",
            "--segment-overlap=300",
            "--calling=thresholds",
            "-o",
        ])
        .arg(&temp_file)
        .status()
        .wrap_err("running rastair2")?;
    assert!(call.success());

    let text = std::fs::read_to_string(&temp_file).wrap_err("read rastair 2 vcf")?;
    text.lines()
        .filter(|line| !line.starts_with("#"))
        .filter_map(|line| line.split("\t").nth(1))
        .filter_map(|x| x.parse::<u32>().ok())
        .try_fold(BTreeSet::new(), |mut set, position| {
            let is_new = set.insert(position);
            if is_new { Ok(set) } else { Err(eyre!("Duplicate position found: {}", position)) }
        })?;

    Ok(())
}
