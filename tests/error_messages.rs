mod utils;
use utils::*;

#[test]
fn missing_bam() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file", "test_data/test.fasta.gz",
        "/path/to/nonexistent/file.bam"
    ]), @r#"
    success: false
    exit_code: 2
    ----- stdout -----

    ----- stderr -----
    error: Invalid value for --fasta-file <FASTA_FILE>: Invalid path "test_data/test.fasta.gz": No such file or directory (os error 2)

    Usage: rastair2 call [OPTIONS] --fasta-file <FASTA_FILE> <BAM_FILE>

    For more information, try '--help'.
    "#);

    Ok(())
}

#[test]
fn missing_fasta() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file", "tests/data/test_which_doesnt_exist.fasta",
        "tests/data/test.bam"
    ]), @r#"
    success: false
    exit_code: 2
    ----- stdout -----

    ----- stderr -----
    error: Invalid value for --fasta-file <FASTA_FILE>: Invalid path "tests/data/test_which_doesnt_exist.fasta": No such file or directory (os error 2)

    Usage: rastair2 call [OPTIONS] --fasta-file <FASTA_FILE> <BAM_FILE>

    For more information, try '--help'.
    "#);

    Ok(())
}

#[test]
fn validates_region_arg() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(rastair().args([
        "call",
        "-r",
        "tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "-l",
        "chr19:6105700-xxx",
    ]), @r"
    success: false
    exit_code: 2
    ----- stdout -----

    ----- stderr -----
    error: invalid value 'chr19:6105700-xxx' for '--region <REGION>': Invalid region string:
    chr19:6105700-xxx
                 ^


    For more information, try '--help'.
    ");

    Ok(())
}
