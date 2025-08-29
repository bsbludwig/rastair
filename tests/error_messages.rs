mod utils;
use utils::*;

#[test]
fn missing_bam() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=test_data/test.fasta.gz",
        "/path/to/nonexistent/file.bam",
        "--vcf"
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
        "--fasta-file=tests/data/test_which_doesnt_exist.fasta",
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
fn missing_output_choice() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
    ]), @r"
    success: false
    exit_code: 1
    ----- stdout -----

    ----- stderr -----
    Error: 
       0: [91mNo output specified. Please specify at least one of `--vcf[=<PATH>]` or `--bed[=<PATH>]`.[0m

    Location:
       [35msrc/call.rs[0m:[35m90[0m

    Backtrace omitted. Run with RUST_BACKTRACE=1 environment variable to display it.
    Run with RUST_BACKTRACE=full to include source snippets.
    ");

    Ok(())
}

#[test]
fn validates_region_arg() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "--region=chr19:6105700-xxx",
        "--vcf",
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

#[test]
fn different_paths_for_bed_and_vcf() -> Result<()> {
    apply_common_filters!();
    assert_cmd_snapshot!(rastair().args([
        "call",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "--region=chr19:6105700",
        "--vcf=foo",
        "--bed=foo",
    ]), @r"
    success: false
    exit_code: 1
    ----- stdout -----

    ----- stderr -----
    Error: 
       0: [91mCan't write both VCF and BED output to the same file. Please specify different output files.[0m

    Location:
       [35msrc/call.rs[0m:[35m94[0m

    Backtrace omitted. Run with RUST_BACKTRACE=1 environment variable to display it.
    Run with RUST_BACKTRACE=full to include source snippets.
    ");

    Ok(())
}
