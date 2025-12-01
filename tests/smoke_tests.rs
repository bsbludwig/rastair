mod utils;
use utils::*;

const CALL_TEST_BAM: [&str; 3] =
    ["call", "--fasta-file=tests/data/test.fasta.gz", "tests/data/test.bam"];
const CHR19: &str = "--region=chr19";
const NO_ML: &str = "--thresholds"; // disable ML for faster tests

#[test]
fn count_bed_lines() -> Result<()> {
    apply_common_filters!();

    let call = rastair().args(CALL_TEST_BAM).args([CHR19, NO_ML, "-c"]).output()?;

    assert_snapshot!(call.stderr(), @r#"
    [TIME] INFO rastair::call::writer: Wrote BED output file="-"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    let stdout = call.stdout();
    assert_snapshot!(stdout.lines().count(), @"1790");

    Ok(())
}

#[test]
fn count_vcf_lines() -> Result<()> {
    apply_common_filters!();

    let call = rastair().args(CALL_TEST_BAM).args([CHR19, NO_ML, "-c", "--vcf"]).output()?;

    assert_snapshot!(call.stderr(), @r#"
    [TIME] INFO rastair::call::writer: Wrote VCF output file="-"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    let stdout = call.stdout();
    assert_snapshot!(stdout.lines().count(), @"1497");

    Ok(())
}
