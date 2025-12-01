mod utils;
use rustc_hash::FxHashSet;
use utils::*;

const CALL_TEST_BAM: [&str; 3] =
    ["call", "--fasta-file=tests/data/test.fasta.gz", "tests/data/test.bam"];
const CHR19: &str = "--region=chr19";
const NO_ML: &str = "--no-ml"; // disable ML for faster tests

#[test]
fn count_bed_lines() -> Result<()> {
    apply_common_filters!();

    let call = rastair().args(CALL_TEST_BAM).args([CHR19, NO_ML, "-c"]).output()?;

    assert_snapshot!(call.stderr(), @r#"
    [TIME] INFO rastair::call::writer: Wrote BED output file="-"
    [TIME] INFO rastair: Call finished [DURATION]
    "#);

    let stdout = call.stdout();
    let cpgs = stdout
        .lines()
        .filter(|l| !l.starts_with('#'))
        .filter(|l| l.trim().ends_with("REF"))
        .count();
    let denovos = stdout
        .lines()
        .filter(|l| !l.starts_with('#'))
        .filter(|l| l.trim().ends_with("NEW"))
        .count();
    assert_snapshot!(cpgs, @"1448");
    assert_snapshot!(denovos, @"341");

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
    let cpgs = stdout
        .lines()
        .filter(|l| !l.starts_with('#'))
        .filter(|l| l.contains("CPG") && !l.contains("CPGnovo"))
        .count();
    let denovos =
        stdout.lines().filter(|l| !l.starts_with('#')).filter(|l| l.contains("CPGnovo")).count();

    assert_snapshot!(cpgs, @"1444");
    assert_snapshot!(denovos, @"6");

    Ok(())
}

#[test]
fn every_bed_line_should_be_in_vcf() -> Result<()> {
    apply_common_filters!();

    let mut bed_call = rastair().args(CALL_TEST_BAM).args([CHR19, NO_ML, "-c"]).output()?;
    bed_call.succeeds()?;

    // bed line are like `chr19	start	end … REF|NEW`
    // we already filter by chromosome in call, so let's collect all start positions
    let (bed_cpg, bed_de_novo): (FxHashSet<u32>, FxHashSet<u32>) =
        bed_call.stdout().lines().filter(|line| !line.starts_with('#')).fold(
            (FxHashSet::default(), FxHashSet::default()),
            |(mut cpg_set, mut denovo_set), line| {
                let fields: Vec<&str> = line.trim().split('\t').collect();
                let pos: u32 = fields[2].parse().unwrap(); // take end position since it's 1-based
                let de_novo: bool = line.ends_with("NEW");
                if de_novo {
                    denovo_set.insert(pos);
                } else {
                    cpg_set.insert(pos);
                }
                (cpg_set, denovo_set)
            },
        );

    let mut vcf_call =
        rastair().args(CALL_TEST_BAM).args([CHR19, NO_ML, "-c", "--vcf"]).output()?;
    vcf_call.succeeds()?;

    let (vcf_cpg, vcf_de_novo): (FxHashSet<u32>, FxHashSet<u32>) =
        vcf_call.stdout().lines().filter(|line| !line.starts_with('#')).fold(
            (FxHashSet::default(), FxHashSet::default()),
            |(mut cpg_set, mut denovo_set), line| {
                let fields: Vec<&str> = line.trim().split('\t').collect();
                let pos: u32 = fields[1].parse().unwrap();
                let info_field = fields[7];
                let is_cpg = info_field.ends_with("CPG");
                let is_denovo = info_field.ends_with("CPGnovo");
                if is_cpg {
                    cpg_set.insert(pos);
                }
                if is_denovo {
                    denovo_set.insert(pos);
                }
                (cpg_set, denovo_set)
            },
        );

    assert!(
        vcf_cpg.is_subset(&bed_cpg),
        "Not all VCF CpG calls are in BED, e.g. {:?}",
        vcf_cpg.difference(&bed_cpg).take(5).collect::<Vec<&u32>>()
    );
    assert!(
        bed_cpg.is_subset(&vcf_cpg),
        "Not all BED CpG calls are in VCF, e.g. {:?}",
        bed_cpg.difference(&vcf_cpg).take(5).collect::<Vec<&u32>>()
    );
    assert!(
        vcf_de_novo.is_subset(&bed_de_novo),
        "Not all VCF de-novo CpG calls are in BED, e.g. {:?}",
        vcf_de_novo.difference(&bed_de_novo).take(5).collect::<Vec<&u32>>()
    );
    assert!(
        bed_de_novo.is_subset(&vcf_de_novo),
        "Not all BED de-novo CpG calls are in VCF, e.g {:?}",
        bed_de_novo.difference(&vcf_de_novo).take(5).collect::<Vec<&u32>>()
    );

    Ok(())
}
