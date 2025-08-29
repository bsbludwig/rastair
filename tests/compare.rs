mod utils;
use utils::*;

#[test]
#[ignore = "Depends on bcftools which we don't have on CI"]
fn compare_rastair_1_and_2() -> Result<()> {
    apply_common_filters!();
    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.vcf");

    let call = rastair()
        .args([
            "call",
            "-r",
            "tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "-l",
            "chr19",
            "-o",
        ])
        .arg(&temp_file)
        .status()
        .wrap_err("running rastair")?;
    assert!(call.success());

    // bcftools query -i "REF=='C' && CPG==1" -f "%CHROM\\t%POS0\\t%POS\\t%REF\\t%ALT\\t%AS_SB[\\t%M5mC\\t%GT\\t%DP]\n" tmp/test.bcf
    let rastair2 = Command::new("bcftools")
        .args([
            "query",
            "-i",
            "CPG==1",
            "-f",
            "%CHROM\t%POS0\t%POS\t%REF\t%ALT\t%AS_SB[\t%M5mC\t%GT\t%DP]\n",
        ])
        .arg(&temp_file)
        .output()
        .wrap_err("running bcftools query")?
        .stdout;

    let rastair1 = {
        let text =
            std::fs::read_to_string("./tests/data/rastair1.vcf").wrap_err("read rastair 1 vcf")?;
        text.lines()
            .filter(|line| !line.starts_with("#"))
            .filter_map(|line| line.split("\t").nth(1))
            .filter_map(|x| x.parse::<u32>().ok())
            .collect::<BTreeSet<u32>>()
    };
    let rastair2 = {
        // let text = std::fs::read_to_string(&temp_file).wrap_err("read rastair 2 vcf")?;
        // let text =
        //     std::fs::read_to_string("./tests/data/rastair2.vcf").wrap_err("read rastair 2 vcf")?;
        let text = String::from_utf8(rastair2).wrap_err("convert rastair 2 output to string")?;

        text.lines()
            .filter(|line| !line.starts_with("#"))
            .filter_map(|line| line.split("\t").nth(1))
            .filter_map(|x| x.parse::<u32>().ok())
            .collect::<BTreeSet<u32>>()
    };

    assert_debug_snapshot!(
        "in_1_but_not_2",
        rastair1.difference(&rastair2).map(|x| x + 1).collect::<Vec<_>>()
    );
    assert_debug_snapshot!(
        "in_2_but_not_1",
        rastair2.difference(&rastair1).map(|x| x + 1).collect::<Vec<_>>()
    );

    Ok(())
}
