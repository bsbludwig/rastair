mod utils;
use utils::*;

#[test]
fn convert_vcf_to_vcf_just_copies_the_file() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let vcf = temp_dir.path().join("test.vcf");
    let vcf2 = temp_dir.path().join("test2.vcf");

    rastair()
        .args(CALL_TEST_BAM)
        .args([NO_ML, CHR19_SMALL, "-o"])
        .arg(&vcf)
        .silent()
        .succeeds()
        .wrap_err("Failed to run rastair call")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&vcf)
        .arg("--output")
        .arg(&vcf2)
        .silent()
        .succeeds()
        .wrap_err("Failed to convert to vcf")?;

    // assert that the files are identical
    let vcf_content = std::fs::read(&vcf)?;
    let vcf2_content = std::fs::read(&vcf2)?;
    assert_eq!(vcf_content, vcf2_content, "VCF files are not identical after conversion");

    Ok(())
}

#[test]
fn convert_vcf_to_bcf() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let vcf = temp_dir.path().join("test.vcf");
    let bcf = temp_dir.path().join("test.bcf");

    rastair()
        .args(CALL_TEST_BAM)
        .args([NO_ML, CHR19_SMALL, "-o"])
        .arg(&vcf)
        .silent()
        .succeeds()
        .wrap_err("Failed to run rastair call")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&vcf)
        .arg("--output")
        .arg(&bcf)
        .succeeds()
        .wrap_err("Failed to convert to bcf")?;

    Ok(())
}

#[test]
fn convert_bcf_to_vcf() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let bcf = temp_dir.path().join("test.bcf");
    let vcf = temp_dir.path().join("test.vcf");

    rastair()
        .args(CALL_TEST_BAM)
        .args([NO_ML, CHR19_SMALL, "-o"])
        .arg(&bcf)
        .silent()
        .succeeds()
        .wrap_err("Failed to run rastair call")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&bcf)
        .arg("--output")
        .arg(&vcf)
        .succeeds()
        .wrap_err("Failed to convert to vcf")?;

    Ok(())
}

#[test]
fn convert_bcf_to_vcf_stdout() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let bcf = temp_dir.path().join("test.bcf");

    rastair()
        .args(CALL_TEST_BAM)
        .args([NO_ML, CHR19_SMALL, "-o"])
        .arg(&bcf)
        .silent()
        .succeeds()
        .wrap_err("Failed to run rastair call")?;

    rastair()
        .args(["convert", "--input"])
        .arg(&bcf)
        .args(["--output-format=vcf", "-o", "-"])
        .silent()
        .succeeds()
        .wrap_err("Failed to convert to vcf")?;

    Ok(())
}

#[test]
fn convert_from_mpk() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let mpk = temp_dir.path().join("test.mpk.lz4");

    rastair()
        .args(CALL_TEST_BAM)
        .args([NO_ML, CHR19_SMALL, "-o"])
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

    rastair()
        .args(CALL_TEST_BAM)
        .args([NO_ML, CHR19_SMALL])
        .arg("-o")
        .arg(&vcf)
        .silent()
        .succeeds()
        .wrap_err("Failed to run rastair call to vcf")?;

    rastair()
        .args(CALL_TEST_BAM)
        .args([NO_ML, CHR19_SMALL])
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

    let call_args =
        CALL_TEST_BAM.iter().chain(&[CHR19_SMALL, NO_ML]).copied().collect::<Vec<&str>>().join(" ");

    let mut cmd = Command::new(insta_cmd::get_cargo_bin("/bin/bash"));
    cmd.arg("-c");
    cmd.env("NO_COLOR", "1");
    cmd.arg(format!("cargo run -q -- {call_args} --vcf | head -n1000 | cargo run -q -- convert -f bcf -F bed | head -n5"));

    assert_cmd_snapshot!(cmd);

    Ok(())
}

/// Indel lines as `(pos, ref, alt, filter)`.
fn indel_calls(vcf_text: &str) -> Vec<(String, String, String, String)> {
    vcf_content_lines(vcf_text)
        .filter_map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            let (pos, r, alt, filter) = (f.get(1)?, f.get(3)?, f.get(4)?, f.get(6)?);
            (r.len() > 1 || alt.len() > 1)
                .then(|| (pos.to_string(), r.to_string(), alt.to_string(), filter.to_string()))
        })
        .collect()
}

/// Round-tripping through `.mpk` must preserve the hard-filter verdicts.
///
/// `convert` always hands the VCF layer an `ml_threshold` (it has a default), so
/// rendering that keys off the threshold rather than off the verdict silently turns
/// every `indel_strand` / `indel_hom_ref` allele into a PASS on this path.
#[test]
fn convert_preserves_indel_hard_filter_verdicts() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let mpk = temp_dir.path().join("indels.mpk.lz4");
    let converted = temp_dir.path().join("indels.vcf");

    let mut direct = rastair()
        .args(CALL_TEST_BAM)
        .args(["--experimental-indels", "--all"])
        .output()
        .wrap_err("Failed to run rastair call")?;
    direct.succeeds()?;
    let direct_calls = indel_calls(&direct.stdout());
    assert!(!direct_calls.is_empty(), "test BAM should produce indel calls");
    assert!(
        direct_calls.iter().any(|(_, _, _, filter)| filter != "PASS"),
        "the comparison is only meaningful if some allele fails a hard filter"
    );

    rastair()
        .args(CALL_TEST_BAM)
        .args(["--experimental-indels", "-o"])
        .arg(&mpk)
        .silent()
        .succeeds()
        .wrap_err("Failed to run rastair call to mpk")?;

    rastair()
        .args(["convert", "--all", "--input"])
        .arg(&mpk)
        .arg("--output")
        .arg(&converted)
        .silent()
        .succeeds()
        .wrap_err("Failed to convert mpk to vcf")?;

    let converted_calls = indel_calls(&std::fs::read_to_string(&converted)?);
    assert_eq!(converted_calls, direct_calls);

    Ok(())
}
