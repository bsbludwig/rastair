mod utils;
use utils::*;

#[cfg(feature = "experimental-seqair")]
use std::collections::BTreeMap;
#[cfg(feature = "experimental-seqair")]
use utils::scenario::{Hap, Scenario};

#[cfg(feature = "experimental-seqair")]
fn diagnose_scenario(scenario: &Scenario, extra: &[&str]) -> Result<(String, String)> {
    let dir = TempDir::new()?;
    let fasta = scenario.write_fasta(dir.path())?;
    let bam = dir.path().join("reads.bam");
    scenario.write_bam(&bam)?;
    let pairs = dir.path().join("pairs.tsv");
    let output = rastair()
        .args(["ml", "haplotypes", "--fasta-file"])
        .arg(fasta)
        .arg(bam)
        .arg(scenario.region())
        .arg("--pair-output")
        .arg(&pairs)
        .args(extra)
        .output()?;
    ensure!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    Ok((String::from_utf8(output.stdout)?, std::fs::read_to_string(pairs)?))
}

#[cfg(feature = "experimental-seqair")]
fn allele_row(tsv: &str, pos: usize, alt: &str) -> Result<BTreeMap<String, String>> {
    let mut lines = tsv.lines().filter(|line| !line.starts_with('#'));
    let header: Vec<_> = lines.next().ok_or_else(|| eyre!("missing header"))?.split('\t').collect();
    for line in lines {
        let row: BTreeMap<_, _> = header
            .iter()
            .copied()
            .zip(line.split('\t'))
            .map(|(key, value)| (key.to_string(), value.to_string()))
            .collect();
        if row.get("pos").is_some_and(|value| value == &(pos + 1).to_string())
            && row.get("alt").is_some_and(|value| value == alt)
        {
            return Ok(row);
        }
    }
    bail!("missing {pos}/{alt} row")
}

#[test]
#[cfg(feature = "experimental-seqair")]
fn full_read_pairs_support_competing_haplotype_likelihoods() -> Result<()> {
    let scenario = Scenario::new(600)
        .reference_at(119, b"AAA")
        .reference_at(419, b"AAA")
        .het(120, b'C', Hap::One)
        .het(420, b'G', Hap::Two)
        .balanced(100, 400, 6);
    let (tsv, pairs) = diagnose_scenario(&scenario, &["--non-taps", "--window-radius=400"])?;
    let row = allele_row(&tsv, 120, "C")?;
    ensure!(row.get("p_het").ok_or_else(|| eyre!("p_het"))?.parse::<f64>()? > 0.999);
    ensure!(row.get("context_sites").is_some_and(|value| value == "2"));
    let best = pairs
        .lines()
        .skip(1)
        .filter_map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            if fields.get(1)? != &"121" {
                return None;
            }
            Some((fields.get(7)?.parse::<f64>().ok()?, fields))
        })
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .ok_or_else(|| eyre!("pair"))?;
    assert!(best.0 > 0.999);
    assert_eq!(best.1.get(2), Some(&"121,421"));
    let a = best.1.get(3).ok_or_else(|| eyre!("haplotype a"))?;
    let b = best.1.get(4).ok_or_else(|| eyre!("haplotype b"))?;
    assert!((a == &"CA" && b == &"AG") || (a == &"AG" && b == &"CA"));
    Ok(())
}

#[test]
#[cfg(feature = "experimental-seqair")]
fn taps_methylation_is_distinct_from_an_ordinary_wgs_snp() -> Result<()> {
    let scenario = Scenario::new(500)
        .reference_at(120, b"CG")
        .methylated(Hap::One)
        .methylated(Hap::Two)
        .balanced(100, 350, 6);
    let (taps, _) = diagnose_scenario(&scenario, &["--max-sites=1"])?;
    let (wgs, _) = diagnose_scenario(&scenario, &["--max-sites=1", "--non-taps"])?;
    let taps = allele_row(&taps, 120, "T")?;
    let wgs = allele_row(&wgs, 120, "T")?;
    ensure!(taps.get("posterior").ok_or_else(|| eyre!("posterior"))?.parse::<f64>()? < 0.01);
    ensure!(wgs.get("posterior").ok_or_else(|| eyre!("posterior"))?.parse::<f64>()? > 0.99);
    Ok(())
}

#[test]
#[cfg(feature = "experimental-seqair")]
fn diagnostic_windows_are_independent_of_input_chunks() -> Result<()> {
    let scenario = Scenario::new(600)
        .reference_at(119, b"AAA")
        .reference_at(169, b"AAA")
        .het(120, b'C', Hap::One)
        .het(170, b'G', Hap::Two)
        .balanced(100, 400, 5);
    let (whole, _) = diagnose_scenario(&scenario, &["--non-taps", "--chunk-length=600"])?;
    let (chunks, _) = diagnose_scenario(&scenario, &["--non-taps", "--chunk-length=150"])?;
    assert_eq!(whole, chunks);
    Ok(())
}

#[test]
#[cfg(feature = "experimental-seqair")]
fn overlapping_mates_are_two_observations_of_one_fragment() -> Result<()> {
    let scenario =
        Scenario::new(400).reference_at(119, b"AAA").het(120, b'C', Hap::One).balanced(100, 110, 5);
    let (tsv, _) = diagnose_scenario(&scenario, &["--non-taps"])?;
    let row = allele_row(&tsv, 120, "C")?;
    assert_eq!(row.get("fragments"), Some(&"20".to_string()));
    assert_eq!(row.get("observations"), Some(&"40".to_string()));
    Ok(())
}

#[test]
#[cfg(feature = "experimental-seqair")]
fn a_denovo_hypothesis_is_scored_without_an_observed_c() -> Result<()> {
    let scenario = Scenario::new(400).reference_at(120, b"TG").balanced(100, 250, 3);
    let (tsv, _) = diagnose_scenario(&scenario, &["--max-sites=1"])?;
    let row = allele_row(&tsv, 120, "C")?;
    assert_eq!(row.get("ref"), Some(&"T".to_string()));
    ensure!(row.get("posterior").ok_or_else(|| eyre!("posterior"))?.parse::<f64>()?.is_finite());
    Ok(())
}

#[test]
#[cfg(feature = "experimental-seqair")]
fn first_and_last_reference_positions_have_the_correct_coordinates() -> Result<()> {
    for (position, left, right) in [(0usize, 0, 250), (399, 150, 300)] {
        let context = position.saturating_sub(1);
        let scenario = Scenario::new(400)
            .reference_at(context, b"AAA")
            .het(position, b'C', Hap::One)
            .balanced(left, right, 4);
        let (tsv, _) = diagnose_scenario(&scenario, &["--non-taps"])?;
        let row = allele_row(&tsv, position, "C")?;
        ensure!(row.get("p_het").ok_or_else(|| eyre!("p_het"))?.parse::<f64>()? > 0.99);
    }
    Ok(())
}

#[test]
#[cfg(feature = "experimental-seqair")]
fn explicit_targets_include_zero_coverage_without_inventing_evidence() -> Result<()> {
    let dir = TempDir::new()?;
    let bed = dir.path().join("sites.bed");
    std::fs::write(&bed, "scenario\t220\t221\n")?;
    let scenario = Scenario::new(400).balanced(100, 250, 3);
    let (tsv, _) = diagnose_scenario(
        &scenario,
        &["--non-taps", "--sites-bed", bed.to_str().ok_or_else(|| eyre!("BED path"))?],
    )?;
    let mut lines = tsv.lines().filter(|line| !line.starts_with('#'));
    let header: Vec<_> = lines.next().ok_or_else(|| eyre!("header"))?.split('\t').collect();
    let posterior =
        header.iter().position(|key| key == &"posterior").ok_or_else(|| eyre!("posterior"))?;
    let fragments =
        header.iter().position(|key| key == &"fragments").ok_or_else(|| eyre!("fragments"))?;
    let data: Vec<_> = lines.collect();
    assert_eq!(data.len(), 3);
    for line in data {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.get(fragments), Some(&"0"));
        let p = fields.get(posterior).ok_or_else(|| eyre!("posterior value"))?.parse::<f64>()?;
        assert!((p - 0.001 * 4.0 / 9.0).abs() < 1e-12);
    }
    Ok(())
}

#[test]
#[cfg(feature = "experimental-seqair")]
fn diagnostics_reject_an_output_that_aliases_an_input() -> Result<()> {
    let dir = TempDir::new()?;
    let scenario = Scenario::new(400).balanced(100, 250, 3);
    let fasta = scenario.write_fasta(dir.path())?;
    let bam = dir.path().join("reads.bam");
    scenario.write_bam(&bam)?;
    let before = std::fs::read(&bam)?;
    let output = rastair()
        .args(["ml", "haplotypes", "--fasta-file"])
        .arg(fasta)
        .arg(&bam)
        .arg(scenario.region())
        .arg("--output")
        .arg(&bam)
        .output()?;
    ensure!(!output.status.success());
    ensure!(String::from_utf8_lossy(&output.stderr).contains("overwrite an input file"));
    assert_eq!(before, std::fs::read(bam)?);
    Ok(())
}

#[test]
#[cfg(not(feature = "experimental-seqair"))]
fn diagnostics_refuse_a_backend_without_full_read_collection() -> Result<()> {
    let output = rastair()
        .args([
            "ml",
            "haplotypes",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            CHR19_SMALL,
        ])
        .output()?;
    ensure!(!output.status.success());
    ensure!(
        String::from_utf8_lossy(&output.stderr)
            .contains("ml haplotypes needs the `experimental-seqair` backend")
    );
    Ok(())
}
