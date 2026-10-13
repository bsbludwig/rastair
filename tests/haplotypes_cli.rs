mod utils;
use utils::*;

#[cfg(feature = "experimental-seqair")]
use std::collections::BTreeMap;
#[cfg(feature = "experimental-seqair")]
use utils::scenario::{Hap, Scenario};

#[cfg(feature = "experimental-seqair")]
fn diagnose_scenario(scenario: &Scenario, extra: &[&str]) -> Result<(String, String)> {
    diagnose_scenario_at(scenario, extra, None)
}

#[cfg(feature = "experimental-seqair")]
fn diagnose_scenario_at(
    scenario: &Scenario,
    extra: &[&str],
    target: Option<usize>,
) -> Result<(String, String)> {
    let dir = TempDir::new()?;
    let fasta = scenario.write_fasta(dir.path())?;
    let bam = dir.path().join("reads.bam");
    scenario.write_bam(&bam)?;
    let pairs = dir.path().join("pairs.tsv");
    let mut command = rastair();
    command
        .args(["ml", "haplotypes", "--fasta-file"])
        .arg(&fasta)
        .arg(bam)
        .arg(scenario.region())
        .arg("--pair-output")
        .arg(&pairs)
        .args(extra);
    if let Some(pos) = target {
        let bed = dir.path().join("sites.bed");
        std::fs::write(&bed, format!("scenario\t{pos}\t{}\n", pos + 1))?;
        command.arg("--sites-bed").arg(bed);
    }
    let output = command.output()?;
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
fn singleton_alternate_pairs_are_one_observation_through_bam_collection() -> Result<()> {
    use utils::scenario::Strand::{Bottom, Top};
    for (context, ref_ot, alt_ot, ref_ob, alt, non_taps) in [
        (b"ACG".as_slice(), 12, 2, 10, b'T', false),
        (b"CAG".as_slice(), 5, 0, 3, b'G', false),
        (b"CAT".as_slice(), 5, 0, 7, b'T', true),
    ] {
        let base = Scenario::new(500)
            .reference_at(119, context)
            .het(120, alt, Hap::One)
            .pairs_to(Hap::Two, Top, 100, 350, ref_ot)
            .pairs_to(Hap::One, Top, 100, 350, alt_ot)
            .pairs_to(Hap::Two, Bottom, 100, 350, ref_ob);
        for engine in ["aligned", "pair-hmm"] {
            let mut args = vec!["--max-sites=1", "--likelihood-engine", engine];
            if non_taps {
                args.push("--non-taps");
            }
            let single = base.clone().pairs_to(Hap::One, Bottom, 100, 350, 1);
            let paired = base.clone().pairs_to(Hap::One, Bottom, 100, 100, 1);
            let alt = char::from(alt).to_string();
            let (single, _) = diagnose_scenario_at(&single, &args, Some(120))?;
            let (paired, _) = diagnose_scenario_at(&paired, &args, Some(120))?;
            let single = allele_row(&single, 120, &alt)?;
            let paired = allele_row(&paired, 120, &alt)?;
            let number = |row: &BTreeMap<String, String>, field: &str| -> Result<f64> {
                Ok(row.get(field).ok_or_else(|| eyre!("missing {field}"))?.parse()?)
            };
            assert_eq!(number(&paired, "fragments")?, (ref_ot + alt_ot + ref_ob + 1) as f64);
            assert_eq!(number(&paired, "observations")?, number(&paired, "fragments")?);
            assert_eq!(number(&paired, "read_observations")?, number(&paired, "fragments")? + 1.0);
            assert_eq!(number(&paired, "overlapping_fragments")?, 1.0);
            ensure!(
                (number(&paired, "log_bayes_factor")? - number(&single, "log_bayes_factor")?).abs()
                    < 1e-8,
                "{engine} {context:?}: single={single:?}, paired={paired:?}"
            );
            let probabilities = paired
                .get("genotype_posteriors")
                .ok_or_else(|| eyre!("genotypes"))?
                .split(',')
                .map(str::parse::<f64>)
                .collect::<std::result::Result<Vec<_>, _>>()?;
            ensure!((probabilities.iter().sum::<f64>() - 1.0).abs() < 5e-12);
        }
    }
    Ok(())
}

#[test]
#[cfg(feature = "experimental-seqair")]
fn discovery_only_matches_likelihood_targets_and_refuses_incompatible_outputs() -> Result<()> {
    let scenario =
        Scenario::new(500).reference_at(119, b"AAA").het(120, b'C', Hap::One).balanced(100, 350, 4);
    let dir = TempDir::new()?;
    let fasta = scenario.write_fasta(dir.path())?;
    let bam = dir.path().join("reads.bam");
    scenario.write_bam(&bam)?;
    for engine in ["aligned", "pair-hmm"] {
        let base = || {
            let mut command = rastair();
            command
                .args(["ml", "haplotypes", "--fasta-file"])
                .arg(&fasta)
                .arg(&bam)
                .arg(scenario.region())
                .args(["--non-taps", "--likelihood-engine", engine]);
            command
        };
        let output = base().arg("--candidates-only").output()?;
        ensure!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
        let candidates = String::from_utf8(output.stdout)?;
        ensure!(candidates.lines().any(|line| line == "scenario\t121\tA"));
        ensure!(!candidates.contains("genotype_posteriors"));
        let (scored, _) =
            diagnose_scenario(&scenario, &["--non-taps", "--likelihood-engine", engine])?;
        let discovered: std::collections::BTreeSet<_> = candidates
            .lines()
            .filter(|line| line.starts_with("scenario\t"))
            .filter_map(|line| line.split('\t').nth(1))
            .collect();
        let scored: std::collections::BTreeSet<_> = scored
            .lines()
            .filter(|line| line.starts_with("scenario\t"))
            .filter_map(|line| line.split('\t').nth(1))
            .collect();
        assert_eq!(discovered, scored);
        ensure!(
            !base()
                .args(["--candidates-only", "--pair-output"])
                .arg(dir.path().join("pairs.tsv"))
                .output()?
                .status
                .success()
        );
    }
    Ok(())
}

#[test]
#[cfg(feature = "experimental-seqair")]
fn disagreeing_bam_mates_reduce_confidence_instead_of_supporting_two_alleles() -> Result<()> {
    use rust_htslib::bam::{self, Read as _};
    let scenario = Scenario::new(500).reference_at(119, b"AAA").het(120, b'G', Hap::One).pairs_to(
        Hap::One,
        utils::scenario::Strand::Bottom,
        100,
        100,
        1,
    );
    let dir = TempDir::new()?;
    let fasta = scenario.write_fasta(dir.path())?;
    let original = dir.path().join("agree.bam");
    scenario.write_bam(&original)?;
    let bed = dir.path().join("sites.bed");
    std::fs::write(&bed, "scenario\t120\t121\n")?;
    let mut bams = vec![(None, original.clone())];
    for quality in [30u8, 40] {
        let mut reader = bam::Reader::from_path(&original)?;
        let path = dir.path().join(format!("conflict_{quality}.bam"));
        let mut writer = bam::Writer::from_path(
            &path,
            &bam::Header::from_template(reader.header()),
            bam::Format::Bam,
        )?;
        for record in reader.records() {
            let mut record = record?;
            if record.is_first_in_template() {
                let mut sequence = record.seq().as_bytes();
                let mut qualities = record.qual().to_vec();
                *sequence.get_mut(20).ok_or_else(|| eyre!("focal base"))? = b'A';
                *qualities.get_mut(20).ok_or_else(|| eyre!("focal quality"))? = quality;
                let name = record.qname().to_vec();
                let cigar = record.cigar().take();
                record.set(&name, Some(&cigar), &sequence, &qualities);
            }
            writer.write(&record)?;
        }
        drop(writer);
        bam::index::build(&path, None, bam::index::Type::Bai, 1)?;
        bams.push((Some(quality), path));
    }
    for engine in ["aligned", "pair-hmm"] {
        let mut agreement = None;
        for (quality, bam) in &bams {
            let output = rastair()
                .args(["ml", "haplotypes", "--fasta-file"])
                .arg(&fasta)
                .arg(bam)
                .arg(scenario.region())
                .args(["--max-sites=1", "--non-taps", "--likelihood-engine", engine, "--sites-bed"])
                .arg(&bed)
                .output()?;
            ensure!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
            let row = allele_row(&String::from_utf8(output.stdout)?, 120, "G")?;
            let p = row.get("posterior").ok_or_else(|| eyre!("posterior"))?.parse::<f64>()?;
            assert_eq!(row.get("observations"), Some(&"1".to_owned()));
            assert_eq!(row.get("read_observations"), Some(&"2".to_owned()));
            match quality {
                None => agreement = Some(p),
                Some(q) => {
                    assert_eq!(row.get("discordant_fragments"), Some(&"1".to_owned()));
                    ensure!(p < agreement.ok_or_else(|| eyre!("agreement baseline"))?);
                    ensure!(p < if *q == 40 { 0.001 } else { 0.02 }, "{engine} Q{q}: {row:?}");
                }
            }
        }
    }
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
fn overlapping_mates_are_one_molecular_observation() -> Result<()> {
    let scenario =
        Scenario::new(400).reference_at(119, b"AAA").het(120, b'C', Hap::One).balanced(100, 110, 5);
    let (tsv, _) = diagnose_scenario(&scenario, &["--non-taps"])?;
    let row = allele_row(&tsv, 120, "C")?;
    assert_eq!(row.get("fragments"), Some(&"20".to_string()));
    assert_eq!(row.get("observations"), Some(&"20".to_string()));
    assert_eq!(row.get("read_observations"), Some(&"40".to_string()));
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

#[test]
#[cfg(feature = "experimental-seqair")]
fn pair_hmm_keeps_methylation_and_overlapping_mates_in_the_molecular_model() -> Result<()> {
    let scenario = Scenario::new(400)
        .reference_at(120, b"CG")
        .methylated(Hap::One)
        .methylated(Hap::Two)
        .balanced(100, 110, 5);
    let (tsv, _) =
        diagnose_scenario(&scenario, &["--likelihood-engine=pair-hmm", "--max-sites=1"])?;
    let row = allele_row(&tsv, 120, "T")?;
    ensure!(row.get("posterior").ok_or_else(|| eyre!("posterior"))?.parse::<f64>()? < 0.01);
    assert_eq!(row.get("realigned_fragments"), Some(&"20".to_owned()));
    assert_eq!(row.get("aligned_fallback_fragments"), Some(&"0".to_owned()));
    assert_eq!(row.get("observations"), Some(&"20".to_owned()));
    assert_eq!(row.get("read_observations"), Some(&"40".to_owned()));
    Ok(())
}

#[test]
#[cfg(feature = "experimental-seqair")]
fn pair_hmm_windows_and_indel_proposals_are_independent_of_input_chunks() -> Result<()> {
    let scenario = Scenario::new(600)
        .reference_at(119, b"AAA")
        .reference_at(169, b"AAA")
        .het(120, b'C', Hap::One)
        .het(170, b'G', Hap::Two)
        .balanced(100, 400, 4);
    let (whole, pairs_whole) = diagnose_scenario(
        &scenario,
        &["--non-taps", "--likelihood-engine=pair-hmm", "--chunk-length=600"],
    )?;
    let (chunks, pairs_chunks) = diagnose_scenario(
        &scenario,
        &["--non-taps", "--likelihood-engine=pair-hmm", "--chunk-length=150"],
    )?;
    assert_eq!(whole, chunks);
    assert_eq!(pairs_whole, pairs_chunks);
    Ok(())
}

#[test]
#[cfg(feature = "experimental-seqair")]
fn simplified_cpg_models_retain_neighbors_and_survive_chunking_through_bam() -> Result<()> {
    let scenario = Scenario::new(600)
        .reference_at(119, b"ACG")
        .reference_at(169, b"ACG")
        .het(120, b'T', Hap::One)
        .het(170, b'T', Hap::Two)
        .methylated(Hap::One)
        .methylated(Hap::Two)
        .balanced(100, 400, 6);
    for model in ["neutral", "shrinkage"] {
        for engine in ["aligned", "pair-hmm"] {
            let common = ["--methylation-model", model, "--likelihood-engine", engine];
            let mut whole_args = common.to_vec();
            whole_args.push("--chunk-length=600");
            let mut chunk_args = common.to_vec();
            chunk_args.push("--chunk-length=150");
            let (whole, pairs_whole) = diagnose_scenario_at(&scenario, &whole_args, Some(120))?;
            let (chunks, pairs_chunks) = diagnose_scenario_at(&scenario, &chunk_args, Some(120))?;
            assert_eq!(whole, chunks, "{model} {engine}");
            assert_eq!(pairs_whole, pairs_chunks, "{model} {engine}");
            let row = allele_row(&whole, 120, "T")?;
            assert_eq!(row.get("cpg_context_limited"), Some(&"false".to_owned()));
            assert_eq!(row.get("context_sites"), Some(&"3".to_owned()));
            ensure!(
                row.get("p_het").ok_or_else(|| eyre!("p_het"))?.parse::<f64>()? > 0.99,
                "{model} {engine}: {row:?}"
            );
        }
    }
    Ok(())
}

#[test]
#[cfg(feature = "experimental-seqair")]
fn a_cigar_indel_hypothesis_explains_a_false_snp_on_misaligned_reads() -> Result<()> {
    use rust_htslib::bam::{
        self, Read as _,
        record::{Cigar, CigarString},
    };
    let dir = TempDir::new()?;
    let scenario = Scenario::new(500).reference_at(120, b"ACTAGTCATG").balanced(100, 350, 5);
    let fasta = scenario.write_fasta(dir.path())?;
    let original = dir.path().join("original.bam");
    scenario.write_bam(&original)?;
    let mut reader = bam::Reader::from_path(&original)?;
    let header = bam::Header::from_template(reader.header());
    let bam = dir.path().join("insertion.bam");
    let mut writer = bam::Writer::from_path(&bam, &header, bam::Format::Bam)?;
    let mut changed = 0;
    let mut correct_names = std::collections::BTreeSet::new();
    for record in reader.records() {
        let mut record = record?;
        if record.pos() == 100 {
            let mut sequence = record.seq().as_bytes();
            sequence.insert(23, b'G');
            let cigar = if changed < 2 {
                CigarString(vec![Cigar::Match(23), Cigar::Ins(1), Cigar::Match(77)].into())
            } else {
                CigarString(vec![Cigar::Match(101)].into())
            };
            let name = record.qname().to_vec();
            if changed < 2 {
                correct_names.insert(name.clone());
            }
            record.set(&name, Some(&cigar), &sequence, &vec![40; sequence.len()]);
            changed += 1;
        }
        writer.write(&record)?;
    }
    drop(writer);
    bam::index::build(&bam, None, bam::index::Type::Bai, 1)?;
    let bed = dir.path().join("sites.bed");
    std::fs::write(&bed, "scenario\t123\t124\n")?;
    let run = |engine: &str, chunk: &str| -> Result<String> {
        let output = rastair()
            .args(["ml", "haplotypes", "--fasta-file"])
            .arg(&fasta)
            .arg(&bam)
            .args([
                "--region=scenario:1-500",
                "--non-taps",
                "--max-sites=1",
                "--min-alt-observations=2",
                "--sites-bed",
            ])
            .arg(&bed)
            .args(["--likelihood-engine", engine, "--chunk-length", chunk])
            .output()?;
        ensure!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
        Ok(String::from_utf8(output.stdout)?)
    };
    let aligned = run("aligned", "500")?;
    let hmm = run("pair-hmm", "500")?;
    assert_eq!(hmm, run("pair-hmm", "125")?);
    let aligned = allele_row(&aligned, 123, "G")?;
    let row = allele_row(&hmm, 123, "G")?;
    ensure!(aligned.get("posterior").ok_or_else(|| eyre!("posterior"))?.parse::<f64>()? > 0.99);
    ensure!(
        row.get("posterior").ok_or_else(|| eyre!("posterior"))?.parse::<f64>()? < 0.01,
        "{row:?}"
    );
    assert_eq!(row.get("haplotypes"), Some(&"8".to_owned()));
    assert_eq!(row.get("realigned_fragments"), Some(&"20".to_owned()));
    let latent_bam = dir.path().join("latent.bam");
    let mut reader = bam::Reader::from_path(&bam)?;
    let mut writer = bam::Writer::from_path(&latent_bam, &header, bam::Format::Bam)?;
    for record in reader.records() {
        let record = record?;
        if correct_names.contains(record.qname()) {
            writer.write(&record)?;
        }
    }
    drop(writer);
    bam::index::build(&latent_bam, None, bam::index::Type::Bai, 1)?;
    let output = rastair()
        .args(["ml", "haplotypes", "--fasta-file"])
        .arg(&fasta)
        .arg(&latent_bam)
        .args([
            "--region=scenario:123-123",
            "--likelihood-engine=pair-hmm",
            "--max-sites=1",
            "--min-alt-observations=2",
        ])
        .output()?;
    ensure!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    let latent = String::from_utf8(output.stdout)?;
    ensure!(
        allele_row(&latent, 122, "C")?
            .get("posterior")
            .ok_or_else(|| eyre!("posterior"))?
            .parse::<f64>()?
            .is_finite()
    );
    Ok(())
}
