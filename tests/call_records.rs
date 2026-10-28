//! What `call` writes into a record's fields, end to end, on simulated reads
//! whose truth is known.

mod utils;
use color_eyre::eyre::OptionExt as _;
use utils::scenario::{Hap, Scenario, Strand};
use utils::*;

/// The data line at `pos` (0-based).
fn record_at(vcf: &str, pos: usize) -> Option<&str> {
    let pos = (pos + 1).to_string();
    vcf.lines()
        .filter(|line| !line.starts_with('#'))
        .find(|line| line.split('\t').nth(1) == Some(&pos))
}

/// An INFO value of a VCF data line, by key.
fn info_field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    line.split('\t').nth(7)?.split(';').find_map(|entry| entry.strip_prefix(key)?.strip_prefix('='))
}

fn integers(value: &str) -> Result<Vec<u32>> {
    value
        .split(',')
        .map(|v| v.parse().wrap_err_with(|| format!("not an integer: {value}")))
        .collect()
}

/// Two homozygous deletions 18 bp apart are one pair-HMM cluster. Each has
/// reads of its own, which end before the other or start after it, and a few
/// reads span both. A read over one deletion only says as much about it as a
/// read over both, so each record's `AD` and `DP` count every read over its
/// deletion.
#[test]
fn indel_depth_counts_reads_over_one_event_of_a_cluster() -> Result<()> {
    const A: usize = 300;
    const B: usize = 318;
    const ONLY_A: usize = 8;
    const ONLY_B: usize = 8;
    const BOTH: usize = 2;
    let pairs = |scenario: Scenario, left: usize, copies: usize| {
        scenario.pairs(Hap::One, Strand::Top, left, copies).pairs(
            Hap::Two,
            Strand::Bottom,
            left,
            copies,
        )
    };
    // Each read is 100 bases: from 215 it ends at 316, past A's deletion and
    // before B's anchor; from 305 it starts past A's deletion.
    let mut scenario = Scenario::new(800);
    for hap in [Hap::One, Hap::Two] {
        scenario = scenario.deletion(A, 2, hap).deletion(B, 2, hap);
    }
    let scenario = pairs(pairs(pairs(scenario, 215, ONLY_A / 2), 305, ONLY_B / 2), 260, BOTH / 2);
    let vcf = call_scenario(&scenario, &[])?;

    for (pos, over) in [(A, ONLY_A + BOTH), (B, ONLY_B + BOTH)] {
        let line = record_at(&vcf, pos).ok_or_else(|| eyre!("no record at {pos}:\n{vcf}"))?;
        let info = |key| info_field(line, key).ok_or_else(|| eyre!("no {key} in {line}"));
        assert_eq!(format_field(line, "GT"), Some("1/1"), "{line}");
        let depth: u32 = info("DP")?.parse()?;
        let ad = integers(info("AD")?)?;
        assert_eq!(
            (depth, ad.as_slice()),
            (u32::try_from(over)?, [0, u32::try_from(over)?].as_slice()),
            "DP and AD at {pos} must count the {over} reads over its deletion: {line}"
        );
    }
    Ok(())
}

/// Where `CCG` becomes `CGG` on the only haplotype the reads come from, the
/// middle base is the `C` of the reference CpG and the `G` of a de-novo one.
/// No read shows the reference CpG, and the bottom-strand reads show the
/// de-novo one fully methylated.
fn reference_and_denovo_cpg_at(site: usize) -> Scenario {
    Scenario::new(1000)
        .reference_at(site - 1, b"CCG")
        .het(site, b'G', Hap::One)
        .methylated(Hap::One)
        .pairs(Hap::One, Strand::Top, site - 50, 8)
        .pairs(Hap::One, Strand::Bottom, site - 50, 8)
}

/// `M5mC`, `DPM5mC` and `ADM5mC` hold one value per CpG a site takes part
/// in, the reference CpG first and the de-novo one second
/// (`docs/src/calling/methylation.md`), whether a read shows it or not.
#[test]
fn methylation_fields_keep_a_slot_per_cpg() -> Result<()> {
    const SITE: usize = 300;
    let vcf = call_scenario(&reference_and_denovo_cpg_at(SITE), &[])?;
    let line = record_at(&vcf, SITE).ok_or_else(|| eyre!("no record at {SITE}:\n{vcf}"))?;
    let field = |key| format_field(line, key).ok_or_else(|| eyre!("no {key} in {line}"));
    assert_eq!(line.split('\t').nth(3), Some("C"), "{line}");
    assert_eq!(field("GT")?, "1/1", "{line}");

    let slots = |key| -> Result<Vec<String>> {
        Ok(field(key)?.split(',').map(str::to_owned).collect::<Vec<_>>())
    };
    let [m5mc, depth, modified] = ["M5mC", "DPM5mC", "ADM5mC"].map(slots);
    let (m5mc, depth, modified) = (m5mc?, depth?, modified?);
    assert_eq!(m5mc.len(), 2, "M5mC: one slot per CpG: {line}");
    assert_eq!(m5mc.get(1).map(String::as_str), Some("1"), "de-novo β second: {line}");
    assert_eq!(depth.first().map(String::as_str), Some("0"), "{line}");
    assert_eq!(depth.len(), 2, "{line}");
    assert_eq!(modified.len(), 2, "{line}");
    Ok(())
}

/// A CpG no read shows has no β: `M5mC` is missing (`.`), not 0, and the
/// BED's `beta_est` is empty (`docs/src/formats/bed.md`), so that averaging
/// does not count it unmethylated. `--all` writes uncovered reference CpGs.
#[test]
fn a_cpg_without_evidence_has_no_beta() -> Result<()> {
    const SITE: usize = 300;
    const UNCOVERED: usize = 700;
    let scenario = reference_and_denovo_cpg_at(SITE).reference_at(UNCOVERED, b"CG");
    let dir = TempDir::new()?;
    let bed = dir.path().join("out.bed");
    let bed_arg = bed.to_str().ok_or_eyre("a UTF-8 path")?;
    let vcf = call_scenario_in(dir.path(), &scenario, &["--all", "--bed", bed_arg])?;
    let bed = std::fs::read_to_string(&bed)?;

    for pos in [UNCOVERED, UNCOVERED + 1] {
        let line = record_at(&vcf, pos).ok_or_else(|| eyre!("no record at {pos}:\n{vcf}"))?;
        assert_eq!(format_field(line, "DPM5mC"), Some("0"), "{line}");
        assert_eq!(format_field(line, "M5mC"), Some("."), "{line}");
        let start = pos.to_string();
        let row = bed
            .lines()
            .find(|row| row.split('\t').nth(1) == Some(&start))
            .ok_or_else(|| eyre!("no BED row at {pos}:\n{bed}"))?;
        assert_eq!(row.split('\t').nth(4), Some(""), "beta_est: {row}");
    }
    // The reference CpG's slot at the de-novo site, which no read shows.
    let line = record_at(&vcf, SITE).ok_or_else(|| eyre!("no record at {SITE}:\n{vcf}"))?;
    let m5mc = format_field(line, "M5mC").ok_or_else(|| eyre!("no M5mC in {line}"))?;
    assert_eq!(m5mc.split(',').next(), Some("."), "{line}");
    Ok(())
}

/// A methylation candidate whose conversion more reads show than
/// `--m-max-coverage` is `m_highDp`, not called: a pileup that deep is a
/// repeat or a collapsed duplication. At the cap it is judged as usual.
#[test]
fn a_methylation_candidate_over_the_cap_is_m_high_dp() -> Result<()> {
    const SITE: usize = 300;
    // Eight top-strand reads show the methylated C as T.
    let scenario = Scenario::new(800)
        .reference_at(SITE, b"CG")
        .methylated(Hap::One)
        .pairs(Hap::One, Strand::Top, SITE - 50, 8)
        .pairs(Hap::One, Strand::Bottom, SITE - 50, 8);
    let filter_of_t = |cap: &str| -> Result<String> {
        let vcf = call_scenario(&scenario, &["--all", "--m-max-coverage", cap])?;
        let pos = (SITE + 1).to_string();
        let line = vcf
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .find(|cols| cols.get(1) == Some(&pos.as_str()) && cols.get(4) == Some(&"T"))
            .ok_or_else(|| eyre!("no C>T record at {SITE}:\n{vcf}"))?;
        Ok(line.get(6).ok_or_eyre("a FILTER column")?.to_string())
    };
    assert_eq!(filter_of_t("7")?, "m_highDp");
    assert_ne!(filter_of_t("8")?, "m_highDp");
    Ok(())
}
