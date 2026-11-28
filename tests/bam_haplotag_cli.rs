//! Haplotype tags from `bam`, end to end, on scenarios whose haplotypes are
//! known: `call` phases the scenario, `bam` tags its reads from that VCF, and
//! every tagged read is checked against the haplotype it was written from.

mod utils;
use rust_htslib::bam::{self, Read as _, record::Aux};
use std::collections::BTreeMap;
use utils::{
    scenario::{Hap, Scenario, Strand},
    *,
};

/// The calls of a scenario, in the formats `bam` reads.
struct Called {
    dir: TempDir,
    scenario: Scenario,
}

impl Called {
    fn new(scenario: Scenario) -> Result<Self> {
        let dir = TempDir::new()?;
        let fasta = scenario.write_fasta(dir.path())?;
        scenario.write_bam(&dir.path().join("reads.bam"))?;
        let called = Self { dir, scenario };
        rastair()
            .arg("call")
            .arg("--fasta-file")
            .arg(&fasta)
            .arg(called.path("reads.bam"))
            .arg(called.scenario.region())
            .arg("--vcf")
            .arg(called.path("calls.vcf.gz"))
            .arg("--bed")
            .arg(called.path("calls.bed.gz"))
            .silent()
            .succeeds()
            .wrap_err("call")?;
        Ok(called)
    }

    fn path(&self, name: &str) -> std::path::PathBuf {
        self.dir.path().join(name)
    }

    /// Run `bam <mode>` with `calls` and return its records. With `batch`,
    /// the input is cut into batches of that many records (debug builds).
    fn rewrite(&self, mode: &str, calls: &str, batch: Option<usize>) -> Result<Vec<bam::Record>> {
        let out = self.path("out.bam");
        let mut command = rastair();
        if let Some(batch) = batch {
            command.env("RASTAIR_BAM_BATCH_RECORDS", batch.to_string());
        }
        command
            .args(["bam", mode, "--fasta-file"])
            .arg(self.path("ref.fa"))
            .arg(self.path("reads.bam"))
            .arg(self.path(calls))
            .arg("-o")
            .arg(&out)
            .silent()
            .succeeds()
            .wrap_err("bam")?;
        let mut reader = bam::Reader::from_path(&out)?;
        reader.records().map(|record| record.wrap_err("read record")).collect()
    }

    /// The tags `bam` gives every read from the VCF, checked against the
    /// scenario: within a phase set, `HP` 1 holds reads of one haplotype only
    /// and `HP` 2 reads of the other only, and mates agree.
    fn tags(&self, mode: &str, batch: Option<usize>) -> Result<Tags> {
        Tags::of(&self.scenario, &self.rewrite(mode, "calls.vcf.gz", batch)?)
    }
}

fn int_tag(record: &bam::Record, tag: &[u8]) -> Result<Option<i64>> {
    Ok(Some(match record.aux(tag) {
        Ok(Aux::U8(value)) => i64::from(value),
        Ok(Aux::I8(value)) => i64::from(value),
        Ok(Aux::U16(value)) => i64::from(value),
        Ok(Aux::I16(value)) => i64::from(value),
        Ok(Aux::U32(value)) => i64::from(value),
        Ok(Aux::I32(value)) => i64::from(value),
        Err(bam::BamError::AuxTagNotFound) => return Ok(None),
        other => bail!("{} is not an integer tag: {other:?}", tag.escape_ascii()),
    }))
}

/// The `(PS, HP)` of each tagged read pair, and the names of the untagged.
#[derive(Debug, PartialEq, Eq)]
struct Tags {
    /// Per phase set (1-based, as `PS` is written), the scenario haplotype
    /// `HP` 1 and `HP` 2 turned out to be.
    sets: BTreeMap<i64, BTreeMap<i64, Hap>>,
    tagged: BTreeMap<String, (i64, i64)>,
    untagged: BTreeSet<String>,
}

impl Tags {
    fn of(scenario: &Scenario, records: &[bam::Record]) -> Result<Self> {
        let truth = scenario.read_haplotypes();
        ensure!(records.len() == 2 * truth.len(), "{} records", records.len());

        let mut tags =
            Self { sets: BTreeMap::new(), tagged: BTreeMap::new(), untagged: BTreeSet::new() };
        for record in records {
            let name = String::from_utf8(record.qname().to_vec())?;
            let hap = *truth.get(&name).ok_or_else(|| eyre!("unknown read {name}"))?;
            let tag = match (int_tag(record, b"PS")?, int_tag(record, b"HP")?) {
                (Some(set), Some(hp)) => Some((set, hp)),
                (None, None) => None,
                other => bail!("{name} has one of HP and PS only: {other:?}"),
            };
            let Some((set, hp)) = tag else {
                ensure!(!tags.tagged.contains_key(&name), "only one mate of {name} is tagged");
                tags.untagged.insert(name);
                continue;
            };
            ensure!(hp == 1 || hp == 2, "{name} has HP {hp}");
            ensure!(!tags.untagged.contains(&name), "only one mate of {name} is tagged");
            if let Some(mate) = tags.tagged.insert(name.clone(), (set, hp)) {
                ensure!(mate == (set, hp), "the mates of {name} disagree: {mate:?}, ({set}, {hp})");
            }
            let seen = *tags.sets.entry(set).or_default().entry(hp).or_insert(hap);
            ensure!(seen == hap, "PS {set} HP {hp} holds reads of both haplotypes ({name})");
        }
        for (set, haps) in &tags.sets {
            ensure!(
                haps.get(&1) != haps.get(&2),
                "PS {set} has one haplotype's reads under both HP"
            );
        }
        Ok(tags)
    }

    /// The scenario haplotype tagged `HP` 1 and `HP` 2 in `set`.
    fn haplotypes(&self, set: usize) -> Option<[Hap; 2]> {
        let haps = self.sets.get(&i64::try_from(set).ok()?)?;
        Some([*haps.get(&1)?, *haps.get(&2)?])
    }
}

/// Two hets under the left read, none under its mate.
const LEFT: usize = 100;
const RIGHT: usize = 350;
const SITE_A: usize = 120;
const SITE_B: usize = 180;

fn two_sites_under_the_left_read() -> Scenario {
    Scenario::new(1000).snv(SITE_A, Hap::One).snv(SITE_B, Hap::Two).balanced(LEFT, RIGHT, 5)
}

/// Every read pair is tagged with the block's `PS`, the mate over no site
/// included, and `HP` 1 is the haplotype left of the `|`: the block opens
/// `0|1` at `SITE_A`, whose alt is on haplotype one, so `HP` 1 is haplotype
/// two.
#[test]
fn read_pairs_are_tagged_with_their_haplotype_and_phase_set() -> Result<()> {
    let called = Called::new(two_sites_under_the_left_read())?;
    for mode in ["standard", "legacy"] {
        let tags = called.tags(mode, None)?;
        assert_eq!(tags.untagged, BTreeSet::new(), "{mode}");
        assert_eq!(
            tags.sets.keys().copied().collect::<Vec<_>>(),
            [i64::try_from(SITE_A + 1)?],
            "{mode}"
        );
        assert_eq!(tags.haplotypes(SITE_A + 1), Some([Hap::Two, Hap::One]), "{mode}");
    }
    Ok(())
}

/// A BED file has no phase: reads get no haplotype tags from it, and the
/// methylation tags are the ones the VCF gives.
#[test]
fn bed_calls_give_the_same_reads_without_haplotype_tags() -> Result<()> {
    let scenario = two_sites_under_the_left_read()
        .reference_at(140, b"CG")
        .reference_at(400, b"CG")
        .methylated(Hap::One);
    let called = Called::new(scenario)?;

    for mode in ["standard", "legacy"] {
        let from_bed = called.rewrite(mode, "calls.bed.gz", None)?;
        let mut from_vcf = called.rewrite(mode, "calls.vcf.gz", None)?;
        ensure!(from_vcf.iter().any(|record| record.aux(b"HP").is_ok()), "no read was tagged");
        let methylation_tag: &[u8] = if mode == "standard" { b"MM" } else { b"XM" };
        ensure!(
            from_bed.iter().any(|record| record.aux(methylation_tag).is_ok()),
            "no methylation tags"
        );

        for record in &from_bed {
            ensure!(record.aux(b"HP").is_err() && record.aux(b"PS").is_err(), "tagged from a BED");
        }
        for record in &mut from_vcf {
            record.remove_aux(b"HP")?;
            record.remove_aux(b"PS")?;
        }
        ensure!(from_bed == from_vcf, "{mode}: records differ beyond HP and PS");
    }
    Ok(())
}

/// A `C>T` het at a methylated CpG: on the top strand the reference
/// haplotype's methylated `C` reads `T`, like the alt. Those reads may not
/// vote there, so no read lands on the wrong haplotype (which [`Tags::of`]
/// checks), and the pairs are still tagged through the other site.
#[test]
fn a_converted_read_is_not_tagged_as_the_alt() -> Result<()> {
    const CPG: usize = 150;
    let scenario = Scenario::new(1000)
        .reference_at(CPG, b"CG")
        .het(CPG, b'T', Hap::One)
        .snv(SITE_B, Hap::Two)
        .methylated(Hap::Two)
        .balanced(LEFT, RIGHT, 8);
    let called = Called::new(scenario)?;

    let tags = called.tags("standard", None)?;
    assert_eq!(tags.untagged, BTreeSet::new());
    assert_eq!(tags.sets.len(), 1, "{:?}", tags.sets);
    Ok(())
}

/// Read pairs whose only het is an indel are tagged from it, whichever
/// haplotype carries it.
#[test]
fn reads_are_tagged_from_an_indel() -> Result<()> {
    const INSERTION: usize = 130;
    const DELETION: usize = 400;
    let scenario = Scenario::new(1000)
        .insertion(INSERTION, b"TTA", Hap::One)
        .deletion(DELETION, 2, Hap::Two)
        .balanced(LEFT, RIGHT, 8);
    let called = Called::new(scenario)?;

    let tags = called.tags("legacy", None)?;
    assert_eq!(tags.untagged, BTreeSet::new());
    assert_eq!(tags.haplotypes(INSERTION + 1), Some([Hap::Two, Hap::One]), "{:?}", tags.sets);
    Ok(())
}

/// A read pair is tagged the same however the input is cut into batches. In
/// batches of seven records the left reads, which cover no het, come batches
/// before their mates, and the reads further on let the first batches be
/// written while later ones still wait.
#[test]
fn tags_do_not_depend_on_the_batches() -> Result<()> {
    const START: usize = 900;
    let mate = START + 250;
    let scenario = Scenario::new(3000)
        .snv(mate + 20, Hap::One)
        .snv(mate + 70, Hap::Two)
        .balanced(START, mate, 6)
        // Elsewhere, so each segment of `call` has a depth like the run's.
        .balanced(2000, 2300, 6)
        .snv(2020, Hap::One)
        .snv(2070, Hap::Two);
    let called = Called::new(scenario)?;

    for mode in ["standard", "legacy"] {
        let whole = called.tags(mode, None)?;
        assert_eq!(whole.untagged, BTreeSet::new(), "{mode}");
        for batch in [1, 7, 24] {
            assert_eq!(whole, called.tags(mode, Some(batch))?, "{mode}, batches of {batch}");
        }
    }
    Ok(())
}

/// Read pairs over no phased site are left untagged.
#[test]
fn a_read_pair_over_no_phased_site_has_no_tag() -> Result<()> {
    let scenario = two_sites_under_the_left_read().pairs(Hap::One, Strand::Top, 600, 3);
    let called = Called::new(scenario)?;

    let tags = called.tags("standard", None)?;
    let truth = called.scenario.read_haplotypes();
    assert_eq!(tags.untagged.len(), 3, "{:?}", tags.untagged);
    assert_eq!(tags.tagged.len(), truth.len() - 3);
    Ok(())
}

/// An uncompressed VCF has no index to read a segment's calls with.
#[test]
fn an_uncompressed_vcf_is_refused_with_a_suggestion() -> Result<()> {
    let called = Called::new(two_sites_under_the_left_read())?;
    std::fs::write(called.path("calls.vcf"), "##fileformat=VCFv4.5\n")?;

    let output = rastair()
        .args(["bam", "standard", "--fasta-file"])
        .arg(called.path("ref.fa"))
        .arg(called.path("reads.bam"))
        .arg(called.path("calls.vcf"))
        .arg("-o")
        .arg(called.path("out.bam"))
        .output()?;
    ensure!(!output.status.success(), "an uncompressed VCF was accepted");
    let stderr = output.stderr();
    ensure!(stderr.contains("uncompressed VCF"), "{stderr}");
    ensure!(stderr.contains("bcftools index"), "{stderr}");
    Ok(())
}
