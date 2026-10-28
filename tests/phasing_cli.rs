//! Phasing in `call`, end to end, on scenarios whose haplotypes are known: every
//! run is checked against the truth by [`Scenario::check_phasing`], and each
//! test then states which sites had to end up phased together.

mod utils;
use std::collections::BTreeMap;
use utils::scenario::Strand;
use utils::{
    scenario::{Hap, Phasing, Scenario},
    *,
};

/// Two sites read pairs link, a left read over one and its mate over the
/// other.
const LEFT: usize = 100;
const RIGHT: usize = 350;
const SITE_A: usize = 120;
const SITE_B: usize = 420;

fn two_sites_in_trans() -> Scenario {
    Scenario::new(1000).snv(SITE_A, Hap::One).snv(SITE_B, Hap::Two).balanced(LEFT, RIGHT, 5)
}

fn blocks<const N: usize>(blocks: [&[usize]; N]) -> BTreeMap<usize, Vec<usize>> {
    blocks.into_iter().filter_map(|sites| Some((*sites.first()?, sites.to_vec()))).collect()
}

/// Two heterozygous SNVs that only read pairs connect come out as one block,
/// filed under the first site's position.
#[test]
fn read_pairs_link_two_hets_into_one_block() -> Result<()> {
    let scenario = two_sites_in_trans();
    let vcf = call_scenario(&scenario, &[])?;

    let phasing = scenario.check_phasing(&vcf)?;
    assert_eq!(phasing, Phasing { blocks: blocks([&[SITE_A, SITE_B]]), unphased: vec![] });
    assert!(vcf.contains("##FORMAT=<ID=PS,"), "PS must be declared");
    Ok(())
}

/// A block is solved on the whole segment, overlap included, and then cut to
/// the core. With 1000-base segments `HUB` is in the first segment's core and
/// the second's overlap, and `SPOKE_1` and `SPOKE_2`, in the second's core, are
/// linked only through it. So the second segment phases the spokes as one
/// block, and the hub, alone in its core, stays unphased. Phasing after the
/// trim would leave all three unphased. The overlap starts before the hub's
/// reads do: the SNV engine genotypes a site only on reads inside the segment.
#[test]
fn a_block_is_solved_across_the_overlap_and_cut_to_the_core() -> Result<()> {
    const SHIFT: usize = 700;
    const HUB: usize = SHIFT + 120;
    const SPOKE_1: usize = SHIFT + 370;
    const SPOKE_2: usize = SHIFT + 480;
    let scenario = Scenario::new(2700)
        .snv(HUB, Hap::One)
        .snv(SPOKE_1, Hap::Two)
        .snv(SPOKE_2, Hap::One)
        .balanced(SHIFT + LEFT, SHIFT + 350, 5)
        .balanced(SHIFT + LEFT, SHIFT + 460, 5)
        // Reads at the far end, outside the second segment, enough to bring
        // the run's reads per base near the second segment's: with fewer, the
        // depth guard keeps the hub, at twice the spokes' depth, from bridging
        // them.
        .balanced(2560, 2600, 25);

    // The overlap reaches back past the hub's reads at `SHIFT + LEFT`: the
    // pair-HMM skips a read that leaves the segment, and a het none of whose
    // reads it scores casts no votes.
    let vcf = call_scenario(&scenario, &["--segment-max-length=1000", "--segment-overlap=250"])?;

    let phasing = scenario.check_phasing(&vcf)?;
    assert_eq!(phasing, Phasing { blocks: blocks([&[SPOKE_1, SPOKE_2]]), unphased: vec![HUB] });
    Ok(())
}

/// Two well-linked pairs of hets, `[A, B]` and `[C, D]`, and `bridges` read
/// pairs from haplotype one reaching from `B` to `C`.
fn two_blocks_bridged(bridges: usize) -> (Scenario, [usize; 4]) {
    const SITES: [usize; 4] = [120, 420, 720, 970];
    let [a, b, c, d] = SITES;
    let scenario = Scenario::new(1300)
        .snv(a, Hap::One)
        .snv(b, Hap::Two)
        .snv(c, Hap::One)
        .snv(d, Hap::One)
        .balanced(100, 350, 5)
        .balanced(700, 950, 5)
        .pairs_to(Hap::One, Strand::Top, 400, 700, bridges);
    (scenario, SITES)
}

/// One read pair is a link wrong 0.26 % of the time, too weak to join two
/// blocks on; two agreeing pairs are enough.
#[test]
fn a_block_is_cut_where_one_read_pair_is_its_only_link() -> Result<()> {
    let (scenario, [a, b, c, d]) = two_blocks_bridged(1);
    let vcf = call_scenario(&scenario, &[])?;
    let phasing = scenario.check_phasing(&vcf)?;
    assert_eq!(phasing, Phasing { blocks: blocks([&[a, b], &[c, d]]), unphased: vec![] });
    assert!(vcf.contains("##FORMAT=<ID=PQ,"), "PQ must be declared");

    let (scenario, sites) = two_blocks_bridged(2);
    let vcf = call_scenario(&scenario, &[])?;
    let phasing = scenario.check_phasing(&vcf)?;
    assert_eq!(phasing, Phasing { blocks: blocks([&sites]), unphased: vec![] });
    Ok(())
}

/// At Q20 a single read pair is enough to keep two blocks together.
#[test]
fn a_lower_switch_quality_keeps_a_one_read_pair_link() -> Result<()> {
    let (scenario, sites) = two_blocks_bridged(1);
    let vcf = call_scenario(&scenario, &["--phase-min-switch-quality", "20"])?;
    let phasing = scenario.check_phasing(&vcf)?;
    assert_eq!(phasing, Phasing { blocks: blocks([&sites]), unphased: vec![] });
    Ok(())
}

/// Six hets within 50 bases, between two well-linked blocks and tied to both
/// by many read pairs: too dense to be trusted with a bridge, so the blocks
/// stay apart and each dense site joins one of them at most.
#[test]
fn a_dense_run_of_hets_does_not_bridge_two_blocks() -> Result<()> {
    const LEFT_BLOCK: [usize; 2] = [120, 300];
    const DENSE: [usize; 6] = [500, 510, 520, 530, 540, 550];
    const RIGHT_BLOCK: [usize; 2] = [760, 900];
    let haps = [Hap::One, Hap::Two, Hap::Two, Hap::One, Hap::Two, Hap::One, Hap::One, Hap::Two];
    let scenario = LEFT_BLOCK
        .into_iter()
        .chain(DENSE)
        .chain(RIGHT_BLOCK)
        .zip(haps.into_iter().cycle())
        .fold(Scenario::new(1100), |scenario, (pos, hap)| scenario.snv(pos, hap))
        .balanced(100, 280, 5)
        .balanced(280, 480, 5)
        .balanced(480, 740, 5)
        .balanced(740, 880, 5);

    let vcf = call_scenario(&scenario, &[])?;
    let phasing = scenario.check_phasing(&vcf)?;

    let block_of = |site: usize| phasing.blocks.values().position(|block| block.contains(&site));
    assert_eq!(phasing.blocks.len(), 2, "{phasing:?}");
    assert!(
        block_of(LEFT_BLOCK[0]).is_some() && block_of(LEFT_BLOCK[0]) == block_of(LEFT_BLOCK[1])
    );
    assert!(
        block_of(RIGHT_BLOCK[0]).is_some() && block_of(RIGHT_BLOCK[0]) == block_of(RIGHT_BLOCK[1])
    );
    assert_ne!(block_of(LEFT_BLOCK[0]), block_of(RIGHT_BLOCK[0]), "{phasing:?}");
    Ok(())
}

/// The SNV the methylation-confounded site is linked to.
const SNV: usize = 130;
/// A C/T heterozygote, or the T/C of a de-novo CpG.
const TRANSITION: usize = 160;

/// Haplotype one carries a methylated `C` in a CpG at `TRANSITION`, so its
/// top-strand reads show `T` there — the allele haplotype two carries. Those
/// reads say nothing about the haplotype and must cast no vote: the `C`
/// haplotype, its CpG profiled at the strand's own conversion level, explains
/// them as well as the `T` one does. Scored against a `C` at beta 1/2 instead,
/// each would vote `T` at a cost of 3, so `converted` of them against the
/// eight bottom-strand reads at 26 flip the link once there are more than
/// seventy.
fn methylated_cpg_het(scenario: Scenario, converted: usize) -> Scenario {
    scenario
        .snv(SNV, Hap::One)
        .methylated(Hap::One)
        .pairs(Hap::One, Strand::Top, LEFT, converted)
        .pairs(Hap::One, Strand::Bottom, LEFT, 4)
        .pairs(Hap::Two, Strand::Bottom, LEFT, 4)
        // Balances the SNV without reaching the CpG.
        .pairs(Hap::Two, Strand::Top, 40, converted)
}

/// Decisive: with the CpG profile off (`with_profiled_cpgs`) the eighty
/// converted reads flip the link into a switch error.
#[test]
fn converted_reads_at_a_reference_cpg_het_are_not_evidence() -> Result<()> {
    let scenario = methylated_cpg_het(
        Scenario::new(1000).reference_at(TRANSITION, b"CG").het(TRANSITION, b'T', Hap::Two),
        80,
    );
    let vcf = call_scenario(&scenario, &[])?;

    let phasing = scenario.check_phasing(&vcf)?;
    assert_eq!(phasing, Phasing { blocks: blocks([&[SNV, TRANSITION]]), unphased: vec![] });
    Ok(())
}

/// The `C` alt before a reference `G` makes a CpG on haplotype one only, and
/// its methylated top-strand reads show the reference `T`. Those reads count
/// against the `C`'s allele fraction, so the thresholds allow only a few, too
/// few to flip the link at a cost of 3 each: this one passes with the CpG
/// profile off too, and guards only that the site is phased at all.
#[test]
fn converted_reads_at_a_denovo_cpg_het_are_not_evidence() -> Result<()> {
    let scenario = methylated_cpg_het(
        Scenario::new(1000).reference_at(TRANSITION, b"TG").het(TRANSITION, b'C', Hap::One),
        8,
    );
    let vcf = call_scenario(&scenario, &[])?;

    let phasing = scenario.check_phasing(&vcf)?;
    assert_eq!(phasing, Phasing { blocks: blocks([&[SNV, TRANSITION]]), unphased: vec![] });
    Ok(())
}

/// Outside a CpG a top-strand `T` is the `T` allele. Here it is the only link:
/// haplotype one's reads cover one site or the other, never both, and so do
/// haplotype two's bottom-strand reads.
#[test]
fn converted_reads_outside_a_cpg_are_evidence() -> Result<()> {
    let scenario = Scenario::new(1000)
        .reference_at(TRANSITION, b"CA")
        .snv(SNV, Hap::One)
        .het(TRANSITION, b'T', Hap::Two)
        .pairs(Hap::Two, Strand::Top, LEFT, 4)
        .pairs(Hap::Two, Strand::Bottom, 40, 4)
        .pairs(Hap::Two, Strand::Bottom, 140, 4)
        .pairs(Hap::One, Strand::Top, 40, 4)
        .pairs(Hap::One, Strand::Bottom, 40, 4)
        .pairs(Hap::One, Strand::Top, 140, 4)
        .pairs(Hap::One, Strand::Bottom, 140, 4);
    let vcf = call_scenario(&scenario, &[])?;

    let phasing = scenario.check_phasing(&vcf)?;
    assert_eq!(phasing, Phasing { blocks: blocks([&[SNV, TRANSITION]]), unphased: vec![] });
    Ok(())
}

/// Only the indel links the two SNVs. Run once with an insertion and once with
/// a deletion, each placed so it has one alignment.
#[test]
fn an_indel_phased_by_the_pair_hmm_bridges_two_snvs() -> Result<()> {
    // Too far apart for any template to reach both SNVs; templates from each
    // side reach the indel.
    const LEFT_SNV: usize = 120;
    const INDEL: usize = 370;
    const RIGHT_SNV: usize = 620;
    let around = |indel: fn(Scenario) -> Scenario| {
        indel(Scenario::new(1000))
            .snv(LEFT_SNV, Hap::One)
            .snv(RIGHT_SNV, Hap::One)
            .balanced(100, 350, 5)
            .balanced(340, 590, 5)
    };
    let insertion = around(|s| s.reference_at(INDEL, b"AT").insertion(INDEL, b"GG", Hap::Two));
    let deletion = around(|s| s.reference_at(INDEL, b"AGTC").deletion(INDEL, 2, Hap::Two));

    for scenario in [insertion, deletion] {
        let vcf = call_scenario(&scenario, &[])?;
        let phasing = scenario.check_phasing(&vcf)?;
        assert_eq!(
            phasing,
            Phasing { blocks: blocks([&[LEFT_SNV, INDEL, RIGHT_SNV]]), unphased: vec![] },
            "{vcf}"
        );
    }
    Ok(())
}

/// The oracle itself: a correctly phased file with one genotype swapped is a
/// switch error, and the same file with the whole block swapped no longer
/// opens on `0|1`.
#[test]
fn check_phasing_rejects_a_switch_and_an_unnormalised_block() -> Result<()> {
    let scenario = two_sites_in_trans();
    let vcf = call_scenario(&scenario, &[])?;
    scenario.check_phasing(&vcf)?;

    let swap = |vcf: &str, pos: usize| -> String {
        let pos = (pos + 1).to_string();
        vcf.lines()
            .map(|line| {
                if line.split('\t').nth(1) == Some(&pos) {
                    line.replace("0|1", "x").replace("1|0", "0|1").replace('x', "1|0")
                } else {
                    line.to_owned()
                }
            })
            .collect::<Vec<_>>()
            .join("\n")
    };
    let switched = swap(&vcf, SITE_B);
    assert!(scenario.check_phasing(&switched).is_err(), "a switch error must fail");
    let both = swap(&switched, SITE_A);
    assert!(scenario.check_phasing(&both).is_err(), "a block opening on 1|0 must fail");
    Ok(())
}

/// CpG-only output drops the variant sites phasing links, so the run neither
/// phases nor calls indels, and is not refused for it.
#[test]
fn cpg_only_output_is_not_phased() -> Result<()> {
    let output = rastair().args(CALL_TEST_BAM).args([CHR19_SMALL, "--cpgs-only"]).output()?;
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    Ok(())
}

/// How phasing on real reads scores against a phased truth set.
#[derive(Debug, Default)]
struct Accuracy {
    /// Phased heterozygous SNVs inside the confident regions.
    phased: usize,
    /// Of those, the ones the truth does not have as the same het.
    phased_wrongly: usize,
    /// Adjacent pairs of truly heterozygous sites inside one block.
    pairs: usize,
    /// Of those, the pairs whose relative phase disagrees with the truth.
    switches: usize,
}

/// The columns of a VCF data line accuracy is scored on.
struct Genotyped<'a> {
    /// As written, 1-based.
    pos: usize,
    alleles: Vec<&'a str>,
    gt: &'a str,
    ps: Option<&'a str>,
    filter: &'a str,
}

impl<'a> Genotyped<'a> {
    fn parse(line: &'a str) -> Option<Self> {
        let columns: Vec<&str> = line.split('\t').collect();
        let mut alleles = vec![*columns.get(3)?];
        alleles.extend(columns.get(4)?.split(','));
        let keys: Vec<&str> = columns.get(8)?.split(':').collect();
        let values: Vec<&str> = columns.get(9)?.split(':').collect();
        let value = |key: &str| values.get(keys.iter().position(|k| *k == key)?).copied();
        Some(Self {
            pos: columns.get(1)?.parse().ok()?,
            alleles,
            gt: value("GT")?,
            ps: value("PS"),
            filter: columns.get(6)?,
        })
    }

    /// The two allele strings a phased GT puts on haplotypes 1 and 2, when
    /// they differ.
    fn haplotypes(&self) -> Option<[&'a str; 2]> {
        let (one, two) = self.gt.split_once('|')?;
        let allele = |index: &str| self.alleles.get(index.parse::<usize>().ok()?).copied();
        let pair = [allele(one)?, allele(two)?];
        (pair[0] != pair[1]).then_some(pair)
    }
}

fn score_against_truth(vcf: &str, truth: &str, confident: &str) -> Result<Accuracy> {
    let intervals: Vec<(usize, usize)> = confident
        .lines()
        .filter_map(|line| {
            let mut fields = line.split('\t').skip(1);
            Some((fields.next()?.parse().ok()?, fields.next()?.parse().ok()?))
        })
        .collect();
    let is_confident = |pos: usize| intervals.iter().any(|&(start, end)| start < pos && pos <= end);
    let truth: BTreeMap<usize, [&str; 2]> = truth
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(Genotyped::parse)
        .filter_map(|record| Some((record.pos, record.haplotypes()?)))
        .collect();

    let mut accuracy = Accuracy::default();
    let mut blocks: BTreeMap<&str, Vec<bool>> = BTreeMap::new();
    for record in vcf.lines().filter(|line| !line.starts_with('#')).filter_map(Genotyped::parse) {
        let snv = record.alleles.iter().all(|allele| allele.len() == 1);
        let Some(called) = record.haplotypes().filter(|_| snv && record.filter == "PASS") else {
            continue;
        };
        if !is_confident(record.pos) {
            continue;
        }
        accuracy.phased += 1;
        let Some(&[first, second]) = truth.get(&record.pos) else {
            accuracy.phased_wrongly += 1;
            continue;
        };
        let same = called == [first, second];
        if !same && called != [second, first] {
            accuracy.phased_wrongly += 1;
            continue;
        }
        let set = record.ps.ok_or_else(|| eyre!("phased without PS at {}", record.pos))?;
        blocks.entry(set).or_default().push(same);
    }
    for orientations in blocks.values() {
        accuracy.pairs += orientations.len().saturating_sub(1);
        accuracy.switches +=
            orientations.iter().zip(orientations.iter().skip(1)).filter(|(a, b)| a != b).count();
    }
    Ok(accuracy)
}

/// 200 kb of NA12878 chr12 TAPS reads (26x) against the Platinum Genomes
/// phased truth, with the bundled model. Floors and ceilings sit just outside
/// what the code measured when they were set; moving them is a decision.
///
/// The reads are not in git: `scripts/make_phasing_slice.py` writes the
/// fixture, and `RASTAIR_PHASING_SLICE` names its directory. Without it the
/// test says so and passes.
#[test]
fn phasing_accuracy_on_an_na12878_slice() -> Result<()> {
    let Some(data) = std::env::var_os("RASTAIR_PHASING_SLICE") else {
        eprintln!("skipping: set RASTAIR_PHASING_SLICE to the make_phasing_slice.py output");
        return Ok(());
    };
    let data = Path::new(&data);
    let temp_dir = TempDir::new()?;
    let vcf = temp_dir.path().join("out.vcf");
    rastair()
        .arg("call")
        .arg("--fasta-file")
        .arg(data.join("ref.fa"))
        .arg(data.join("reads.bam"))
        .args(["--vcf"])
        .arg(&vcf)
        .succeeds()?;

    let accuracy = score_against_truth(
        &std::fs::read_to_string(&vcf)?,
        &std::fs::read_to_string(data.join("truth.vcf"))?,
        &std::fs::read_to_string(data.join("confident.bed"))?,
    )?;
    eprintln!("{accuracy:?}");
    // Measured: phased 320, wrongly 0, pairs 237, switches 0. The scorer agrees
    // with `whatshap compare` on pairs and switches.
    assert!(accuracy.phased >= 320, "fewer sites phased: {accuracy:?}");
    assert!(accuracy.pairs >= 237, "shorter blocks: {accuracy:?}");
    assert!(accuracy.switches == 0, "switch errors: {accuracy:?}");
    assert!(accuracy.phased_wrongly == 0, "non-hets phased: {accuracy:?}");
    Ok(())
}
