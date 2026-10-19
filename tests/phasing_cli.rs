//! `call --phase` end to end, on scenarios whose haplotypes are known: every
//! run is checked against the truth by [`Scenario::check_phasing`], and each
//! test then states which sites had to end up phased together.

mod utils;
use color_eyre::eyre::OptionExt as _;
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
    let vcf = call_scenario(&scenario, &["--phase"])?;

    let phasing = scenario.check_phasing(&vcf)?;
    assert_eq!(phasing, Phasing { blocks: blocks([&[SITE_A, SITE_B]]), unphased: vec![] });
    assert!(vcf.contains("##FORMAT=<ID=PS,"), "PS must be declared");
    assert!(vcf.contains(r#""phase":true"#), "the config names the flag");
    Ok(())
}

#[test]
fn masked_sites_are_called_but_unphased() -> Result<()> {
    let dir = TempDir::new()?;
    let mask = dir.path().join("centromere.txt");
    std::fs::write(
        &mask,
        format!("Chr\tStart\tEnd\tGapType\nscenario\t{SITE_A}\t{}\tcentromere\n", SITE_A + 1),
    )?;
    let scenario = two_sites_in_trans();
    let plain = call_scenario_in(dir.path(), &scenario, &["--phase"])?;
    let masked = call_scenario_in(
        dir.path(),
        &scenario,
        &["--phase", "--phase-mask", mask.to_str().ok_or_eyre("Mask path is not UTF-8")?],
    )?;
    assert_eq!(
        scenario.check_phasing(&masked)?,
        Phasing { blocks: BTreeMap::new(), unphased: vec![SITE_A, SITE_B] }
    );
    assert_nonphase_fields_equal(&plain, &masked);
    Ok(())
}

#[test]
fn a_mask_is_a_barrier_even_when_fragments_bridge_it() -> Result<()> {
    let dir = TempDir::new()?;
    let mask = dir.path().join("blacklist.bed");
    std::fs::write(&mask, "scenario\t250\t251\n")?;
    let scenario = Scenario::new(1000)
        .snv(SITE_A, Hap::One)
        .snv(140, Hap::Two)
        .snv(380, Hap::One)
        .snv(SITE_B, Hap::Two)
        .balanced(LEFT, RIGHT, 5);
    let plain = call_scenario_in(dir.path(), &scenario, &["--phase"])?;
    assert_eq!(scenario.check_phasing(&plain)?.blocks, blocks([&[SITE_A, 140, 380, SITE_B]]));
    let masked = call_scenario_in(
        dir.path(),
        &scenario,
        &["--phase", "--phase-mask", mask.to_str().ok_or_eyre("Mask path is not UTF-8")?],
    )?;
    assert_eq!(
        scenario.check_phasing(&masked)?,
        Phasing { blocks: blocks([&[SITE_A, 140], &[380, SITE_B]]), unphased: vec![] }
    );
    assert_nonphase_fields_equal(&plain, &masked);
    Ok(())
}

fn assert_nonphase_fields_equal(a: &str, b: &str) {
    let canonical = |vcf: &str| {
        vcf_content_lines(vcf)
            .map(|line| {
                let fields: Vec<_> = line.split('\t').collect();
                let mut format: BTreeMap<_, _> = fields
                    .get(8)
                    .into_iter()
                    .flat_map(|v| v.split(':'))
                    .zip(fields.get(9).into_iter().flat_map(|v| v.split(':')))
                    .filter(|(key, _)| *key != "PS")
                    .map(|(key, value)| (key.to_owned(), value.to_owned()))
                    .collect();
                if let Some(gt) = format.get_mut("GT") {
                    let mut alleles: Vec<_> = gt.split(['/', '|']).collect();
                    alleles.sort_unstable();
                    *gt = alleles.join("/");
                }
                (fields.iter().take(8).map(|v| (*v).to_owned()).collect::<Vec<_>>(), format)
            })
            .collect::<Vec<_>>()
    };
    assert_eq!(canonical(a), canonical(b));
}

/// A block is solved on the whole segment, overlap included, and then cut to
/// the core. With 300-base segments `HUB` is in the first segment's core and
/// the second's overlap, and `SPOKE_1` and `SPOKE_2`, in the second's core, are
/// linked only through it. So the second segment phases the spokes as one
/// block, and the hub, alone in its core, stays unphased. Phasing after the
/// trim would leave all three unphased.
#[test]
fn a_block_is_solved_across_the_overlap_and_cut_to_the_core() -> Result<()> {
    const HUB: usize = 120;
    const SPOKE_1: usize = 370;
    const SPOKE_2: usize = 480;
    let scenario = Scenario::new(1000)
        .snv(HUB, Hap::One)
        .snv(SPOKE_1, Hap::Two)
        .snv(SPOKE_2, Hap::One)
        .balanced(LEFT, 350, 5)
        .balanced(LEFT, 460, 5);

    let vcf = call_scenario(
        &scenario,
        &["--phase", "--segment-max-length=300", "--segment-overlap=200"],
    )?;

    let phasing = scenario.check_phasing(&vcf)?;
    assert_eq!(phasing, Phasing { blocks: blocks([&[SPOKE_1, SPOKE_2]]), unphased: vec![HUB] });
    Ok(())
}

/// The SNV the methylation-confounded site is linked to.
const SNV: usize = 130;
/// A C/T heterozygote, or the T/C of a de-novo CpG.
const TRANSITION: usize = 160;

/// Haplotype one carries a methylated `C` in a CpG at `TRANSITION`, so its
/// top-strand reads show `T` there — the allele haplotype two carries. Those
/// reads say nothing about the haplotype and must be dropped. Counted, they
/// would cancel the bottom-strand evidence exactly and leave both sites
/// unphased.
fn methylated_cpg_het(scenario: Scenario) -> Scenario {
    scenario
        .snv(SNV, Hap::One)
        .methylated(Hap::One)
        .pairs(Hap::One, Strand::Top, LEFT, 8)
        .pairs(Hap::One, Strand::Bottom, LEFT, 4)
        .pairs(Hap::Two, Strand::Bottom, LEFT, 4)
        // Balances the SNV without reaching the CpG.
        .pairs(Hap::Two, Strand::Top, 40, 8)
}

#[test]
fn converted_reads_at_a_reference_cpg_het_are_not_evidence() -> Result<()> {
    let scenario = methylated_cpg_het(Scenario::new(1000).reference_at(TRANSITION, b"CG").het(
        TRANSITION,
        b'T',
        Hap::Two,
    ));
    let vcf = call_scenario(&scenario, &["--phase"])?;

    let phasing = scenario.check_phasing(&vcf)?;
    assert_eq!(phasing, Phasing { blocks: blocks([&[SNV, TRANSITION]]), unphased: vec![] });
    Ok(())
}

/// The `C` alt before a reference `G` makes a CpG on haplotype one only, and
/// its methylated top-strand reads show the reference `T`.
#[test]
fn converted_reads_at_a_denovo_cpg_het_are_not_evidence() -> Result<()> {
    let scenario = methylated_cpg_het(Scenario::new(1000).reference_at(TRANSITION, b"TG").het(
        TRANSITION,
        b'C',
        Hap::One,
    ));
    let vcf = call_scenario(&scenario, &["--phase"])?;

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
    let vcf = call_scenario(&scenario, &["--phase"])?;

    let phasing = scenario.check_phasing(&vcf)?;
    assert_eq!(phasing, Phasing { blocks: blocks([&[SNV, TRANSITION]]), unphased: vec![] });
    Ok(())
}

/// Only the indel links the two SNVs. Indel records are not phased, so the two
/// stay apart with and without indel calling, and the called indel is written
/// unphased. Run once with an insertion and once with a deletion, each placed
/// so it has one alignment.
#[test]
fn an_indel_called_by_the_pair_hmm_is_not_phased_and_does_not_bridge_snvs() -> Result<()> {
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
        let vcf = call_scenario(&scenario, &["--phase"])?;
        let phasing = scenario.check_phasing(&vcf)?;
        assert_eq!(
            phasing,
            Phasing { blocks: BTreeMap::new(), unphased: vec![LEFT_SNV, RIGHT_SNV] },
            "{vcf}"
        );

        let vcf = call_scenario(&scenario, &["--phase", "--experimental-indels"])?;
        let phasing = scenario.check_phasing(&vcf)?;
        assert_eq!(
            phasing,
            Phasing { blocks: BTreeMap::new(), unphased: vec![LEFT_SNV, INDEL, RIGHT_SNV] },
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
    let vcf = call_scenario(&scenario, &["--phase"])?;
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

/// Without `--phase` no record is phased and the config does not mention it.
/// The `PS` header line is declared either way.
#[test]
fn an_unphased_run_mentions_no_phasing() -> Result<()> {
    let scenario = two_sites_in_trans();
    let vcf = call_scenario(&scenario, &[])?;

    let phasing = scenario.check_phasing(&vcf)?;
    assert_eq!(phasing, Phasing { blocks: BTreeMap::new(), unphased: vec![SITE_A, SITE_B] });
    assert!(!vcf.contains(r#""phase""#), "no phasing config");
    Ok(())
}

/// CpG-only output drops the variant sites phasing links, so the run is
/// refused rather than paying for observations it throws away.
#[test]
fn phasing_refuses_cpg_only_output() -> Result<()> {
    let output = rastair()
        .args(CALL_TEST_BAM)
        .args([CHR19_SMALL, NO_ML, "--phase", "--cpgs-only"])
        .output()?;
    assert!(!output.status.success(), "--phase --cpgs-only must fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("--phase has nothing to work with under CpG-only output"), "{stderr}");
    Ok(())
}

#[test]
fn pair_hmm_expansion_phases_a_chain_that_requires_pruning() -> Result<()> {
    let positions: Vec<_> = (0..10).map(|index| 120 + index * 150).collect();
    let mut scenario = Scenario::new(2200);
    for (index, &pos) in positions.iter().enumerate() {
        scenario = scenario.snv(pos, if index % 2 == 0 { Hap::One } else { Hap::Two });
    }
    for pair in positions.windows(2) {
        let [left, right] = pair else { continue };
        scenario = scenario.balanced(left - 20, right - 20, 3);
    }
    let vcf = call_scenario(&scenario, &["--phase"])?;
    let phasing = scenario.check_phasing(&vcf)?;
    assert_eq!(phasing, Phasing { blocks: blocks([positions.as_slice()]), unphased: vec![] });
    Ok(())
}

#[test]
fn fetching_the_same_fragment_in_byte_budget_subdivisions_does_not_change_phase() -> Result<()> {
    let scenario = two_sites_in_trans();
    let ordinary = scenario.check_phasing(&call_scenario(&scenario, &["--phase"])?)?;
    let subdivided = scenario
        .check_phasing(&call_scenario(&scenario, &["--phase", "--segment-max-bytes=1"])?)?;
    assert_eq!(ordinary, subdivided);
    Ok(())
}
