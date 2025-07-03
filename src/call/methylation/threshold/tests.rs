// todo: test this!
// - create builder for records, maybe not here
// - test: only consider CpG sites
// - test: only consider T base in alts
// - test: check VAF threshold
// - test: check read depth threshold
// - test: check A and G bases vs C and T bases

use super::*;
use crate::call::{test_helpers::variant_pileup, variant_calling::VariantCallingParams};
use color_eyre::Result;
use insta::assert_debug_snapshot;

#[test]
fn cpg_c_methylated() -> Result<()> {
    let pileup = variant_pileup("chr19", 6105084)?;
    let config = ThresholdParams::default();
    let mut metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
    call(&config, &mut metrics, None, None)?;

    assert_debug_snapshot!((
        pileup.chrom(),
        pileup.pos,
        pileup.reference_base,
        &metrics.info.in_cp_g,
        &metrics.info.allele_specific_strand_bias,
        &metrics.samples[0].methylated,
    ), @r#"
    (
        "chr19",
        6105084,
        C,
        CpG::C,
        AlleleSpecificStrandBias(
            [
                ByStrand {
                    base: C,
                    ot: 4,
                    ob: 19,
                },
                ByStrand {
                    base: T,
                    ot: 14,
                    ob: 0,
                },
            ],
        ),
        Methylated::OriginalCpG(
            0.7777777777777778,
        ),
    )
    "#);
    Ok(())
}

#[test]
fn cpg_g_methylated() -> Result<()> {
    let pileup = variant_pileup("chr19", 6105085)?;
    let config = ThresholdParams::default();
    let mut metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
    call(&config, &mut metrics, None, None)?;

    assert_debug_snapshot!((
        pileup.chrom(),
        pileup.pos,
        pileup.reference_base,
        &metrics.info.in_cp_g,
        &metrics.info.allele_specific_strand_bias,
        &metrics.samples[0].methylated,
    ), @r#"
    (
        "chr19",
        6105085,
        G,
        CpG::G,
        AlleleSpecificStrandBias(
            [
                ByStrand {
                    base: G,
                    ot: 18,
                    ob: 4,
                },
                ByStrand {
                    base: A,
                    ot: 0,
                    ob: 15,
                },
            ],
        ),
        Methylated::OriginalCpG(
            0.7894736842105263,
        ),
    )
    "#);
    Ok(())
}

#[test]
fn c_but_not_cpg() -> Result<()> {
    let pileup = variant_pileup("chr19", 6105197)?;
    let config = ThresholdParams::default();
    let mut metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
    call(&config, &mut metrics, None, None)?;

    assert_debug_snapshot!((
        pileup.chrom(),
        pileup.pos,
        pileup.reference_base,
        &metrics.info.in_cp_g,
        &metrics.info.allele_specific_strand_bias,
        &metrics.samples[0].methylated,
    ), @r#"
    (
        "chr19",
        6105197,
        C,
        NoCpg,
        AlleleSpecificStrandBias(
            [
                ByStrand {
                    base: C,
                    ot: 20,
                    ob: 13,
                },
                ByStrand {
                    base: T,
                    ot: 1,
                    ob: 0,
                },
            ],
        ),
        Methylated::Unknown,
    )
    "#);
    Ok(())
}

#[test]
fn random_other_variant() -> Result<()> {
    let pileup = variant_pileup("chr19", 6105114)?;
    let config = ThresholdParams::default();
    let mut metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
    call(&config, &mut metrics, None, None)?;

    assert_debug_snapshot!((
        pileup.chrom(),
        pileup.pos,
        pileup.reference_base,
        &metrics.info.in_cp_g,
        &metrics.info.allele_specific_strand_bias,
        &metrics.samples[0].methylated,
    ), @r#"
    (
        "chr19",
        6105114,
        A,
        NoCpg,
        AlleleSpecificStrandBias(
            [
                ByStrand {
                    base: A,
                    ot: 20,
                    ob: 18,
                },
                ByStrand {
                    base: G,
                    ot: 0,
                    ob: 1,
                },
            ],
        ),
        Methylated::Unknown,
    )
    "#);
    Ok(())
}

#[test]
fn methylatable_position_not_methylated() -> Result<()> {
    let pileup = variant_pileup("chr19", 6115809)?;
    let config = ThresholdParams::default();
    let mut metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
    call(&config, &mut metrics, None, None)?;

    assert_debug_snapshot!((
        pileup.chrom(),
        pileup.pos,
        pileup.reference_base,
        &metrics.info.in_cp_g,
        &metrics.info.allele_specific_strand_bias,
        &metrics.samples[0].methylated,
    ), @r#"
    (
        "chr19",
        6115809,
        C,
        CpG::C,
        AlleleSpecificStrandBias(
            [
                ByStrand {
                    base: C,
                    ot: 5,
                    ob: 4,
                },
                ByStrand {
                    base: A,
                    ot: 1,
                    ob: 0,
                },
            ],
        ),
        Methylated::NoEvidence,
    )
    "#);

    Ok(())
}
