use crate::{
    call::{ml::MachineLearningParams, process},
    sequence::ReaderParams,
    utils::default,
};
use color_eyre::{
    Result,
    eyre::{Context as _, ContextCompat as _},
};
use rastair_types::Probability;

#[test]
fn test_cpg_detection() -> Result<()> {
    let ml_threshold = Some(Probability::new_panicky(0.8));
    let mut readers = ReaderParams::test_data().around("chr19", 6105711).readers()?;
    let chunk = readers.segments(1000, 0)?.next().wrap_err("failed to fetch segment")?;

    let pileup_mapping_params = process::PileupMappingParams { variant_calling: default() };
    let (segment, pileups) = process::get_pileups(&mut readers, &chunk, &pileup_mapping_params)
        .wrap_err("failed to process region")?;
    let mut pileups = process::calculate_pileup_metrics(
        pileups,
        &segment,
        &process::PileupMetricsParams { variant_calling: default(), methylation: default() },
    )
    .collect::<Result<Vec<_>>>()
    .wrap_err("failed to calculate pileup metrics")?;

    process::apply_threshold_filters(
        &mut pileups,
        &process::ThresholdFilterParams {
            variant_calling: default(),
            methylation: default(),
            denovo_cpg: default(),
        },
    )
    .wrap_err("Failed to apply threshold filters")?;

    process::propagate_cpg_pass_flags(&mut pileups, ml_threshold)
        .wrap_err("Failed to propagate CpG pass flags")?;

    // get pileups for a CpG site
    let ref_c =
        pileups.iter().find(|p| p.pileup.pos == 6105711).wrap_err("Could not find C pileup")?;
    let ref_g =
        pileups.iter().find(|p| p.pileup.pos == 6105712).wrap_err("Could not find G pileup")?;

    // dbg!(ref_c.pos_filters.other_pos_in_cpg_passes);
    // dbg!(ref_g.pos_filters.other_pos_in_cpg_passes);

    // dbg!(ref_c.pos_filters.len());
    // ref_c.alts.iter().for_each(|alt| {
    //     dbg!(alt.filters.filters.other_pos_in_cpg_passes);
    //     dbg!(alt.filters.filters.len());
    // });
    // dbg!(ref_g.pos_filters.len());
    // ref_g.alts.iter().for_each(|alt| {
    //     dbg!(alt.filters.filters.other_pos_in_cpg_passes);
    //     dbg!(alt.filters.filters.len());
    // });

    assert!(ref_c.pass(ml_threshold));
    assert!(ref_g.pass(ml_threshold));

    Ok(())
}

#[test]
fn set_filters() -> Result<()> {
    let ml_threshold = Some(Probability::new_panicky(0.8));
    let mut readers = ReaderParams::test_data().around("chr19", 6105742).readers()?;
    let chunk = readers.segments(1000, 0)?.next().wrap_err("failed to fetch segment")?;

    let pileup_mapping_params = process::PileupMappingParams { variant_calling: default() };
    let (segment, pileups) = process::get_pileups(&mut readers, &chunk, &pileup_mapping_params)
        .wrap_err("failed to process region")?;
    let mut pileups = process::calculate_pileup_metrics(
        pileups,
        &segment,
        &process::PileupMetricsParams { variant_calling: default(), methylation: default() },
    )
    .collect::<Result<Vec<_>>>()
    .wrap_err("failed to calculate pileup metrics")?;

    process::apply_threshold_filters(
        &mut pileups,
        &process::ThresholdFilterParams {
            variant_calling: default(),
            methylation: default(),
            denovo_cpg: default(),
        },
    )
    .wrap_err("Failed to apply threshold filters")?;

    process::propagate_cpg_pass_flags(&mut pileups, ml_threshold)
        .wrap_err("Failed to propagate CpG pass flags")?;

    process::add_ml_metrics(&mut pileups, &MachineLearningParams::default().init()?)?;

    // get pileups for a CpG site
    let low_dp_on_a =
        pileups.iter().find(|p| p.pileup.pos == 6105742).wrap_err("Could not find C pileup")?;

    // dbg!(&low_dp_on_a.pos_filters);
    // dbg!(low_dp_on_a.alts.iter().map(|alt| &alt.filters.filters).collect::<Vec<_>>());

    assert!(!low_dp_on_a.pass(ml_threshold));

    Ok(())
}

// #[test]
// fn test_allele_specific_strand_bias_1() -> Result<()> {
//     let pileup = variant_pileup("bacteriophage_lambda_CpG", 2636)?;
//     assert_debug_snapshot!((
//             pileup.contig(),
//             pileup.pos,
//             pileup.reference_base,
//             &pileup.reads,
//             pileup.allele_specific_strand_bias()
//         ), @r#"
//         (
//             "bacteriophage_lambda_CpG",
//             2636,
//             C,
//             [
//                 C OB Q32 MQ60,
//                 C OB Q36 MQ60,
//                 T OT Q36 MQ60,
//                 C OB Q36 MQ60,
//                 T OT Q36 MQ60,
//                 T OT Q36 MQ60,
//                 T OT Q36 MQ60,
//                 T OT Q36 MQ60,
//                 C OB Q36 MQ60,
//             ],
//             AlleleSpecificStrandBias(
//                 [
//                     ByStrand {
//                         base: C,
//                         ot: 0,
//                         ob: 4,
//                     },
//                     ByStrand {
//                         base: T,
//                         ot: 5,
//                         ob: 0,
//                     },
//                 ],
//             ),
//         )
//         "#);
//     Ok(())
// }

// #[test]
// fn test_in_cpg() -> Result<()> {
//     // a CpG site
//     let pileup = variant_pileup("chr19", 6105084)?;
//     assert_debug_snapshot!((pileup.contig(), pileup.reference_base, pileup.pos, pileup.in_cpg()), @r#"
//         (
//             "chr19",
//             C,
//             6105084,
//             CpG::C,
//         )
//         "#);

//     // a C variant, followed by a C
//     let pileup = variant_pileup("chr19", 6104589)?;
//     assert_debug_snapshot!((pileup.contig(), pileup.reference_base, pileup.pos, pileup.in_cpg()), @r#"
//         (
//             "chr19",
//             C,
//             6104589,
//             NoCpg,
//         )
//         "#);

//     // some random variant with base G
//     let pileup = variant_pileup("chr19", 6105116)?;
//     assert_debug_snapshot!((pileup.contig(), pileup.reference_base, pileup.pos, pileup.in_cpg()), @r#"
//         (
//             "chr19",
//             G,
//             6105116,
//             NoCpg,
//         )
//         "#);
//     Ok(())
// }
