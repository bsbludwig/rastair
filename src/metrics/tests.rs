use crate::{
    call::{ml::MachineLearningParams, process},
    metrics::PileupMetrics,
    sequence::{ReaderParams, Segment},
    utils::{PileupMetricsIteratorExt, default},
};
use color_eyre::{
    Result,
    eyre::{Context as _, ContextCompat as _},
};
use seqair_types::Probability;

#[cfg(not(feature = "experimental-seqair"))]
fn build_test_pileups(
    pileups: impl Iterator<Item = crate::call::pileup::Pileup>,
    segment: &Segment,
    threshold_filters: &process::ThresholdFilterParams,
) -> Vec<PileupMetrics> {
    process::calculate_pileup_metrics(pileups, segment)
        .map(|x| {
            let mut x = x.unwrap();
            process::apply_threshold_filters(&mut x, threshold_filters).unwrap();
            x
        })
        .collect()
}

#[cfg(feature = "experimental-seqair")]
fn build_test_pileups(
    pileups: impl Iterator<Item = PileupMetrics>,
    _segment: &Segment,
    threshold_filters: &process::ThresholdFilterParams,
) -> Vec<PileupMetrics> {
    pileups
        .map(|mut x| {
            process::apply_threshold_filters(&mut x, threshold_filters).unwrap();
            x
        })
        .collect()
}

#[test]
fn test_cpg_detection() -> Result<()> {
    let ml_threshold = Some(Probability::new_panicky(0.5));
    let mut readers = ReaderParams::test_data().around("chr19", 6105711).pileup_readers()?;
    let chunk = readers.segments(1000, 0)?.next().wrap_err("failed to fetch segment")?;

    let pileup_mapping_params = process::PileupMappingParams::default();
    let (segment, pileups) = process::get_pileups(&mut readers, &chunk, &pileup_mapping_params)
        .wrap_err("failed to process region")?;

    let threshold_filters = process::ThresholdFilterParams {
        variant_calling: default(),
        methylation: default(),
        denovo_cpg: default(),
    };

    let pileups = build_test_pileups(pileups, &segment, &threshold_filters)
        .into_iter()
        .map_surrounding(|b, c, a| process::propagate_denovo_pass_flags(b, c, a, ml_threshold))
        .collect::<Result<Vec<_>>>()?;

    // get pileups for a CpG site
    let ref_c = pileups.iter().find(|p| p.pos == 6105711).wrap_err("Could not find C pileup")?;
    let ref_g = pileups.iter().find(|p| p.pos == 6105712).wrap_err("Could not find G pileup")?;

    assert!(*ref_c.pos_metrics.cpg);
    assert!(*ref_g.pos_metrics.cpg);

    Ok(())
}

#[test]
fn set_filters() -> Result<()> {
    let ml_threshold = Some(Probability::new_panicky(0.5));
    let mut readers = ReaderParams::test_data().around("chr19", 6105742).pileup_readers()?;
    let chunk = readers.segments(1000, 0)?.next().wrap_err("failed to fetch segment")?;

    let pileup_mapping_params = process::PileupMappingParams::default();
    let (segment, pileups) = process::get_pileups(&mut readers, &chunk, &pileup_mapping_params)
        .wrap_err("failed to process region")?;

    let threshold_filters = process::ThresholdFilterParams {
        variant_calling: default(),
        methylation: default(),
        denovo_cpg: default(),
    };

    let ml = MachineLearningParams::default().init()?;

    let pileups = build_test_pileups(pileups, &segment, &threshold_filters)
        .into_iter()
        .map_surrounding(|b, c, a| process::propagate_denovo_pass_flags(b, c, a, ml_threshold))
        .collect::<Result<Vec<_>>>()?;

    let pileups = pileups
        .into_iter()
        .map_surrounding(|b, c, a| process::add_ml_metrics(b, c, a, &ml, true))
        .collect::<Result<Vec<_>>>()?;

    // get pileups for a CpG site
    let low_dp_on_a =
        pileups.iter().find(|p| p.pos == 6105742).wrap_err("Could not find C pileup")?;

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
