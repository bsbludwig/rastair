// use crate::call::test_helpers::variant_pileup;
// use color_eyre::Result;
// use insta::assert_debug_snapshot;

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
