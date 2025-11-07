use super::utils::*;
use crate::{
    utils::Base,
    vcf::{DeNovoCpGCandidate, Record, utils::NoStrandBiasForBaseErrorExt},
};
use ndarray::{Array1, Array2, array};
use tracing::{debug, instrument};

/// Extract feature parameters from a VCF record for CpG classification
#[instrument(level = "debug", skip_all, fields(pos = record.main.pos))]
pub fn params_from_record(
    record: &Record,
    before: Option<&Record>,
    after: Option<&Record>,
) -> Array2<f64> {
    let ref_base = Base::from(&record.main.r#ref);
    let depth = *record.info.read_depth as f64;

    let mapq = *record.info.mapping_quality;
    let num_mapq0 = *record.info.mapping_quality0 as f64;
    let region_entropy = record.info.entropy[0];

    // Extract sequence context and one-hot encode
    let seq_ctx = &record.info.sequence_context;
    let (p1a, p1c, p1g, p1t) = one_hot_encode_base(seq_ctx.before_2);
    let (p2a, p2c, p2g, p2t) = one_hot_encode_base(seq_ctx.before_1);
    let (p4a, p4c, p4g, p4t) = one_hot_encode_base(seq_ctx.after_1);
    let (p5a, p5c, p5g, p5t) = one_hot_encode_base(seq_ctx.after_2);

    // One-hot encode ref
    let (ref_a, ref_c, ref_g, ref_t) = one_hot_encode_base(ref_base);

    // Use the DeNovoCpGCandidate enum to get denovo CpG information
    let (target_alt_base, alt_index) = match record.info.de_novo_cp_g_candidate {
        DeNovoCpGCandidate::Candidate { alt_base, alt_index, .. } => (alt_base, alt_index),
        // TODO: Maybe handle adjecent positions differently?
        DeNovoCpGCandidate::NotCandidate | DeNovoCpGCandidate::Adjecent { .. } => {
            debug!("Not a denovo CpG candidate");
            // Return default values for non-candidates
            return array![[0.0; 54]];
        }
    };

    // One-hot encode alt allele
    let (alt_a, alt_c, alt_g, alt_t) = one_hot_encode_base(target_alt_base);

    // Extract normalized allele depths
    let ad_ref = record.info.allele_read_depth.first().copied().unwrap_or(0) as f64;
    let ad_alt = record.info.allele_read_depth.get(alt_index + 1).copied().unwrap_or(0) as f64;

    // Extract normalized strand bias counts
    let ref_strand = record.strand_count(ref_base).or_empty();
    let alt_strand = record.strand_count(target_alt_base).or_empty();

    let sb_ot_ref = f64::from(ref_strand.ot);
    let sb_ob_ref = f64::from(ref_strand.ob);
    let sb_ot_alt = f64::from(alt_strand.ot);
    let sb_ob_alt = f64::from(alt_strand.ob);

    // Calculate alt_score based on target alt allele
    let alt_score = if target_alt_base == Base::C {
        // For C alt alleles: use "ob" (original bottom) strand data
        let bq_ob_alt = get_strand_base_quality(record, target_alt_base).ob;
        let bq_ob_ref = get_strand_base_quality(record, ref_base).ob;
        (sb_ob_alt * bq_ob_alt + 1.0).log2() - (sb_ob_ref * bq_ob_ref + 1.0).log2()
    } else {
        // For G alt alleles: use "ot" (original top) strand data
        let bq_ot_alt = get_strand_base_quality(record, target_alt_base).ot;
        let bq_ot_ref = get_strand_base_quality(record, ref_base).ot;
        (sb_ot_alt * bq_ot_alt + 1.0).log2() - (sb_ot_ref * bq_ot_ref + 1.0).log2()
    };

    // Extract base quality metrics
    let bq_ref = record.info.allele_base_quality.first().copied().unwrap_or(0.0);
    let bq_alt = record.info.allele_base_quality.get(alt_index + 1).copied().unwrap_or(0.0);
    let bq_ot_ref = get_strand_base_quality(record, ref_base).ot;
    let bq_ob_ref = get_strand_base_quality(record, ref_base).ob;
    let bq_ot_alt = get_strand_base_quality(record, target_alt_base).ot;
    let bq_ob_alt = get_strand_base_quality(record, target_alt_base).ob;

    // Extract mapping quality metrics
    let mq_ref = record.info.allele_map_quality.first().copied().unwrap_or(0.0);
    let mq_alt = record.info.allele_map_quality.get(alt_index + 1).copied().unwrap_or(0.0);
    let mq_ot_ref = get_strand_map_quality(record, ref_base).ot;
    let mq_ob_ref = get_strand_map_quality(record, ref_base).ob;
    let mq_ot_alt = get_strand_map_quality(record, target_alt_base).ot;
    let mq_ob_alt = get_strand_map_quality(record, target_alt_base).ob;

    // Extract other metrics
    let pos_in_read_ref = record.info.position_in_read.first().copied().unwrap_or(0.0);
    let pos_in_read_alt = record.info.position_in_read.get(alt_index + 1).copied().unwrap_or(0.0);
    let num_aligned_bases_ref = record.info.num_aligned_bases.first().copied().unwrap_or(0.0);
    let num_aligned_bases_alt =
        record.info.num_aligned_bases.get(alt_index + 1).copied().unwrap_or(0.0);
    let num_indels_ref = record.info.num_indels.first().copied().unwrap_or(0.0);
    let num_indels_alt = record.info.num_indels.get(alt_index + 1).copied().unwrap_or(0.0);

    // Calculate adjacent position features specific to denovo CpGs
    let (beta_ratio, ad_alt_adj, alt_score_adj, sb_adj) =
        calculate_denovo_adjacent_features(record, before, after);

    array![[
        ad_alt_adj,
        alt_score_adj,
        sb_adj,
        ref_a,
        ref_c,
        ref_g,
        ref_t,
        alt_a,
        alt_c,
        alt_g,
        alt_t,
        *mapq,
        num_mapq0,
        p1a,
        p1c,
        p1g,
        p1t,
        p2a,
        p2c,
        p2g,
        p2t,
        p4a,
        p4c,
        p4g,
        p4t,
        p5a,
        p5c,
        p5g,
        p5t,
        region_entropy,
        ad_ref / depth,
        ad_alt / depth,
        sb_ot_ref / depth,
        sb_ob_ref / depth,
        sb_ot_alt / depth,
        sb_ob_alt / depth,
        alt_score,
        bq_ref,
        bq_alt,
        *bq_ot_ref,
        *bq_ob_ref,
        *bq_ot_alt,
        *bq_ob_alt,
        *mq_ot_ref,
        *mq_ob_ref,
        *mq_ot_alt,
        *mq_ob_alt,
        mq_ref,
        mq_alt,
        pos_in_read_ref,
        pos_in_read_alt,
        num_aligned_bases_ref,
        num_aligned_bases_alt,
        num_indels_ref,
        num_indels_alt,
        beta_ratio
    ]]
}

#[instrument(level = "debug", skip_all)]
fn calculate_denovo_adjacent_features(
    record: &Record,
    before: Option<&Record>,
    after: Option<&Record>,
) -> (f64, f64, f64, f64) {
    // Use the DeNovoCpGCandidate enum to determine the adjacent position logic
    match record.info.de_novo_cp_g_candidate {
        DeNovoCpGCandidate::Candidate { alt_base: Base::C, .. } => {
            let beta_center = {
                let c_count = record.strand_count(Base::C).or_empty().ot;
                let t_count = record.strand_count(Base::T).or_empty().ot;
                if c_count + t_count == 0 {
                    0.0
                } else {
                    f64::from(t_count) / (f64::from(t_count) + f64::from(c_count))
                }
            };
            // For C alt alleles (creating CpG with next G): look for G→A at position-1
            if let Some(after) = after
                && let Some(alt_index) = after.main.alt.iter().position(|a| a == "A")
            {
                assert!(after.main.r#ref == "G", "De-novo CpG not followed by G"); // this should always be the case!
                let ad_alt =
                    after.info.allele_read_depth.get(alt_index + 1).copied().unwrap_or(0) as f64;
                let depth = *after.info.read_depth as f64;
                let ad_alt_norm = ad_alt / depth;

                // Calculate alt_score for G→A using bottom strand
                let alt_strand = after.strand_count(Base::A).or_empty();
                let ref_strand = after.strand_count(Base::G).or_empty();
                let bq_alt = get_strand_base_quality(after, Base::A);
                let bq_ref = get_strand_base_quality(after, Base::G);
                let alt_score = (f64::from(alt_strand.ot) * bq_alt.ot + 1.0).log2()
                    - (f64::from(ref_strand.ot) * bq_ref.ot + 1.0).log2();

                // Calculate beta ratio: beta at the center vs beta at the adjacent position
                let beta_after = {
                    let g_count = alt_strand.ob; // FIXME: bug! this is flipped!
                    let a_count = ref_strand.ob;
                    if g_count + a_count == 0 {
                        0.0
                    } else {
                        f64::from(a_count) / (f64::from(a_count) + f64::from(g_count))
                    }
                };
                let beta_ratio = (beta_center + 1.0).log2() - (beta_after + 1.0).log2();

                let sb_adj = f64::from(alt_strand.ob + 1) / f64::from(alt_strand.ot + 1);
                (beta_ratio, ad_alt_norm, alt_score, sb_adj)
            } else {
                let beta_ratio = (beta_center + 1.0).log2();
                (beta_ratio, 0.0, 0.0, 0.0)
            }
        }
        DeNovoCpGCandidate::Candidate { alt_base: Base::G, .. } => {
            let beta_center = {
                let g_count = record.strand_count(Base::G).or_empty().ob;
                let a_count = record.strand_count(Base::A).or_empty().ob;
                if g_count + a_count == 0 {
                    0.0
                } else {
                    f64::from(a_count) / (f64::from(a_count) + f64::from(g_count))
                }
            };
            // For G alt alleles (creating CpG with prev C): look for C→T at position-1
            if let Some(before) = before
                && let Some(alt_index) = before.main.alt.iter().position(|a| a == "T")
            {
                assert!(before.main.r#ref == "C", "De-novo CpG not preceded by C`"); // This should always be the case
                let ad_alt =
                    before.info.allele_read_depth.get(alt_index + 1).copied().unwrap_or(0) as f64;
                let depth = *before.info.read_depth as f64;
                let ad_alt_norm = ad_alt / depth;

                // Calculate alt_score for C→T using ot strand
                let alt_strand = before.strand_count(Base::T).or_empty();
                let ref_strand = before.strand_count(Base::C).or_empty();
                let bq_alt = get_strand_base_quality(before, Base::T);
                let bq_ref = get_strand_base_quality(before, Base::C);
                let alt_score = (f64::from(alt_strand.ob) * bq_alt.ob + 1.0).log2()
                    - (f64::from(ref_strand.ob) * bq_ref.ob + 1.0).log2();

                // Calculate beta ratio: beta at the center vs beta at the adjacent position
                let beta_before = {
                    let c_count = ref_strand.ot;
                    let t_count = alt_strand.ot;
                    if c_count + t_count == 0 {
                        0.0
                    } else {
                        f64::from(t_count) / (f64::from(t_count) + f64::from(c_count))
                    }
                };
                let beta_ratio = (beta_center + 1.0).log2() - (beta_before + 1.0).log2();

                let sb_adj = f64::from(alt_strand.ot + 1) / f64::from(alt_strand.ob + 1);
                (beta_ratio, ad_alt_norm, alt_score, sb_adj)
            } else {
                let beta_ratio = (beta_center + 1.0).log2();
                (beta_ratio, 0.0, 0.0, 0.0)
            }
        }
        _ => {
            // Not a denovo CpG candidate or unexpected alt base
            debug!("No denovo CpG context found for adjacent feature calculation");
            (0.0, 0.0, 0.0, 0.0)
        }
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::{
//         call::{
//             process::PileupMappingParams, test_helpers::variant_pileup,
//             variant_calling::VariantCallingParams,
//         },
//         sequence::ReaderParams,
//     };
//     use color_eyre::Result;

//     #[test]
//     #[ignore = "needs big test file"]
//     fn test_denovo_cpg_params() -> Result<()> {
//         let reader = ReaderParams::test_with(
//             "tmp/taps/NA12878_aa_chr12.bam",
//             "tmp/na12878/GRCh38_full_analysis_set_plus_decoy_hla.fa",
//         );

//         // Test positions that match the expected output
//         {
//             let record =
//                 reader.pileup("chr12", 10601)?.variant_metrics(&VariantCallingParams::default())?;
//             let fields = params_from_record(&record, None, None);
//             let chr = record.main.chrom;
//             let pos = record.main.pos + 1; // Convert to 1-based position
//             let ref_base = record.main.r#ref;
//             let alt_base =
//                 record.info.de_novo_cp_g_candidate.alt_base().map(|b| b.as_str()).unwrap_or("?");
//             let tsv = to_tsv(fields);
//             eprintln!("{chr}:{pos}_{ref_base}>{alt_base}\t{tsv}\tREF");
//         }
//         {
//             let record =
//                 reader.pileup("chr12", 10619)?.variant_metrics(&VariantCallingParams::default())?;
//             let fields = params_from_record(&record, None, None);
//             let chr = record.main.chrom;
//             let pos = record.main.pos + 1; // Convert to 1-based position
//             let ref_base = record.main.r#ref;
//             let alt_base =
//                 record.info.de_novo_cp_g_candidate.alt_base().map(|b| b.as_str()).unwrap_or("?");
//             let tsv = to_tsv(fields);
//             eprintln!("{chr}:{pos}_{ref_base}>{alt_base}\t{tsv}\tREF");
//         }

//         Ok(())
//     }

//     fn to_tsv(fields: Array2<f64>) -> String {
//         let values: Vec<String> = fields.row(0).iter().map(|f| f.to_string()).collect();

//         values.join("\t")
//     }
// }
