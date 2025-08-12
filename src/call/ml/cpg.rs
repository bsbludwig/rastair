use super::utils::*;
use crate::{
    utils::Base,
    vcf::{ByStrand, Record, utils::NoStrandBiasForBaseErrorExt},
};
use ndarray::{Array1, Array2, array};
use tracing::{debug, instrument};

/// Extract feature parameters from a VCF record for CpG classification
#[instrument(level = "debug", skip_all, fields(pos = record.main.pos))]
pub fn params_from_record(
    record: &Record,
    before: Option<&Record>,
    after: Option<&Record>,
    target_alt: Base,
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

    // One-hot encode ref and alt (fetch specific alt for CpG methylation)
    let (ref_a, ref_c, ref_g, ref_t) = one_hot_encode_base(Some(ref_base));
    let (alt_a, alt_c, alt_g, alt_t) = one_hot_encode_base(Some(target_alt));

    // Extract normalized allele depths
    let ad_ref = record.info.allele_read_depth.first().copied().unwrap_or(0) as f64;
    let alt_index = record.main.alt.iter().position(|a| a == &target_alt).unwrap_or(0);
    let ad_alt = record.info.allele_read_depth.get(alt_index + 1).copied().unwrap_or(0) as f64;

    // Extract normalized strand bias counts
    let ref_strand = record.strand_count(ref_base).or_empty();
    let alt_strand = record.strand_count(target_alt).or_empty();

    let sb_ot_ref = f64::from(ref_strand.ot);
    let sb_ob_ref = f64::from(ref_strand.ob);
    let sb_ot_alt = f64::from(alt_strand.ot);
    let sb_ob_alt = f64::from(alt_strand.ob);

    // Calculate alt_score based on ref base (C vs G)
    let alt_score = if ref_base == Base::C {
        // For C: use "ob" (original bottom) strand data
        let bq_ob_alt = get_strand_base_quality(record, target_alt).ob;
        let bq_ob_ref = get_strand_base_quality(record, ref_base).ob;
        (sb_ob_alt * bq_ob_alt + 1.0) / (sb_ob_ref * bq_ob_ref + 1.0).log2()
    } else {
        // For G: use "ot" (original top) strand data
        let bq_ot_alt = get_strand_base_quality(record, target_alt).ot;
        let bq_ot_ref = get_strand_base_quality(record, ref_base).ot;
        (sb_ot_alt * bq_ot_alt + 1.0) / (sb_ot_ref * bq_ot_ref + 1.0).log2()
    };

    // Extract base quality metrics
    let bq_ref = record.info.allele_base_quality.first().copied().unwrap_or(0.0);
    let bq_alt = record.info.allele_base_quality.get(alt_index + 1).copied().unwrap_or(0.0);
    let bq_ot_ref = get_strand_base_quality(record, ref_base).ot;
    let bq_ob_ref = get_strand_base_quality(record, ref_base).ob;
    let bq_ot_alt = get_strand_base_quality(record, target_alt).ot;
    let bq_ob_alt = get_strand_base_quality(record, target_alt).ob;

    // Extract mapping quality metrics
    let mq_ref = record.info.allele_map_quality.first().copied().unwrap_or(0.0);
    let mq_alt = record.info.allele_map_quality.get(alt_index + 1).copied().unwrap_or(0.0);
    let mq_ot_ref = get_strand_map_quality(record, ref_base).ot;
    let mq_ob_ref = get_strand_map_quality(record, ref_base).ob;
    let mq_ot_alt = get_strand_map_quality(record, target_alt).ot;
    let mq_ob_alt = get_strand_map_quality(record, target_alt).ob;

    // Extract other metrics
    let pos_in_read_ref = record.info.position_in_read.first().copied().unwrap_or(0.0);
    let pos_in_read_alt = record.info.position_in_read.get(alt_index + 1).copied().unwrap_or(0.0);
    let num_aligned_bases_ref = record.info.num_aligned_bases.first().copied().unwrap_or(0.0);
    let num_aligned_bases_alt =
        record.info.num_aligned_bases.get(alt_index + 1).copied().unwrap_or(0.0);
    let num_indels_ref = record.info.num_indels.first().copied().unwrap_or(0.0);
    let num_indels_alt = record.info.num_indels.get(alt_index + 1).copied().unwrap_or(0.0);

    // Calculate adjacent position features
    let (beta_ratio, ad_alt_adj, alt_score_adj) =
        calculate_adjacent_features(record, before, after);

    array![[
        ad_alt_adj,
        alt_score_adj,
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
        bq_ot_ref,
        bq_ob_ref,
        bq_ot_alt,
        bq_ob_alt,
        mq_ot_ref,
        mq_ob_ref,
        mq_ot_alt,
        mq_ob_alt,
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
fn calculate_adjacent_features(
    record: &Record,
    before: Option<&Record>,
    after: Option<&Record>,
) -> (f64, f64, f64) {
    let ref_base = &record.main.r#ref;

    if ref_base == "C"
        && let Some(after) = after
        && after.main.r#ref == "G"
    {
        // Calculate beta of current position
        let beta_center = f64::from(record.strand_count(Base::T).or_empty().ot)
            / (f64::from(record.strand_count(Base::T).or_empty().ot)
                + f64::from(record.strand_count(Base::C).or_empty().ot));
        if let Some(alt_index) = after.main.alt.iter().position(|a| a == "A") {
            // For C positions: look for G→A transitions in the after record
            let ad_alt =
                after.info.allele_read_depth.get(alt_index + 1).copied().unwrap_or(0) as f64;
            let depth = *after.info.read_depth as f64;
            let ad_alt_norm = ad_alt / depth;

            // Calculate alt_score for G→A
            let alt_strand = after.strand_count(Base::A).or_empty();
            let ref_strand = after.strand_count(Base::G).or_empty();
            let bq_ot_alt = get_strand_base_quality(after, Base::A).ot;
            let bq_ot_ref = get_strand_base_quality(after, Base::G).ot;
            let alt_score = (f64::from(alt_strand.ot) * bq_ot_alt + 1.0)
                / (f64::from(ref_strand.ot) * bq_ot_ref + 1.0);
            let beta_after =
                f64::from(alt_strand.ob) / (f64::from(alt_strand.ob) + f64::from(ref_strand.ob));
            let beta_ratio = ((beta_center + 1.0) / (beta_after + 1.0)).log2();
            (beta_ratio, ad_alt_norm, alt_score.log2())
        } else {
            (((beta_center + 1.0) / 1.0).log2(), 0.0, 0.0)
        }
    } else if ref_base == "G"
        && let Some(before) = before
        && before.main.r#ref == "C"
    {
        let beta_center = f64::from(record.strand_count(Base::A).or_empty().ob)
            / (f64::from(record.strand_count(Base::A).or_empty().ob)
                + f64::from(record.strand_count(Base::G).or_empty().ob));
        if let Some(alt_index) = before.main.alt.iter().position(|a| a == "T") {
            // For G positions: look for C→T transitions in the before record
            let ad_alt =
                before.info.allele_read_depth.get(alt_index + 1).copied().unwrap_or(0) as f64;
            let depth = *before.info.read_depth as f64;
            let ad_alt_norm = ad_alt / depth;

            // Calculate alt_score for C→T
            let alt_strand = before.strand_count(Base::T).or_empty();
            let ref_strand = before.strand_count(Base::C).or_empty();
            let bq_ob_alt = get_strand_base_quality(before, Base::T).ob;
            let bq_ob_ref = get_strand_base_quality(before, Base::C).ob;
            let alt_score = (f64::from(alt_strand.ob) * bq_ob_alt + 1.0)
                / (f64::from(ref_strand.ob) * bq_ob_ref + 1.0);
            let beta_before =
                f64::from(alt_strand.ot) / (f64::from(alt_strand.ot) + f64::from(ref_strand.ot));
            let beta_ratio = ((beta_center + 1.0) / (beta_before + 1.0)).log2();
            (beta_ratio, ad_alt_norm, alt_score.log2())
        } else {
            (((beta_center + 1.0) / 1.0).log2(), 0.0, 0.0)
        }
    } else {
        // No adjacent evidence for methylation, return defaults
        debug!(%ref_base, before=%before.is_some(), after=%after.is_some(), "No adjacent evidence for methylation");
        (0.0, 0.0, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        call::{
            process::PileupMappingParams, test_helpers::variant_pileup,
            variant_calling::VariantCallingParams,
        },
        sequence::ReaderParams,
    };
    use color_eyre::Result;

    #[test]
    #[ignore = "needs big test file"]
    fn ch12_10588_c_t() -> Result<()> {
        let reader = ReaderParams::test_with(
            "tmp/taps/NA12878_aa_chr12.bam",
            "tmp/na12878/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        );
        {
            let record =
                reader.pileup("chr12", 10587)?.variant_metrics(&VariantCallingParams::default())?;
            let fields = params_from_record(&record, None, None, Base::T);
            eprintln!(
                "{}:{}_{}\t{}",
                record.main.chrom,
                record.main.pos + 1,
                record.main.r#ref,
                to_tsv(fields)
            );
        }
        {
            let record =
                reader.pileup("chr12", 10601)?.variant_metrics(&VariantCallingParams::default())?;
            let fields = params_from_record(&record, None, None, Base::A);
            eprintln!(
                "{}:{}_{}\t{}",
                record.main.chrom,
                record.main.pos,
                record.main.r#ref,
                to_tsv(fields)
            );
        }
        Ok(())
    }

    fn to_tsv(fields: Array2<f64>) -> String {
        fields.iter().map(|f| f.to_string()).collect::<Vec<_>>().join("\t")
    }
}
