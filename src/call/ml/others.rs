use super::utils::*;
use crate::{
    utils::Base,
    vcf::{ByStrand, DeNovoCpGCandidate, Record, utils::NoStrandBiasForBaseErrorExt},
};
use ndarray::{Array1, Array2, array};
use tracing::{debug, instrument};

/// Extract feature parameters from a VCF record for CpG classification
#[instrument(level = "debug", skip_all, fields(pos = record.main.pos))]
pub fn params_from_record(
    record: &Record,
    before: Option<&Record>,
    after: Option<&Record>,
    alt: Base,
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

    // One-hot encode ref and alt (use the first alt allele)
    let (ref_a, ref_c, ref_g, ref_t) = one_hot_encode_base(Some(ref_base));
    let (alt_a, alt_c, alt_g, alt_t) = one_hot_encode_base(Some(alt));

    // Extract normalized allele depths
    let ad_ref = record.info.allele_read_depth.first().copied().unwrap_or(0) as f64;
    let ad_alt = record.info.allele_read_depth.get(1).copied().unwrap_or(0) as f64;

    // Extract normalized strand bias counts
    let ref_strand = record.strand_count(ref_base).or_empty();
    let alt_strand = record.strand_count(alt).or_empty();

    let sb_ot_ref = f64::from(ref_strand.ot);
    let sb_ob_ref = f64::from(ref_strand.ob);
    let sb_ot_alt = f64::from(alt_strand.ot);
    let sb_ob_alt = f64::from(alt_strand.ob);

    // Calculate strand bias ratios
    let sb_alt = (sb_ot_alt + 1.0) / (sb_ob_alt + 1.0);
    let sb_ref = (sb_ot_ref + 1.0) / (sb_ob_ref + 1.0);

    // Calculate alt_score
    let bq_ref = record.info.allele_base_quality.first().copied().unwrap_or(0.0);
    let bq_alt = record.info.allele_base_quality.get(1).copied().unwrap_or(0.0);
    let alt_score = ((ad_alt * bq_alt + 1.0) / (ad_ref * bq_ref + 1.0)).log2();

    // Extract base quality metrics
    let bq_ot_ref = get_strand_base_quality(record, ref_base).ot;
    let bq_ob_ref = get_strand_base_quality(record, ref_base).ob;
    let bq_ot_alt = get_strand_base_quality(record, alt).ot;
    let bq_ob_alt = get_strand_base_quality(record, alt).ob;

    // Extract mapping quality metrics
    let mq_ref = record.info.allele_map_quality.first().copied().unwrap_or(0.0);
    let mq_alt = record.info.allele_map_quality.get(1).copied().unwrap_or(0.0);
    let mq_ot_ref = get_strand_map_quality(record, ref_base).ot;
    let mq_ob_ref = get_strand_map_quality(record, ref_base).ob;
    let mq_ot_alt = get_strand_map_quality(record, alt).ot;
    let mq_ob_alt = get_strand_map_quality(record, alt).ob;

    // Extract other metrics
    let pos_in_read_ref = record.info.position_in_read.first().copied().unwrap_or(0.0);
    let pos_in_read_alt = record.info.position_in_read.get(1).copied().unwrap_or(0.0);
    let num_aligned_bases_ref = record.info.num_aligned_bases.first().copied().unwrap_or(0.0);
    let num_aligned_bases_alt = record.info.num_aligned_bases.get(1).copied().unwrap_or(0.0);
    let num_indels_ref = record.info.num_indels.first().copied().unwrap_or(0.0);
    let num_indels_alt = record.info.num_indels.get(1).copied().unwrap_or(0.0);

    // Never change the order of these variables, as they are used in the model
    array![[
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
        sb_alt,
        sb_ref,
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
        num_indels_alt
    ]]
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
    fn test_other_snp_extraction() -> Result<()> {
        let reader = ReaderParams::test_with(
            "tmp/taps/NA12878_aa_chr12.bam",
            "tmp/na12878/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        );

        // Test A>G transition at position 10004
        let record =
            reader.pileup("chr12", 10003)?.variant_metrics(&VariantCallingParams::default())?;
        let fields = params_from_record(&record, None, None, Base::G);
        eprintln!(
            "{}:{}_{}\t{}",
            record.main.chrom,
            record.main.pos + 1,
            record.main.r#ref,
            to_tsv(fields)
        );

        Ok(())
    }

    fn to_tsv(fields: Array2<f64>) -> String {
        fields.iter().map(|f| f.to_string()).collect::<Vec<_>>().join("\t")
    }
}
