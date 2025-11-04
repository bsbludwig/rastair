use super::Record;
use crate::{
    metrics2::PileupMetrics,
    utils::conversion::IntoF64 as _,
    vcf::{Filters, Format, Info},
};
use color_eyre::Result;
use rastair_vcf::{
    VcfFixedFields,
    standard_fields::{AlleleReadDepth, BaseQuality},
};
use smallvec::{SmallVec, smallvec_inline};

impl PileupMetrics {
    pub fn to_vcf_record(&self) -> Result<Record> {
        let main = VcfFixedFields {
            chrom: self.contig(),
            pos: self.pos(),
            id: Default::default(),
            r#ref: self.pileup.reference_base.into(),
            alt: self.alts().iter().map(|alt| (*alt).into()).collect::<SmallVec<_, 2>>(),
            qual: Some(0.9),
        };
        let filters = Filters::default();
        let info = Info {
            allele_read_depth: AlleleReadDepth(
                self.ref_alts_metrics().map(|m| m.depth as usize).collect(),
            ),
            base_quality: BaseQuality(self.pos_metrics.baseq),
            read_depth: todo!(),
            mapping_quality: todo!(),
            mapping_quality0: todo!(),
            samples_with_data: todo!(),
            allele_specific_strand_bias: todo!(),
            sequence_context: todo!(),
            allele_frequency: todo!(),
            allele_base_quality: todo!(),
            allele_map_quality: todo!(),
            strand_specific_base_quality: todo!(),
            strand_specific_mapping_quality: todo!(),
            position_in_read: todo!(),
            entropy: todo!(),
            num_aligned_bases: todo!(),
            num_indels: todo!(),
            in_cp_g: todo!(),
            de_novo_cp_g_candidate: todo!(),
        };
        let format_fields = Format {
            genotype: todo!(),
            genotype_likelihood: todo!(),
            genotype_confidence: todo!(),
            sample_read_depth: todo!(),
            methylated: todo!(),
            machine_learning_prediction: todo!(),
        };
        Ok(Record { main, filters, info, samples: smallvec_inline![format_fields] })
    }
}
