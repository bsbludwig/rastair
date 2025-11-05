use super::Record;
use crate::{
    metrics2::PileupMetrics,
    utils::conversion::IntoF64 as _,
    vcf::{
        AlleleBaseQuality, AlleleMapQuality, Entropy, Filters, Format, GenotypeConfidence,
        GenotypeLikelihood, InCpG, Info, MachineLearningPrediction, NumAlignedBases, NumIndels,
        PositionInRead,
    },
};
use color_eyre::Result;
use rastair_types::Phred;
use rastair_vcf::{
    VcfFixedFields,
    standard_fields::{
        AlleleReadDepth, BaseQuality, Genotype, GenotypeAllele, MappingQuality, MappingQuality0,
        PASS, ReadDepth, SampleReadDepth, SamplesWithData,
    },
};
use smallvec::{SmallVec, smallvec, smallvec_inline};

pub struct VcfOutputFilter {
    pub reject_low_quality_variants: bool,
}

impl PileupMetrics {
    pub fn to_vcf_record(&self) -> Result<Record> {
        let main = VcfFixedFields {
            chrom: self.contig(),
            pos: self.pos(),
            id: Default::default(),
            r#ref: self.ref_base().into(),
            alt: self.alts().iter().map(|alt| (*alt).into()).collect::<SmallVec<_, 2>>(),
            qual: self.pileup.qual(),
        };
        let info = Info {
            allele_read_depth: AlleleReadDepth(
                self.ref_alts_metrics().map(|m| m.depth as usize).collect(),
            ),
            base_quality: BaseQuality(self.pos_metrics.baseq),
            read_depth: ReadDepth(self.pos_metrics.read_depth),
            mapping_quality: MappingQuality(self.pos_metrics.mapq),
            mapping_quality0: MappingQuality0(self.pos_metrics.mapq0),
            samples_with_data: SamplesWithData(1),
            allele_specific_strand_bias: self.pileup.allele_specific_strand_bias(),
            sequence_context: self.pos_metrics.sequence_context.clone(),
            allele_frequency: self.pileup.allel_frequency(),
            allele_base_quality: AlleleBaseQuality(
                self.ref_alts_metrics().map(|m| m.baseq.f()).collect(),
            ),
            allele_map_quality: AlleleMapQuality(
                self.ref_alts_metrics().map(|m| m.mapq.f()).collect(),
            ),
            strand_specific_base_quality: self.pileup.strand_specific_base_quality(),
            strand_specific_mapping_quality: self.pileup.strand_specific_mapping_quality(),
            position_in_read: PositionInRead(
                self.ref_alts_metrics().map(|m| m.position_in_read.f()).collect(),
            ),
            entropy: Entropy(smallvec_inline![self.pos_metrics.region_entropy]),
            num_aligned_bases: NumAlignedBases(
                self.ref_alts_metrics().map(|m| m.num_aligned_bases.f()).collect(),
            ),
            num_indels: NumIndels(self.ref_alts_metrics().map(|m| m.num_indels.f()).collect()),
            in_cp_g: InCpG::from(&self.pileup),
            de_novo_cp_g_candidate: self.pos_metrics.de_novo_cpg_candidate.clone(),
        };

        let (genotype, genotype_likelihood, genotype_confidence) =
            if let Some(estimate) = self.pos_metrics.genotype {
                (
                    Genotype(<[GenotypeAllele; 2]>::from(estimate.genotype).into()),
                    GenotypeLikelihood(smallvec_inline![
                        Phred::from_probability(1.0 - estimate.likelihood).ok()
                    ]),
                    GenotypeConfidence(smallvec_inline![
                        Phred::from_probability(1.0 - estimate.confidence).ok()
                    ]),
                )
            } else {
                (
                    Genotype(smallvec![]),
                    GenotypeLikelihood(smallvec_inline![None]),
                    GenotypeConfidence(smallvec_inline![None]),
                )
            };

        let format_fields = Format {
            genotype,
            genotype_likelihood,
            genotype_confidence,
            sample_read_depth: SampleReadDepth(self.pileup.reads.len()),
            methylated: self.pos_metrics.methylated.clone(),
            machine_learning_prediction: MachineLearningPrediction(
                self.alts.iter().map(|alt| *alt.filters.ml.unwrap_or_default()).collect(),
            ),
        };

        // FIXME: Add real filters based on metrics
        let mut filters = Filters::default();
        self.pos_filters.iter().for_each(|f| {
            filters.add(f.clone());
        });
        if filters.pass() {
            filters.add_all(PASS);
        }

        Ok(Record { main, filters, info, samples: smallvec_inline![format_fields] })
    }
}
