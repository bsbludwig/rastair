use crate::{call::pileup::Pileup, vcf::*};
use color_eyre::eyre::Result;
use rastair_types::Phred;
use rastair_vcf::standard_fields::*;

mod params;
pub use params::VariantCallingParams;
mod error_model;
pub use error_model::ErrorModel;
use smallvec::smallvec;
use tracing::instrument;
mod genotype;
pub use genotype::{EstimatedGenotype, GenotypeTag};
mod read_flags;
pub use read_flags::ReadFlags;
mod read_masking;
pub use read_masking::ReadMaskParams;
mod quality_filters;
pub use quality_filters::QualityFilterParams;

impl Pileup {
    #[instrument(level="trace", skip_all, fields(chr = %self.chrom(), pos = self.pos))]
    pub fn calling_metrics(&self, error_model: ErrorModel) -> Result<Format> {
        let (genotype, genotype_likelihood, genotype_confidence) =
            if let Some(estimate) = self.estimate_genotype(error_model) {
                (
                    Genotype(<[GenotypeAllele; 2]>::from(estimate.genotype).into()),
                    GenotypeLikelihood(smallvec![
                        Phred::from_probability(1.0 - estimate.likelihood).ok()
                    ]),
                    GenotypeConfidence(smallvec![
                        Phred::from_probability(1.0 - estimate.confidence).ok()
                    ]),
                )
            } else {
                (
                    Genotype(smallvec![]),
                    GenotypeLikelihood(smallvec![None]),
                    GenotypeConfidence(smallvec![None]),
                )
            };

        Ok(Format {
            genotype,
            genotype_likelihood,
            genotype_confidence,
            sample_read_depth: SampleReadDepth(self.reads.len()),
            methylated: Methylated::Unknown,
            machine_learning_prediction: MachineLearningPrediction(smallvec![]),
        })
    }
}
