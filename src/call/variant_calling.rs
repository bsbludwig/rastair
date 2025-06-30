use crate::{call::variants::VariantCandidatePileup, vcf::*};
use color_eyre::eyre::Result;
use rastair2_vcf::standard_fields::*;

mod params;
pub use params::VariantCallingParams;
mod error_model;
pub use error_model::ErrorModel;
use smallvec::smallvec;
use tracing::instrument;
mod genotype;

impl VariantCandidatePileup {
    #[instrument(level="trace", skip_all, fields(chr = %self.chrom(), pos = self.pos))]
    pub fn calling_metrics(&self, params: &VariantCallingParams) -> Result<Format> {
        let (genotype, genotype_likelihood, genotype_confidence) =
            if let Some(estimate) = self.estimate_genotype(params.error_model) {
                (
                    Genotype(<[GenotypeAllele; 2]>::from(estimate.genotype).into()),
                    GenotypeLikelihood(smallvec![Some(estimate.likelihood)]),
                    GenotypeConfidence(smallvec![Some(estimate.confidence)]),
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
            sample_read_depth: SampleReadDepth(self.bases.len()),
            methylated: Methylated(None),
            de_novo_cpg: DeNovoCpg(None),
        })
    }
}
