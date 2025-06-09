use crate::call::{variants::VariantCandidatePileup, vcf::Format};
use color_eyre::eyre::Result;
use rastair2_vcf::standard_fields::*;

impl VariantCandidatePileup {
    pub fn calling_metrics(&self) -> Result<Format> {
        Ok(Format { sample_read_depth: SampleReadDepth(self.bases.len()) })
    }
}
