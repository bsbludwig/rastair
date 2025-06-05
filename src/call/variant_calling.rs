use crate::call::{variants::VariantCandidatePileup, vcf::Format};
use color_eyre::eyre::Result;
use rastair2_vcf::standard_fields::*;
use smallvec::smallvec;

impl VariantCandidatePileup {
    pub fn calling_metrics(&self) -> Result<Format> {
        Ok(Format { sample_read_depth: SampleReadDepth(smallvec![self.bases.len()]) })
    }
}
