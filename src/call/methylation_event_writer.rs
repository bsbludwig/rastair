use std::collections::BTreeSet;

use super::{scores::VariantCandidatePileupMetrics, variants::VariantCandidatePileup};
use crate::{
    call::vcf::{self, ReadDepth},
    utils::{Phred, RootMeanSquare},
};
use color_eyre::eyre::Result;
use rastair2_vcf::Vcf;
use smallvec::SmallVec;
use tracing::instrument;

/// TODO: Make this a proper VCF writer
/// cf. <https://samtools.github.io/hts-specs/VCFv4.5.pdf>
pub struct MethylationEventWriter<'p, 'm>(
    pub &'p VariantCandidatePileup,
    pub &'m VariantCandidatePileupMetrics,
);

impl MethylationEventWriter<'_, '_> {
    #[instrument(level = "trace", skip_all)]
    pub fn write(&self, w: &mut Vcf<vcf::Record>) -> Result<()> {
        let record = vcf::Record {
            fixed_fields: rastair2_vcf::VcfFixedFields {
                chrom: self.0.chrom.clone(),
                pos: self.0.pos,
                id: BTreeSet::default(),
                r#ref: self.0.reference_base.into(),
                alt: self.0.bases.alts(self.0.reference_base),
                qual: self.qual(),
            },
            // TODO: Add filters
            filters: vcf::Filters::new(),
            info: vcf::Info {
                BaseQuality: vcf::BaseQuality(vec![*RootMeanSquare::new(
                    &self.0.bases.iter().map(|b| b.qual).collect::<SmallVec<u8, 20>>(),
                )]),
                MappingQuality: vcf::MappingQuality(vec![*RootMeanSquare::new(
                    &self.0.bases.iter().map(|b| b.mapq).collect::<SmallVec<u8, 20>>(),
                )]),
            },
            samples: smallvec::smallvec![vcf::Format {
                ReadDepth: ReadDepth(vec![self.0.bases.len()])
            }],
        };
        w.add(record)?;

        Ok(())
    }

    /// Phred-scaled quality score for the assertion made in ALT
    ///
    /// > If ALT is `.` (no variant) then this is $-10\log_{10}$ prob(call in ALT is wrong).
    /// > If ALT is not `.` this is $-10\log_{10}$ prob(no variant).
    /// > If unknown, the MISSING value must be specified. (Float)
    /// >
    /// > [spec](https://github.com/samtools/hts-specs/blob/0d7f8774658f7cee0a4540b0682174e460726432/VCFv4.5.tex#L420C3-L422C59)
    ///
    /// NOTE: we are in the case where we have a variant
    #[allow(clippy::cast_possible_truncation)]
    fn qual(&self) -> f32 {
        let probability_call_wrong = 0.001; // TODO: Calculate this properly based on the pileup
        // self.0.bases.iter().fold(1.0, |acc, b| acc * (f64::from(b.qual) / 10.0).powf(-1.0));
        *Phred::new(probability_call_wrong).expect("probability must be finite") as f32
    }
}
