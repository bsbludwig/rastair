use std::collections::BTreeSet;

use super::{scores::VariantCandidatePileupMetrics, variants::VariantCandidatePileup};
use crate::{
    call::vcf::{self, ReadDepth},
    utils::RootMeanSquare,
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
                alt: self.0.bases.alts(),
                qual: 0., // TODO: Calculate quality score
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
}
