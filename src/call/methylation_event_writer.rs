use super::{scores::VariantCandidatePileupMetrics, variants::VariantCandidatePileup};
use crate::{
    call::{
        variants::SeenBases,
        vcf::{self, SequenceContext},
    },
    utils::{Base, Phred, RootMeanSquare},
};
use color_eyre::eyre::Result;
use rastair2_vcf::{Vcf, VcfFixedFields, standard_fields::*};
use smallvec::{SmallVec, smallvec, smallvec_inline};
use smol_str::{SmolStr, SmolStrBuilder};
use std::collections::BTreeSet;
use tracing::instrument;

pub struct MethylationEventWriter<'p, 'm>(
    pub &'p VariantCandidatePileup,
    pub &'m VariantCandidatePileupMetrics,
);

impl MethylationEventWriter<'_, '_> {
    #[instrument(level = "trace", skip_all)]
    pub fn write(&self, w: &mut Vcf<vcf::Record>) -> Result<()> {
        let filters = vcf::Filters::new().add(PASS);

        let record = vcf::Record {
            fixed_fields: VcfFixedFields {
                chrom: self.0.chrom.clone(),
                pos: self.0.pos,
                id: BTreeSet::default(),
                r#ref: self.0.reference_base.into(),
                alt: self
                    .0
                    .bases
                    .alts(self.0.reference_base)
                    .into_iter()
                    .map(|b| b.into())
                    .collect(),
                qual: self.qual(),
            },
            // TODO: Add filters
            filters,
            info: vcf::Info {
                ReadDepthPerAllel: ReadDepthPerAllel(self.read_depth_per_allele()),
                StrandBias: self.strand_bias(),
                BaseQuality: BaseQuality(smallvec![*RootMeanSquare::new(
                    &self.0.bases.iter().map(|b| b.qual).collect::<SmallVec<u8, 20>>(),
                )]),
                MappingQuality: MappingQuality(smallvec![*RootMeanSquare::new(
                    &self.0.bases.iter().map(|b| b.mapq).collect::<SmallVec<u8, 20>>(),
                )]),
                ReadDepth: ReadDepth(smallvec![self.0.bases.len()]),
                MappingQuality0: MappingQuality0(smallvec![
                    self.0.bases.iter().filter(|b| b.mapq == 0).count()
                ]),
                // by construction, we arrived here because we have at least one base
                SamplesWithData: SamplesWithData(smallvec_inline![1]),
                SequenceContext: SequenceContext(smallvec![self.sequence_context()]),
                AllelFrequency: AllelFrequency(
                    self.0
                        .bases
                        .alts(self.0.reference_base)
                        .iter()
                        .map(|alt| {
                            let count = self.0.bases.iter().filter(|b| b.base == *alt).count();
                            count as f64 / self.0.bases.len() as f64
                        })
                        .collect(),
                ),
            },
            samples: smallvec::smallvec![vcf::Format {
                SampleReadDepth: SampleReadDepth(smallvec![self.0.bases.len()])
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
        // TODO: Calculate this properly based on the pileup
        let probability_call_wrong = 0.001;
        // self.0.bases.iter().fold(1.0, |acc, b| acc * (f64::from(b.qual) / 10.0).powf(-1.0));
        *Phred::new(probability_call_wrong).expect("probability must be finite") as f32
    }

    fn read_depth_per_allele(&self) -> SmallVec<usize, 3> {
        fn count_bases(bases: &SeenBases, base: Base) -> usize {
            bases.iter().filter(|b| b.base == base).count()
        }

        let mut depth = SmallVec::new();
        depth.push(count_bases(&self.0.bases, self.0.reference_base));
        for alt in self.0.bases.alts(self.0.reference_base) {
            depth.push(count_bases(&self.0.bases, alt));
        }
        depth
    }

    fn strand_bias(&self) -> StrandBias {
        StrandBias {
            reads_ref_fwd: self.1.strand_bias.reference_ot,
            reads_ref_rev: self.1.strand_bias.reference_ob,
            reads_alt_fwd: self.1.strand_bias.alt_ot,
            reads_alt_rev: self.1.strand_bias.alt_ob,
        }
    }

    fn sequence_context(&self) -> SmolStr {
        let mut builder = SmolStrBuilder::new();
        for base in self.0.sequence_before.iter() {
            builder.push_str((*base).into());
        }
        builder.push_str(self.0.reference_base.into());
        for base in self.0.sequence_after.iter() {
            builder.push_str((*base).into());
        }
        builder.finish()
    }
}

#[test]
fn record_size() {
    dbg!(size_of::<VcfFixedFields>());
    dbg!(size_of::<vcf::Filters>());
    dbg!(size_of::<vcf::Info>());
    dbg!(size_of::<vcf::Format>());
    dbg!(size_of::<vcf::Record>());
    assert!(size_of::<vcf::Record>() < 1024);
}
