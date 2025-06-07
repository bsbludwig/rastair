use crate::{
    call::{
        variants::{SeenBases, VariantCandidatePileup},
        vcf::*,
    },
    utils::{Base, Counter, Phred, RootMeanSquare},
};
use color_eyre::{Result, eyre::ContextCompat as _};
use rastair2_vcf::standard_fields::*;
use smallvec::{SmallVec, smallvec, smallvec_inline};
use smol_str::{SmolStr, SmolStrBuilder};
use std::collections::BTreeSet;
use tracing::warn;

impl VariantCandidatePileup {
    pub fn fixed_fields(&self) -> rastair2_vcf::VcfFixedFields {
        rastair2_vcf::VcfFixedFields {
            chrom: self.chrom(),
            pos: self.pos,
            id: BTreeSet::default(),
            r#ref: self.reference_base.into(),
            alt: self.bases.alts(self.reference_base).iter().map(|b| (*b).into()).collect(),
            // TODO: Figure out how to handle this. When do we have the data for
            // this? Should we start with `None`?
            qual: self.qual(),
        }
    }

    pub fn metrics(&self) -> Result<Info> {
        Ok(Info {
            read_depth_per_allel: self.read_depth_per_allele(),
            strand_bias: self.strand_bias(),
            base_quality: BaseQuality(smallvec![*RootMeanSquare::new(
                &self.bases.iter().map(|b| b.qual).collect::<SmallVec<u8, 30>>(),
            )]),
            mapping_quality: MappingQuality(smallvec![*RootMeanSquare::new(
                &self.bases.iter().map(|b| b.mapq).collect::<SmallVec<u8, 30>>(),
            )]),
            read_depth: ReadDepth(smallvec![self.bases.len()]),
            mapping_quality0: MappingQuality0(smallvec![
                self.bases.iter().filter(|b| b.mapq == 0).count()
            ]),
            // by construction, we arrived here because we have at least one base
            samples_with_data: SamplesWithData(smallvec_inline![1]),
            sequence_context: SequenceContext(smallvec![self.sequence_context()]),
            allel_frequency: self.allel_frequency(),
            allel_base_quality: self.allel_base_quality(),
            allel_map_quality: self.allel_map_quality(),
            position_in_read: self.position_in_read(),
            entropy: self.entropy(),
            num_aligned_bases: self.num_aligned_bases(),
            num_indels: self.num_indels(),
        })
    }

    fn num_indels(&self) -> NumIndels {
        NumIndels(
            self.bases
                .alleles()
                .iter()
                .map(|alt| {
                    *RootMeanSquare::new(
                        &self
                            .bases
                            .iter()
                            .filter(|b| b.base == *alt)
                            .map(|b| b.indels)
                            .collect::<SmallVec<u32, 20>>(),
                    )
                })
                .collect(),
        )
    }

    fn num_aligned_bases(&self) -> NumAlignedBases {
        NumAlignedBases(
            self.bases
                .alleles()
                .iter()
                .map(|alt| {
                    *RootMeanSquare::new(
                        &self
                            .bases
                            .iter()
                            .filter(|b| b.base == *alt)
                            .map(|b| b.matching_bases)
                            .collect::<SmallVec<u32, 20>>(),
                    )
                })
                .collect(),
        )
    }

    ///  Calculate Shannon entropy for 100bp context around variant position
    fn entropy(&self) -> Entropy {
        let pos = usize::try_from(self.pos).expect("pos fits usize");
        let idx = usize::try_from(self.pos)
            .expect("position fits in usize")
            .checked_sub(usize::try_from(self.segment.range.start).expect("index fits in usize"))
            .wrap_err_with(|| {
                format!(
                    "pile position {} is not in segment range {}..{}",
                    pos, self.segment.range.start, self.segment.range.end
                )
            })
            .expect("valid index");

        let start = idx.saturating_sub(50);
        let end = idx.saturating_add(50);
        let seq_context = self.segment.get(start, end).expect("sequence context indices are valid");
        let counts: Counter = seq_context.iter().filter_map(|&b| Base::try_from(b).ok()).collect();
        let total = seq_context.len() as f64;
        let entropy = counts
            .entries()
            .iter()
            .map(|(_base, count)| {
                let p = *count as f64 / total;
                -p * p.log2()
            })
            .sum::<f64>();

        Entropy(smallvec![entropy])
    }

    fn position_in_read(&self) -> PositionInRead {
        PositionInRead(
            self.bases
                .alleles()
                .iter()
                .map(|alt| {
                    *RootMeanSquare::new(
                        self.bases
                            .iter()
                            .filter(|b| b.base == *alt)
                            .map(|b| f64::from(b.position.pos) / f64::from(b.position.read_length))
                            .collect::<SmallVec<f64, 20>>()
                            .as_slice(),
                    )
                })
                .collect(),
        )
    }

    fn allel_map_quality(&self) -> AllelMapQuality {
        AllelMapQuality(
            self.bases
                .alleles()
                .iter()
                .map(|alt| {
                    *RootMeanSquare::new(
                        &self
                            .bases
                            .iter()
                            .filter(|b| b.base == *alt)
                            .map(|b| b.mapq)
                            .collect::<SmallVec<u8, 20>>(),
                    )
                })
                .collect(),
        )
    }

    fn allel_base_quality(&self) -> AllelBaseQuality {
        AllelBaseQuality(
            self.bases
                .alleles()
                .iter()
                .map(|alt| {
                    *RootMeanSquare::new(
                        &self
                            .bases
                            .iter()
                            .filter(|b| b.base == *alt)
                            .map(|b| b.qual)
                            .collect::<SmallVec<u8, 20>>(),
                    )
                })
                .collect(),
        )
    }

    fn allel_frequency(&self) -> AllelFrequency {
        AllelFrequency(
            self.bases
                .alts(self.reference_base)
                .iter()
                .map(|alt| {
                    let count = self.bases.iter().filter(|b| b.base == *alt).count();
                    count as f64 / self.bases.len() as f64
                })
                .collect(),
        )
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
    fn qual(&self) -> Option<f32> {
        // TODO: Calculate this properly based on the pileup
        let probability_call_wrong = 0.001;
        // self.bases.iter().fold(1.0, |acc, b| acc * (f64::from(b.qual) / 10.0).powf(-1.0));
        match Phred::new(probability_call_wrong) {
            Ok(phred) => Some(*phred as f32),
            Err(error) => {
                warn!(%error, "Failed to calculate quality score for variant");
                None
            }
        }
    }

    fn read_depth_per_allele(&self) -> ReadDepthPerAllel {
        fn count_bases(bases: &SeenBases, base: Base) -> usize {
            bases.iter().filter(|b| b.base == base).count()
        }

        let mut depth = SmallVec::new();
        depth.push(count_bases(&self.bases, self.reference_base));
        for alt in self.bases.alts(self.reference_base) {
            depth.push(count_bases(&self.bases, alt));
        }
        ReadDepthPerAllel(depth)
    }

    fn strand_bias(&self) -> StrandBias {
        let reference_bases = self.bases.iter().filter(|b| b.base == self.reference_base);
        let alt_bases = self.bases.iter().filter(|b| b.base != self.reference_base);

        StrandBias {
            reads_ref_fwd: reference_bases.clone().filter(|b| !b.reverse).count(),
            reads_ref_rev: reference_bases.clone().filter(|b| b.reverse).count(),
            reads_alt_fwd: alt_bases.clone().filter(|b| !b.reverse).count(),
            reads_alt_rev: alt_bases.clone().filter(|b| b.reverse).count(),
        }
    }

    fn sequence_context(&self) -> SmolStr {
        let mut builder = SmolStrBuilder::new();
        for base in self.sequence_before.iter() {
            builder.push_str((*base).into());
        }
        builder.push_str(self.reference_base.into());
        for base in self.sequence_after.iter() {
            builder.push_str((*base).into());
        }
        builder.finish()
    }
}
