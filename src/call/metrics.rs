use crate::{
    call::pileup::{Pileup, SimpleReads},
    utils::{Base, ByStrand, Phred, RootMeanSquare, Strand},
    vcf::*,
};
use color_eyre::Result;
use rastair_vcf::standard_fields::*;
use smallvec::{SmallVec, smallvec_inline};
use std::collections::BTreeSet;
use tracing::warn;

mod entropy;

impl Pileup {
    pub fn fixed_fields(&self) -> rastair_vcf::VcfFixedFields {
        rastair_vcf::VcfFixedFields {
            chrom: self.chrom(),
            pos: self.pos,
            id: BTreeSet::default(),
            r#ref: self.reference_base.into(),
            alt: self.alts().iter().map(|b| (*b).into()).collect(),
            // TODO: Figure out how to handle this. When do we have the data for
            // this? Should we start with `None`?
            qual: self.qual(),
        }
    }

    pub fn metrics(&self) -> Result<Info> {
        Ok(Info {
            allele_read_depth: self.read_depth_per_allele(),
            allele_specific_strand_bias: self.allele_specific_strand_bias(),
            base_quality: BaseQuality(
                self.reads.iter().map(|b| b.qual).collect::<RootMeanSquare>(),
            ),
            mapping_quality: MappingQuality(
                self.reads.iter().map(|b| b.mapq).collect::<RootMeanSquare>(),
            ),
            read_depth: ReadDepth(self.reads.len()),
            mapping_quality0: MappingQuality0(self.reads.iter().filter(|b| b.mapq == 0).count()),
            // by construction, we arrived here because we have at least one base
            samples_with_data: SamplesWithData(1),
            sequence_context: self.sequence_context(),
            allele_frequency: self.allel_frequency(),
            allele_base_quality: self.allele_base_quality(),
            allele_map_quality: self.allele_map_quality(),
            strand_specific_base_quality: self.strand_specific_base_quality(),
            strand_specific_mapping_quality: self.strand_specific_mapping_quality(),
            position_in_read: self.position_in_read(),
            entropy: Entropy(smallvec_inline!(self.entropy())),
            num_aligned_bases: self.num_aligned_bases(),
            num_indels: self.num_indels(),
            in_cp_g: self.in_cpg(),
            de_novo_cp_g_candidate: self.de_novo_cpg(),
        })
    }

    fn num_indels(&self) -> NumIndels {
        NumIndels(
            self.by_allele()
                .iter()
                .map(|seen| *seen.iter().map(|b| b.indels).collect::<RootMeanSquare>())
                .collect(),
        )
    }

    fn num_aligned_bases(&self) -> NumAlignedBases {
        NumAlignedBases(
            self.by_allele()
                .iter()
                .map(|seen| *seen.iter().map(|b| b.matching_bases).collect::<RootMeanSquare>())
                .collect(),
        )
    }

    fn position_in_read(&self) -> PositionInRead {
        PositionInRead(
            self.by_allele()
                .iter()
                .map(|seen| {
                    *seen
                        .value
                        .iter()
                        .map(|b| f64::from(b.position.pos) / f64::from(b.position.read_length))
                        .collect::<RootMeanSquare>()
                })
                .collect(),
        )
    }

    fn allele_map_quality(&self) -> AlleleMapQuality {
        AlleleMapQuality(
            self.by_allele()
                .iter()
                .map(|seen| *seen.iter().map(|b| b.mapq).collect::<RootMeanSquare>())
                .collect(),
        )
    }

    fn allele_base_quality(&self) -> AlleleBaseQuality {
        AlleleBaseQuality(
            self.by_allele()
                .iter()
                .map(|seen| *seen.iter().map(|b| b.qual).collect::<RootMeanSquare>())
                .collect(),
        )
    }

    pub fn allel_frequency(&self) -> AlleleFrequency {
        AlleleFrequency(
            self.alts()
                .iter()
                .map(|alt| {
                    let count = self.reads.iter().filter(|b| b.base == *alt).count();
                    count as f64 / self.reads.len() as f64
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
    // TODO: return Phred<f32>
    pub fn qual(&self) -> Option<f32> {
        // TODO: Calculate this properly based on the pileup
        let probability_call_wrong = 0.001;
        // self.bases.iter().fold(1.0, |acc, b| acc * (f64::from(b.qual) / 10.0).powf(-1.0));
        match Phred::from_probability(probability_call_wrong) {
            Ok(phred) => Some(*phred as f32),
            Err(error) => {
                warn!(
                    error = format!("{error:#}"),
                    "Failed to calculate quality score for variant"
                );
                None
            }
        }
    }

    pub fn read_depth_per_allele(&self) -> AlleleReadDepth {
        fn count_bases(bases: &SimpleReads, base: Base) -> usize {
            bases.iter().filter(|b| b.base == base).count()
        }

        let mut depth = SmallVec::new();
        depth.push(count_bases(&self.reads, self.reference_base));
        for alt in self.alts() {
            depth.push(count_bases(&self.reads, alt));
        }
        AlleleReadDepth(depth)
    }

    pub fn allele_specific_strand_bias(&self) -> AlleleSpecificStrandBias {
        AlleleSpecificStrandBias(
            self.by_allele()
                .iter()
                .map(|seen| {
                    let ots = seen.iter().filter(|b| b.strand == Strand::OT).count();
                    let obs = seen.iter().filter(|b| b.strand == Strand::OB).count();
                    ByStrand {
                        base: seen.base,
                        ot: u32::try_from(ots).expect("count should fit in u32"),
                        ob: u32::try_from(obs).expect("count should fit in u32"),
                    }
                })
                .collect(),
        )
    }

    pub fn strand_specific_base_quality(&self) -> StrandSpecificBaseQuality {
        StrandSpecificBaseQuality(
            self.by_allele()
                .iter()
                .map(|seen| {
                    let ots = seen
                        .iter()
                        .filter(|b| b.strand == Strand::OT)
                        .map(|b| b.qual)
                        .collect::<RootMeanSquare>();
                    let obs = seen
                        .iter()
                        .filter(|b| b.strand == Strand::OB)
                        .map(|b| b.qual)
                        .collect::<RootMeanSquare>();
                    ByStrand { base: seen.base, ot: *ots, ob: *obs }
                })
                .collect(),
        )
    }

    pub fn strand_specific_mapping_quality(&self) -> StrandSpecificMappingQuality {
        StrandSpecificMappingQuality(
            self.by_allele()
                .iter()
                .map(|seen| {
                    let ots = seen
                        .iter()
                        .filter(|b| b.strand == Strand::OT)
                        .map(|b| b.mapq)
                        .collect::<RootMeanSquare>();
                    let obs = seen
                        .iter()
                        .filter(|b| b.strand == Strand::OB)
                        .map(|b| b.mapq)
                        .collect::<RootMeanSquare>();
                    ByStrand { base: seen.base, ot: *ots, ob: *obs }
                })
                .collect(),
        )
    }

    pub fn sequence_context(&self) -> SequenceContext {
        self.context.clone()
    }

    fn in_cpg(&self) -> InCpG {
        InCpG::from(self)
    }

    fn de_novo_cpg(&self) -> DeNovoCpGCandidate {
        DeNovoCpGCandidate::from(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::call::test_helpers::variant_pileup;
    use color_eyre::Result;
    use insta::assert_debug_snapshot;

    #[test]
    fn test_allele_specific_strand_bias_1() -> Result<()> {
        let pileup = variant_pileup("bacteriophage_lambda_CpG", 2636)?;
        assert_debug_snapshot!((
            pileup.chrom(),
            pileup.pos,
            pileup.reference_base,
            &pileup.reads,
            pileup.allele_specific_strand_bias()
        ), @r#"
        (
            "bacteriophage_lambda_CpG",
            2636,
            C,
            [
                C OB Q32 MQ60,
                C OB Q36 MQ60,
                T OT Q36 MQ60,
                C OB Q36 MQ60,
                T OT Q36 MQ60,
                T OT Q36 MQ60,
                T OT Q36 MQ60,
                T OT Q36 MQ60,
                C OB Q36 MQ60,
            ],
            AlleleSpecificStrandBias(
                [
                    ByStrand {
                        base: C,
                        ot: 0,
                        ob: 4,
                    },
                    ByStrand {
                        base: T,
                        ot: 5,
                        ob: 0,
                    },
                ],
            ),
        )
        "#);
        Ok(())
    }

    #[test]
    fn test_in_cpg() -> Result<()> {
        // a CpG site
        let pileup = variant_pileup("chr19", 6105084)?;
        assert_debug_snapshot!((pileup.chrom(), pileup.reference_base, pileup.pos, pileup.in_cpg()), @r#"
        (
            "chr19",
            C,
            6105084,
            CpG::C,
        )
        "#);

        // a C variant, followed by a C
        let pileup = variant_pileup("chr19", 6104589)?;
        assert_debug_snapshot!((pileup.chrom(), pileup.reference_base, pileup.pos, pileup.in_cpg()), @r#"
        (
            "chr19",
            C,
            6104589,
            NoCpg,
        )
        "#);

        // some random variant with base G
        let pileup = variant_pileup("chr19", 6105116)?;
        assert_debug_snapshot!((pileup.chrom(), pileup.reference_base, pileup.pos, pileup.in_cpg()), @r#"
        (
            "chr19",
            G,
            6105116,
            NoCpg,
        )
        "#);
        Ok(())
    }
}
