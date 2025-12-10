use super::{Methylated, Record};
use crate::{
    call::{RecordFilters, pileup::Pileup},
    metrics::{AlleleMetrics, Alt, Filters, PileupMetrics, PositionMetrics},
    utils::{IntoF64 as _, default},
    vcf::{
        AlleleBaseQuality, AlleleMapQuality, AlleleSpecificStrandBias, DeNovoCpGCandidate, Entropy,
        Filters as VcfFilter, Format, GenotypeConfidence, GenotypeLikelihood, Info,
        MachineLearningPrediction, NumAlignedBases, NumIndels, PositionInRead,
        StrandSpecificBaseQuality, StrandSpecificMappingQuality,
    },
};
use color_eyre::Result;
use color_eyre::eyre::ensure;
use rastair_types::{
    Phred, Probability,
    smallvec::{SmallVec, smallvec, smallvec_inline},
};
use rastair_vcf::{
    VcfFilter as _, VcfFixedFields,
    standard_fields::{
        AlleleFrequency, AlleleReadDepth, BaseQuality, Genotype, GenotypeAllele, MappingQuality,
        MappingQuality0, PASS, ReadDepth, SampleReadDepth, SamplesWithData,
    },
};

impl PileupMetrics {
    /// Convert the metrics to VCF records
    ///
    /// We write this
    /// 1. write one VCF record with all passing alts
    ///    1. if there is methylation evidence, we add `.` as an alt to indicate methylation status
    ///    2. if no alts pass, we write a ref-only record (ref->.). people want this for CpG sites.
    /// 2. write additional VCF records for each failing alt (may be filtered out later)
    ///
    // TODO: Handle methylation evidence when both a cpg and de-novo cpg candidate are present
    pub fn to_vcf_records(&self, ml_threshold: Option<Probability>) -> Result<VcfRecordSet> {
        // Validate that all alts have been called
        for alt in &self.alts {
            ensure!(
                !matches!(alt.call, crate::metrics::AltCall::Uncalled),
                "Alt {} at position {} is Uncalled - this should not happen",
                alt.base,
                self.pileup.pos
            );
        }

        // Categorize alts by their call type
        let mut real_variants = Vec::new();
        let mut methylation_evidence = Vec::new();
        let mut read_errors = Vec::new();

        for alt in &self.alts {
            match alt.call {
                crate::metrics::AltCall::RealVariant => real_variants.push(alt),
                crate::metrics::AltCall::MethylationEvidenceOnly { .. } => {
                    methylation_evidence.push(alt)
                }
                crate::metrics::AltCall::ReadError => read_errors.push(alt),
                crate::metrics::AltCall::Uncalled => unreachable!("checked above"),
            }
        }

        // Build main record
        let main = self.build_main_record(&real_variants, &methylation_evidence, ml_threshold)?;

        // Build rejected records (one per methylation evidence, one per read error)
        let mut rejected = Vec::new();
        for alt in methylation_evidence.iter().chain(read_errors.iter()) {
            rejected.push(self.build_rejected_record(alt, ml_threshold)?);
        }

        Ok(VcfRecordSet { main, rejected })
    }

    fn build_main_record(
        &self,
        real_variants: &[&Alt],
        methylation_evidence: &[&Alt],
        ml_threshold: Option<Probability>,
    ) -> Result<Record> {
        // Determine alt alleles for main record
        let has_methylation = !methylation_evidence.is_empty();

        // Build VCF fixed fields
        let main = VcfFixedFields {
            chrom: self.pileup.contig().clone(),
            pos: self.pileup.pos,
            id: default(),
            r#ref: self.pileup.reference_base.into(),
            alt: has_methylation
                .then_some(rastair_types::SmolStr::new_inline("."))
                .into_iter()
                .chain(real_variants.iter().map(|alt| alt.base.into()))
                .collect(),
            qual: Some({
                #[allow(clippy::cast_possible_truncation, reason = "const")]
                {
                    // FIXME: use real quality
                    *Phred::from(Probability::new_panicky(0.001)) as f32
                }
            }),
        };

        // Build info fields
        let info = self.build_info(real_variants);

        // Build format fields
        let format = self.build_format(real_variants, has_methylation);

        // Build filters
        let filters = self.build_filters(real_variants, ml_threshold);

        Ok(Record { main, filters, info, samples: smallvec_inline![format] })
    }

    fn build_rejected_record(
        &self,
        alt: &Alt,
        ml_threshold: Option<Probability>,
    ) -> Result<Record> {
        let main = VcfFixedFields {
            chrom: self.pileup.contig().clone(),
            pos: self.pileup.pos,
            id: default(),
            r#ref: self.pileup.reference_base.into(),
            alt: smallvec![alt.base.into()],
            qual: Some({
                #[allow(clippy::cast_possible_truncation, reason = "const")]
                {
                    *Phred::from(Probability::new_panicky(0.001)) as f32
                }
            }),
        };

        let info = self.build_info(&[alt]);
        let format = self.build_format(&[alt], false);
        let filters = self.build_filters(&[alt], ml_threshold);

        Ok(Record { main, filters, info, samples: smallvec_inline![format] })
    }

    fn build_info(&self, alts: &[&Alt]) -> Info {
        let ref_alts_metrics: SmallVec<&AlleleMetrics, 3> = {
            let mut xs = smallvec![&self.ref_metrics];
            for alt in alts {
                xs.push(&alt.metrics);
            }
            xs
        };
        let alts_metrics = ref_alts_metrics.get(1..).unwrap_or(&[]);

        Info {
            allele_read_depth: AlleleReadDepth(
                ref_alts_metrics.iter().map(|m| m.depth as usize).collect(),
            ),
            base_quality: BaseQuality(self.pos_metrics.baseq),
            read_depth: ReadDepth(self.pos_metrics.depth as usize),
            mapping_quality: MappingQuality(self.pos_metrics.mapq),
            mapping_quality0: MappingQuality0(self.pos_metrics.mapq0 as usize),
            samples_with_data: SamplesWithData(1),
            allele_specific_strand_bias: AlleleSpecificStrandBias(
                ref_alts_metrics.iter().map(|m| m.strand_count).collect(),
            ),
            sequence_context: self.pileup.context.clone(),
            allele_frequency: AlleleFrequency(
                alts_metrics.iter().map(|m| m.allele_frequency.f()).collect(),
            ),
            allele_base_quality: AlleleBaseQuality(
                ref_alts_metrics.iter().map(|m| m.baseq.f()).collect(),
            ),
            allele_map_quality: AlleleMapQuality(
                ref_alts_metrics.iter().map(|m| m.mapq.f()).collect(),
            ),
            strand_specific_base_quality: StrandSpecificBaseQuality(
                ref_alts_metrics.iter().map(|m| m.baseq_s).collect(),
            ),
            strand_specific_mapping_quality: StrandSpecificMappingQuality(
                ref_alts_metrics.iter().map(|m| m.mapq_s).collect(),
            ),
            position_in_read: PositionInRead(
                ref_alts_metrics.iter().map(|m| m.position_in_read.f()).collect(),
            ),
            entropy: Entropy(smallvec_inline![self.pos_metrics.region_entropy]),
            num_aligned_bases: NumAlignedBases(
                ref_alts_metrics.iter().map(|m| m.num_aligned_bases.f()).collect(),
            ),
            num_indels: NumIndels(ref_alts_metrics.iter().map(|m| m.num_indels.f()).collect()),
            in_cp_g: self.pos_metrics.cpg,
            methylation_evidence_strand_info: self.pos_metrics.extended.methylation_strand_info,
            de_novo_cp_g_candidate: {
                if let Some(alt) = alts_metrics.iter().find(|m| *m.denovo) {
                    DeNovoCpGCandidate::Candidate {
                        ref_base: self.pileup.reference_base,
                        alt_base: alt.base,
                    }
                } else if alts.is_empty() && *self.pos_metrics.denovo_adj {
                    DeNovoCpGCandidate::Adjecent { ref_base: self.pileup.reference_base }
                } else {
                    DeNovoCpGCandidate::NotCandidate
                }
            },
        }
    }

    fn build_format(&self, real_variants: &[&Alt], has_methylation: bool) -> Format {
        // Calculate genotype from real variants
        let (genotype, genotype_likelihood, genotype_confidence) =
            calculate_genotype(real_variants, &self.ref_metrics, self.pos_metrics.depth);

        let has_ml = real_variants.iter().any(|alt| alt.filters.ml.is_some());

        // Determine methylation status
        let methylated = if has_methylation || real_variants.iter().any(|alt| *alt.metrics.denovo) {
            self.pos_metrics.methylated.clone()
        } else {
            Methylated::Unknown
        };

        Format {
            genotype,
            genotype_likelihood,
            genotype_confidence,
            sample_read_depth: SampleReadDepth(self.pileup.reads.len()),
            methylated,
            machine_learning_prediction: MachineLearningPrediction(if has_ml {
                real_variants.iter().map(|alt| *alt.filters.ml.unwrap_or_default()).collect()
            } else {
                smallvec![]
            }),
        }
    }

    fn build_filters(&self, alts: &[&Alt], ml_threshold: Option<Probability>) -> super::Filters {
        let mut filters = super::Filters::default();

        let passes = if self.pos_filters.other_pos_in_denovo_passes {
            true
        } else {
            self.pos_filters.pass() && alts.iter().any(|a| a.filters.pass(ml_threshold))
        };

        if passes {
            filters.add(PASS.filter());
        } else {
            // Add position-level filters
            self.pos_filters.iter().for_each(|f| {
                filters.add(f.clone());
            });
            // Add alt-specific filters
            alts.iter().for_each(|alt| {
                alt.filters.filters.iter().for_each(|f| {
                    filters.add(f.clone());
                });
            });
        }

        filters
    }
}

/// Calculate genotype from real variant alts
///
/// Simple genotype caller: if we have real variants, we assume they're real.
/// This will be replaced with ML-based genotyping in the future.
fn calculate_genotype(
    real_variants: &[&Alt],
    _ref_metrics: &AlleleMetrics,
    _total_depth: u32,
) -> (Genotype, GenotypeLikelihood, GenotypeConfidence) {
    use GenotypeAllele::*;

    let gt = match real_variants.len() {
        0 => {
            // No real variants: homozygous reference (0/0)
            smallvec![Unphased(0), Unphased(0)]
        }
        1 => {
            // Single real variant: assume heterozygous (0/1)
            smallvec![Unphased(0), Unphased(1)]
        }
        _ => {
            // Multiple real variants: assume compound heterozygous (1/2)
            smallvec![Unphased(1), Unphased(2)]
        }
    };

    (
        Genotype(gt),
        GenotypeLikelihood(smallvec_inline![None]),
        GenotypeConfidence(smallvec_inline![None]),
    )
}

pub struct VcfRecordSet {
    main: Record,
    rejected: Vec<Record>,
}

impl VcfRecordSet {
    pub fn into_iter(self, filters: &RecordFilters) -> Box<dyn Iterator<Item = Record>> {
        match (filters.vcf_all, filters.cpgs_only) {
            (false, false) => {
                // default behavior: only passing records with alts
                if !self.main.main.alt.is_empty() {
                    Box::new(Some(self.main).into_iter())
                } else {
                    Box::new(std::iter::empty())
                }
            }
            (false, true) => {
                // CpGs with read evidence only
                if *self.main.info.read_depth > 0
                    && (*self.main.info.in_cp_g || *self.main.info.de_novo_cp_g_candidate)
                {
                    Box::new(Some(self.main).into_iter())
                } else {
                    Box::new(std::iter::empty())
                }
            }
            (true, false) => {
                // `--all` means all, do nothing
                Box::new(Some(self.main).into_iter().chain(self.rejected))
            }
            (true, true) => {
                if *self.main.info.in_cp_g {
                    Box::new(Some(self.main).into_iter().chain(self.rejected))
                } else {
                    Box::new(std::iter::empty())
                }
            }
        }
    }
}

#[deprecated]
mod old {
    use super::*;

    /// A subset of metrics needed to write a VCF record.
    /// Each instance represents one row in the VCF output.
    struct MetricsSubset<'m> {
        pub pileup: &'m Pileup,
        pub pos_metrics: &'m PositionMetrics,
        pub pos_filters: &'m Filters,
        pub ref_metrics: &'m AlleleMetrics,
        /// Alt alleles for this row. Empty for ref-only rows (ref->.).
        pub alts: SmallVec<&'m Alt, 2>,
        /// True if this is a ref-only row (ref->.) representing methylation evidence or no variants.
        pub is_ref_only_row: bool,
    }

    impl MetricsSubset<'_> {
        pub fn to_vcf_row(&self, ml_threshold: Option<Probability>) -> Record {
            let main = self.vcf_main();
            let info = self.info();
            let format_fields = self.format();
            let filters = self.filters(ml_threshold);

            Record { main, filters, info, samples: smallvec_inline![format_fields] }
        }

        fn vcf_main(&self) -> VcfFixedFields {
            VcfFixedFields {
                chrom: self.pileup.contig().clone(),
                pos: self.pileup.pos,
                id: default(),
                r#ref: self.pileup.reference_base.into(),
                alt: self.alts.iter().map(|alt| alt.base.into()).collect(),
                qual: Some(
                    #[allow(clippy::cast_possible_truncation, reason = "const")]
                    {
                        // FIXME: use real quality
                        *Phred::from(Probability::new_panicky(0.001)) as f32
                    },
                ),
            }
        }

        pub fn pass(&self, ml_threshold: Option<Probability>) -> bool {
            if self.pos_filters.other_pos_in_denovo_passes {
                return true;
            }

            // For alt rows, check position filters and that at least one alt passes
            self.pos_filters.pass() && self.alts.iter().any(|a| a.filters.pass(ml_threshold))
        }

        fn filters(&self, ml_threshold: Option<Probability>) -> VcfFilter {
            let mut filters = VcfFilter::default();

            if self.pass(ml_threshold) {
                filters.add(PASS.filter());
            } else {
                // Add position-level filters
                self.pos_filters.iter().for_each(|f| {
                    filters.add(f.clone());
                });
                // Add alt-specific filters
                self.alts.iter().for_each(|alt| {
                    alt.filters.filters.iter().for_each(|f| {
                        filters.add(f.clone());
                    });
                });
            }

            filters
        }

        fn info(&self) -> Info {
            let pileup = self.pileup;
            let pos_metrics = self.pos_metrics;
            let ref_alts_metrics: SmallVec<&AlleleMetrics, 3> = {
                let mut xs = smallvec![self.ref_metrics,];
                for alt in &self.alts {
                    xs.push(&alt.metrics);
                }
                xs
            };
            let alts_metrics = &ref_alts_metrics[1..];

            Info {
                allele_read_depth: AlleleReadDepth(
                    ref_alts_metrics.iter().map(|m| m.depth as usize).collect(),
                ),
                base_quality: BaseQuality(pos_metrics.baseq),
                read_depth: ReadDepth(pos_metrics.depth as usize),
                mapping_quality: MappingQuality(pos_metrics.mapq),
                mapping_quality0: MappingQuality0(pos_metrics.mapq0 as usize),
                samples_with_data: SamplesWithData(1),
                allele_specific_strand_bias: AlleleSpecificStrandBias(
                    ref_alts_metrics.iter().map(|m| m.strand_count).collect(),
                ),
                sequence_context: pileup.context.clone(),
                allele_frequency: AlleleFrequency(
                    alts_metrics.iter().map(|m| m.allele_frequency.f()).collect(),
                ),
                allele_base_quality: AlleleBaseQuality(
                    ref_alts_metrics.iter().map(|m| m.baseq.f()).collect(),
                ),
                allele_map_quality: AlleleMapQuality(
                    ref_alts_metrics.iter().map(|m| m.mapq.f()).collect(),
                ),
                strand_specific_base_quality: StrandSpecificBaseQuality(
                    ref_alts_metrics.iter().map(|m| m.baseq_s).collect(),
                ),
                strand_specific_mapping_quality: StrandSpecificMappingQuality(
                    ref_alts_metrics.iter().map(|m| m.mapq_s).collect(),
                ),
                position_in_read: PositionInRead(
                    ref_alts_metrics.iter().map(|m| m.position_in_read.f()).collect(),
                ),
                entropy: Entropy(smallvec_inline![pos_metrics.region_entropy]),
                num_aligned_bases: NumAlignedBases(
                    ref_alts_metrics.iter().map(|m| m.num_aligned_bases.f()).collect(),
                ),
                num_indels: NumIndels(ref_alts_metrics.iter().map(|m| m.num_indels.f()).collect()),
                in_cp_g: self.pos_metrics.cpg,
                methylation_evidence_strand_info: pos_metrics.extended.methylation_strand_info,
                de_novo_cp_g_candidate: {
                    if let Some(alt) = alts_metrics.iter().find(|m| *m.denovo) {
                        DeNovoCpGCandidate::Candidate {
                            ref_base: pileup.reference_base,
                            alt_base: alt.base,
                        }
                    } else if self.is_ref_only_row && *pos_metrics.denovo_adj {
                        DeNovoCpGCandidate::Adjecent { ref_base: pileup.reference_base }
                    } else {
                        DeNovoCpGCandidate::NotCandidate
                    }
                },
            }
        }

        fn format(&self) -> Format {
            let (genotype, genotype_likelihood, genotype_confidence) = if let Some(estimate) =
                self.pos_metrics.genotype
            {
                (
                    Genotype(<[GenotypeAllele; 2]>::from(estimate.genotype).into()),
                    GenotypeLikelihood(smallvec_inline![Some(Phred::from(estimate.likelihood))]),
                    GenotypeConfidence(smallvec_inline![Some(Phred::from(estimate.confidence))]),
                )
            } else {
                (
                    Genotype(smallvec![]),
                    GenotypeLikelihood(smallvec_inline![None]),
                    GenotypeConfidence(smallvec_inline![None]),
                )
            };

            let has_ml = self.alts.iter().any(|alt| alt.filters.ml.is_some());

            let methylated = if self.alts.is_empty() || self.is_ref_only_row {
                // For ref->. rows, indicate methylation status based on position metrics
                self.pos_metrics.methylated.clone()
            } else if self.alts.iter().any(|alt| *alt.metrics.denovo) {
                // For N->C or N->G de-novo CpG candidates, indicate methylation status based on alt metrics
                self.pos_metrics.methylated.clone()
            } else {
                Methylated::Unknown
            };

            Format {
                genotype,
                genotype_likelihood,
                genotype_confidence,
                sample_read_depth: SampleReadDepth(self.pileup.reads.len()),
                methylated,
                machine_learning_prediction: MachineLearningPrediction(if has_ml {
                    self.alts.iter().map(|alt| *alt.filters.ml.unwrap_or_default()).collect()
                } else {
                    smallvec![]
                }),
            }
        }
    }
}
