use super::{MethylationAltDepth, MethylationDepth, Record};
use crate::{
    call::{RecordFilters, pileup::indels::IndelAllele, variant_calling::ErrorModel},
    metrics::{AlleleMetrics, Alt, AltCall, PileupMetrics},
    utils::{IntoF64 as _, default},
    vcf::{
        AlleleBaseQuality, AlleleMapQuality, AlleleSpecificStrandBias, CpgBeta, CpgOrigin,
        DeNovoCpGCandidate, Entropy, Format, GenotypeConfidence, GenotypeLikelihood, Info,
        MachineLearningPrediction, Methylated, NumAlignedBases, NumIndels, PositionInRead,
        StrandSpecificBaseQuality, StrandSpecificMappingQuality, low_ml_score,
    },
};
use color_eyre::eyre::ensure;
use color_eyre::{Result, eyre::Context as _};
use rastair_vcf::{
    VcfFilter as _, VcfFixedFields,
    standard_fields::{
        AlleleFrequency, AlleleReadDepth, BaseQuality, Genotype, GenotypeAllele, MappingQuality,
        MappingQuality0, PASS, ReadDepth, SampleReadDepth, SamplesWithData,
    },
};
use seqair_types::{
    Phred, Probability,
    smallvec::{SmallVec, smallvec, smallvec_inline},
};
use std::num::NonZeroU8;

impl PileupMetrics {
    /// Convert the metrics to VCF records
    ///
    /// We write:
    /// 1. One main VCF record with real variants only
    ///    - If real variants exist: use only real variants
    ///    - If no real variants: use '.' for reference-only CpG tracking
    /// 2. Additional rejected records for methylation evidence and read errors
    ///
    /// Methylation information is preserved in the `M5mC` format field.
    ///
    // TODO: Handle methylation evidence when both a cpg and de-novo cpg candidate are present
    pub fn to_vcf_records(
        &self,
        ml_threshold: Option<Probability>,
        error_model: &ErrorModel,
    ) -> Result<VcfRecordSet<'_>> {
        // Validate that all alts have been called
        for alt in &self.alts {
            ensure!(
                !matches!(alt.call, AltCall::Uncalled),
                "Alt {} at position {} is Uncalled - this should not happen",
                alt.base,
                self.pileup.pos
            );
        }

        // Categorize alts by their call type
        let mut real_variants: SmallVec<_, 2> = SmallVec::new();
        let mut methylation_evidence: SmallVec<_, 2> = SmallVec::new();
        let mut read_errors: SmallVec<_, 2> = SmallVec::new();

        for alt in &self.alts {
            match alt.call {
                AltCall::RealVariant => real_variants.push(alt),
                AltCall::MethylationEvidenceOnly { .. } => methylation_evidence.push(alt),
                AltCall::ReadError => read_errors.push(alt),
                AltCall::Uncalled => unreachable!("checked above"),
            }
        }

        // Build main record
        let main = self.build_main_record(&real_variants, ml_threshold, error_model)?;

        // Build rejected records
        let mut rejected = SmallVec::new();

        // Create rejected records for methylation evidence and read errors
        for alt in methylation_evidence.iter().chain(read_errors.iter()) {
            rejected.push(self.build_rejected_record(alt, ml_threshold)?);
        }

        let indel_records = build_indel_records(self, ml_threshold);

        Ok(VcfRecordSet { pileup: self, main, rejected, indel_records })
    }

    fn build_main_record(
        &self,
        real_variants: &[&Alt],
        ml_threshold: Option<Probability>,
        error_model: &ErrorModel,
    ) -> Result<Record> {
        let depth = self.pos_metrics.depth.max(1).f();

        // Build alt alleles for main record:
        // - If we have real variants: use ONLY real variants (no '.' mixing)
        // - If no real variants: use '.' for reference-only records (needed for valid VCF)
        let alt_alleles: seqair_types::SmallVec<seqair_types::SmolStr, 2> =
            if !real_variants.is_empty() {
                // Case 1: Real variants exist - use only real variants (fixes the C .,G issue)
                real_variants.iter().map(|alt| alt.base.into()).collect()
            } else {
                // Case 2: No real variants - use '.' for valid VCF reference records
                seqair_types::SmallVec::from([".".into()])
            };

        // For info and format fields, use real_variants when they exist
        let main_alts = if !real_variants.is_empty() { real_variants } else { &[] };

        // Build VCF fixed fields
        let main = VcfFixedFields {
            chrom: self.pileup.contig().clone(),
            pos: self.pileup.pos,
            id: default(),
            r#ref: self.pileup.reference_base.into(),
            alt: alt_alleles,
            qual: {
                let ml_qual = if real_variants.is_empty() {
                    // No real variants, VCF spec says: QUAL = -10log10(P(variant))
                    // As there is _no_ evidence for a variant, the Phred should be MAX
                    Some(Phred::from_phred(99_u8).as_int())
                } else {
                    // There are variants, VCF spec says: QUAL = -10log10(P(no variant))
                    // Use the *inverted* maximum ML score from all real variants
                    real_variants
                        .iter()
                        .filter_map(|alt| alt.filters.ml)
                        .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                        .map(|p| Phred::from(p.inverted()).as_int())
                };
                // Fallback: if no ML scores available, use sequencing error rate
                if let Some(ml_qual) = ml_qual {
                    Some(ml_qual)
                } else {
                    Some(
                        Phred::from(
                            Probability::new(*error_model.error_rate() / depth)
                                .wrap_err("Failed to compute QUAL from error rate")?,
                        )
                        .as_int(),
                    )
                }
            },
        };

        // Build info fields
        let info = self.build_info(main_alts);

        // Build format fields
        let format = self.build_format(main_alts, ml_threshold);

        // Build filters
        let filters = {
            let mut filters = super::Filters::default();
            filters.add(PASS.filter());
            filters
        };

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
            qual: alt.filters.ml.map(|ml| Phred::from(ml.inverted()).as_int()),
        };

        let info = self.build_info(&[alt]);
        let format = self.build_format(&[alt], ml_threshold);

        let filters = {
            let mut filters = super::Filters::default();

            if alt.filters.ml < ml_threshold {
                filters.add(low_ml_score.filter());
            }
            // Add position-level filters
            self.pos_filters.iter().for_each(|f| {
                filters.add(f.clone());
            });
            // Add alt-specific filters
            alt.filters.filters.iter().for_each(|f| {
                filters.add(f.clone());
            });

            filters
        };

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
                let t = &self.tags;
                if t.denovo_cpg || t.denovo_cpg_partner {
                    if let Some(alt) = alts_metrics.iter().find(|m| *m.denovo) {
                        DeNovoCpGCandidate::Candidate {
                            ref_base: self.pileup.reference_base,
                            alt_base: alt.base,
                        }
                    } else {
                        DeNovoCpGCandidate::Adjecent { ref_base: self.pileup.reference_base }
                    }
                } else {
                    DeNovoCpGCandidate::NotCandidate
                }
            },
        }
    }

    fn build_format(&self, main_alts: &[&Alt], ml_threshold: Option<Probability>) -> Format {
        // No more "." alt, so no index offset needed
        let vcf_index_offset = 0;

        // Calculate genotype from real variants using estimate_genotype
        let (genotype, genotype_likelihood, genotype_confidence) = if let Some(estimated) =
            self.estimate_genotype(ml_threshold, ErrorModel::default())
        {
            // Build a mapping from self.alts index to main_alts index (position in VCF)
            let mut self_alts_to_vcf_index: SmallVec<_, 2> = smallvec![None; self.alts.len()];
            for (vcf_idx, main_alt) in main_alts.iter().enumerate() {
                // Find this alt in self.alts
                for (self_idx, self_alt) in self.alts.iter().enumerate() {
                    if std::ptr::eq(*main_alt, self_alt) {
                        self_alts_to_vcf_index[self_idx] = Some(vcf_idx + 1 + vcf_index_offset);
                        break;
                    }
                }
            }

            // Remap the genotype indices from self.alts positions to VCF positions
            use crate::call::variant_calling::GenotypeTag;
            #[expect(clippy::cast_possible_truncation, reason = "small indices only")]
            let remapped_genotype = match estimated.genotype {
                GenotypeTag::HomRef => GenotypeTag::HomRef,
                GenotypeTag::RefHet(n) => {
                    let self_idx = n.get() as usize - 1;
                    if let Some(vcf_idx) = self_alts_to_vcf_index.get(self_idx).and_then(|&x| x) {
                        GenotypeTag::ref_het(NonZeroU8::new(vcf_idx as u8).unwrap_or(n))
                    } else {
                        // Alt not in real_variants, default to 0/0
                        GenotypeTag::HomRef
                    }
                }
                GenotypeTag::HomAlt(n) => {
                    let self_idx = n.get() as usize - 1;
                    if let Some(vcf_idx) = self_alts_to_vcf_index.get(self_idx).and_then(|&x| x) {
                        GenotypeTag::hom_alt(NonZeroU8::new(vcf_idx as u8).unwrap_or(n))
                    } else {
                        // Alt not in real_variants, default to 0/0
                        GenotypeTag::HomRef
                    }
                }
                GenotypeTag::AltHet(m, n) => {
                    let self_idx_m = m.get() as usize - 1;
                    let self_idx_n = n.get() as usize - 1;
                    let vcf_idx_m = self_alts_to_vcf_index.get(self_idx_m).and_then(|&x| x);
                    let vcf_idx_n = self_alts_to_vcf_index.get(self_idx_n).and_then(|&x| x);

                    match (vcf_idx_m, vcf_idx_n) {
                        (Some(vm), Some(vn)) => GenotypeTag::alt_het(
                            NonZeroU8::new(vm as u8).unwrap_or(m),
                            NonZeroU8::new(vn as u8).unwrap_or(n),
                        ),
                        (Some(vm), None) => {
                            // Only first alt in real_variants, call as 0/1
                            GenotypeTag::ref_het(NonZeroU8::new(vm as u8).unwrap_or(m))
                        }
                        (None, Some(vn)) => {
                            // Only second alt in real_variants, call as 0/2
                            GenotypeTag::ref_het(NonZeroU8::new(vn as u8).unwrap_or(n))
                        }
                        (None, None) => {
                            // Neither alt in real_variants, default to 0/0
                            GenotypeTag::HomRef
                        }
                    }
                }
            };

            let gt = Genotype::from(remapped_genotype);
            let gl = GenotypeLikelihood(smallvec_inline![Some(Phred::from(estimated.likelihood))]);
            let gc = GenotypeConfidence(smallvec_inline![Some(Phred::from(estimated.confidence))]);
            (gt, gl, gc)
        } else {
            // Fallback to homozygous reference (0/0) when estimate_genotype returns None
            (
                Genotype(smallvec![GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)]),
                GenotypeLikelihood(smallvec_inline![Some(Phred::from(Probability::ZERO))]),
                GenotypeConfidence(smallvec_inline![Some(Phred::from(Probability::ZERO))]),
            )
        };

        let has_ml = main_alts.iter().any(|alt| alt.filters.ml.is_some());

        let is_cpg = matches!(self.pos_metrics.cpg, super::InCpG::C | super::InCpG::G);
        // Original CpG context: either the position itself is a CpG in the reference,
        // or it is the matching partner of a de-novo CpG (denovo_adj flag).
        let is_orig_cpg = is_cpg || *self.pos_metrics.denovo_adj;
        // De-novo CpG context: any alt whose base + adjacent reference base creates a new
        // CpG dinucleotide. FormsDenovo already encodes the adjacency check, so we don't
        // need to replicate it here. We include all call types (not just RealVariant) so
        // that rejected records written with --all also receive a beta value.
        let is_denovo_cpg = self.alts.iter().any(|a| *a.metrics.denovo);
        let methylated = if self.pos_metrics.methylated.is_empty() && (is_orig_cpg || is_denovo_cpg)
        {
            let zero = CpgBeta {
                origin: CpgOrigin::Original,
                beta: Probability::ZERO,
                mod_count: 0,
                total_count: 0,
            };
            let mut betas = rastair_types::SmallVec::new();
            if is_orig_cpg {
                betas.push(zero);
            }
            if is_denovo_cpg {
                betas.push(CpgBeta { origin: CpgOrigin::DeNovo, ..zero });
            }
            Methylated(betas)
        } else {
            self.pos_metrics.methylated.clone()
        };

        let methylation_depth = MethylationDepth::from(&methylated);
        let methylation_alt_depth = MethylationAltDepth::from(&methylated);

        Format {
            genotype,
            genotype_likelihood,
            genotype_confidence,
            sample_read_depth: SampleReadDepth(self.pileup.reads.len()),
            methylated,
            methylation_depth,
            methylation_alt_depth,
            machine_learning_prediction: MachineLearningPrediction(if has_ml {
                main_alts.iter().map(|alt| *alt.filters.ml.unwrap_or_default()).collect()
            } else {
                smallvec![]
            }),
        }
    }
}

pub struct VcfRecordSet<'p> {
    pileup: &'p PileupMetrics,
    main: Record,
    rejected: SmallVec<Record, 2>,
    indel_records: SmallVec<Record, 1>,
}

impl<'p> VcfRecordSet<'p> {
    pub fn to_vec(&self, filters: &RecordFilters) -> SmallVec<&Record, 3> {
        let t = &self.pileup.tags;
        let cpg = t.cpg || t.denovo_cpg || t.denovo_cpg_partner;

        let mut v = match (filters.vcf_all, filters.cpgs_only) {
            (false, false) => {
                if t.covered {
                    smallvec![&self.main]
                } else {
                    smallvec![]
                }
            }
            (false, true) => {
                if t.covered && cpg {
                    smallvec![&self.main]
                } else {
                    smallvec![]
                }
            }
            (true, false) => {
                let mut v = smallvec![&self.main];
                v.extend(&self.rejected);
                v
            }
            (true, true) => {
                if cpg {
                    let mut v = smallvec![&self.main];
                    v.extend(&self.rejected);
                    v
                } else {
                    smallvec![]
                }
            }
        };

        // Indel records are always emitted when present (they already passed filters).
        v.extend(&self.indel_records);
        v
    }
}

fn build_indel_records(
    metrics: &PileupMetrics,
    ml_threshold: Option<Probability>,
) -> SmallVec<Record, 1> {
    metrics
        .indel_calls
        .iter()
        .map(|call| {
            let anchor: seqair_types::SmolStr = metrics.pileup.reference_base.into();

            let (ref_allele, alt_allele) = match &call.allele {
                IndelAllele::Insertion(bases) => {
                    let mut alt = String::with_capacity(1 + bases.len());
                    alt.push_str(&anchor);
                    for b in bases {
                        alt.push_str(b.as_str());
                    }
                    (anchor, alt.into())
                }
                IndelAllele::Deletion(bases) => {
                    let mut refr = String::with_capacity(1 + bases.len());
                    refr.push_str(&anchor);
                    for b in bases {
                        refr.push_str(b.as_str());
                    }
                    (refr.into(), anchor)
                }
            };

            let qual = call
                .ml
                .map(|ml| Phred::from(ml.inverted()).as_int())
                .or_else(|| Some(call.quality.as_int()));

            let main = VcfFixedFields {
                chrom: metrics.pileup.contig(),
                pos: metrics.pileup.pos,
                id: default(),
                r#ref: ref_allele,
                alt: smallvec![alt_allele],
                qual,
            };

            let genotype = Genotype::from(call.genotype);

            let mut filters = super::Filters::default();
            let ml_below_threshold = call.ml.zip(ml_threshold).is_some_and(|(ml, t)| ml < t);
            if ml_below_threshold {
                filters.add(low_ml_score.filter());
            } else {
                filters.add(PASS.filter());
            }

            Record {
                main,
                filters,
                info: Info {
                    read_depth: ReadDepth(call.depth as usize),
                    allele_read_depth: AlleleReadDepth(smallvec![
                        (call.depth.saturating_sub(call.alt_count)) as usize,
                        call.alt_count as usize,
                    ]),
                    ..default()
                },
                samples: smallvec_inline![Format {
                    genotype,
                    sample_read_depth: SampleReadDepth(call.depth as usize),
                    machine_learning_prediction: MachineLearningPrediction(
                        if let Some(ml) = call.ml { smallvec![*ml] } else { smallvec![] }
                    ),
                    ..default()
                }],
            }
        })
        .collect()
}
