use super::{Methylated, Record};
use crate::{
    call::{RecordFilters, pileup::Pileup},
    metrics::{AlleleMetrics, Alt, Filters, PileupMetrics, PositionMetrics},
    utils::{IntoF64 as _, default},
    vcf::{
        AlleleBaseQuality, AlleleMapQuality, AlleleSpecificStrandBias, DeNovoCpGCandidate, Entropy,
        Filters as VcfFilter, Format, GenotypeConfidence, GenotypeLikelihood, Info,
        MachineLearningPrediction, NumAlignedBases, NumIndels, PositionInRead,
        StrandSpecificBaseQuality, StrandSpecificMappingQuality, dnCpG_adj,
    },
};
use color_eyre::Result;
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

pub struct VcfRecordSet {
    reference_line: Option<Record>,
    alt_lines: Vec<Record>,
    is_cpg_site: bool,
}

impl VcfRecordSet {
    pub fn iter(self) -> impl Iterator<Item = Record> {
        // If there is a reference line, yield it first
        // Then yield all alt lines
        self.reference_line.into_iter().chain(self.alt_lines)
    }

    pub fn filter(&mut self, filters: &RecordFilters) {
        match (filters.vcf_all, filters.cpgs_only) {
            (false, false) => {
                // default behavior: only passing records with alts
                self.reference_line = self.reference_line.take().filter(|r| r.filters.pass());
                self.alt_lines.retain(|r| r.filters.pass());
            }
            (false, true) => {
                self.reference_line = self.reference_line.take().filter(|r| r.filters.pass());
                // passing CpGs (in ref line, untouched) and passing de-novo CpG candidates
                self.alt_lines.retain(|r| *r.info.de_novo_cp_g_candidate && r.filters.pass());
            }
            (true, false) => {
                // `--all` means all, do nothing
            }
            (true, true) => {
                // `--all -c` means all CpGs and de-novo CpG candidates
                self.alt_lines.retain(|r| *r.info.de_novo_cp_g_candidate);
            }
        }
    }
}

impl PileupMetrics {
    /// Convert the metrics to VCF records
    ///
    /// Generates multiple VCF rows from a single pileup position:
    /// - Writes ref->. when there are no alts or when methylation evidence exists (PASS)
    /// - Writes separate ref->alt rows for each alt allele with their respective filters
    /// - Methylation evidence = C->T or G->A transitions with low ML scores
    pub fn to_vcf_records(&self, ml_threshold: Option<Probability>) -> Result<VcfRecordSet> {
        let mut res = VcfRecordSet {
            reference_line: None,
            alt_lines: Vec::new(),
            is_cpg_site: *self.pos_metrics.cpg,
        };
        let pass_filters = {
            let mut filters = Filters::default();
            filters.add(PASS, || true);
            filters
        };

        // Write ref->. (no alt) row when:
        if *self.pos_metrics.cpg {
            // 1. There might be methylation evidence (the ref base is methylated)
            let row = MetricsSubset {
                pileup: &self.pileup,
                pos_metrics: &self.pos_metrics,
                // This row is PASS since it represent the reference allele.
                pos_filters: &pass_filters,
                ref_metrics: &self.ref_metrics,
                alts: smallvec![],
                is_ref_only_row: true,
            };
            let row = row.to_vcf_row(ml_threshold);
            res.reference_line = Some(row);
        } else if *self.pos_metrics.denovo_adj {
            // 2. This is the matching position of a passing de-novo CpG
            let row = MetricsSubset {
                pileup: &self.pileup,
                pos_metrics: &self.pos_metrics,
                // This row is PASS when the other position in the de-novo CpG passes filters.
                pos_filters: &{
                    let mut filters = self.pos_filters.clone();
                    filters.other_pos_in_denovo_passes =
                        self.pos_filters.other_pos_in_denovo_passes;
                    filters.add(dnCpG_adj, || !self.pos_filters.other_pos_in_denovo_passes);
                    filters
                },
                ref_metrics: &self.ref_metrics,
                alts: smallvec![],
                is_ref_only_row: true,
            };
            let row = row.to_vcf_row(ml_threshold);
            res.reference_line = Some(row);
        } else if self.alts.is_empty() {
            // 3. There are no alts at this position
            let row = MetricsSubset {
                pileup: &self.pileup,
                pos_metrics: &self.pos_metrics,
                // This row is PASS since it represent the reference allele.
                pos_filters: &pass_filters,
                ref_metrics: &self.ref_metrics,
                alts: smallvec![],
                is_ref_only_row: true,
            };
            let row = row.to_vcf_row(ml_threshold);
            res.reference_line = Some(row);
        }

        // Write separate ref->alt row for each alt allele.
        // Each alt gets its own row with its specific filters applied.
        for alt in &self.alts {
            let row = MetricsSubset {
                pileup: &self.pileup,
                pos_metrics: &self.pos_metrics,
                pos_filters: &self.pos_filters,
                ref_metrics: &self.ref_metrics,
                alts: smallvec![alt],
                is_ref_only_row: false,
            };
            let row = row.to_vcf_row(ml_threshold);
            res.alt_lines.push(row);
        }

        Ok(res)
    }
}

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
        let (genotype, genotype_likelihood, genotype_confidence) =
            if let Some(estimate) = self.pos_metrics.genotype {
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
        } else {
            // For alt rows, no methylation info
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
