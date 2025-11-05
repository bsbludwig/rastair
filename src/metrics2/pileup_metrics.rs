use crate::{
    call::{
        variant_calling::EstimatedGenotype,
        variants::{SimpleRead, VariantCandidatePileup},
    },
    utils::ByStrand,
    vcf::{DeNovoCpGCandidate, InCpG, Methylated, SequenceContext},
};
use better_default::Default;
use color_eyre::{
    Result,
    eyre::{Context, bail},
};
use rastair_types::{Base, Probability, RootMeanSquare, Strand, rms::RootMeanSquareExt};
use smallvec::SmallVec;
use smol_str::SmolStr;
use thiserror::Error;
use tracing::trace;

pub struct PileupMetrics {
    /// The underlying pileup
    pub pileup: VariantCandidatePileup,
    /// Metrics about the position itself
    pub pos_metrics: PositionMetrics,
    /// Filters that apply to the entire pileup
    pub pos_filters: SmallVec<SmolStr, 5>,
    /// Metrics for the reference allele
    pub ref_metrics: AlleleMetrics,
    /// Metrics and filters for each alternative allele
    pub alts: SmallVec<Alt, 2>,
}

pub struct Alt {
    pub base: Base,
    pub metrics: AlleleMetrics,
    pub filters: AltFilters,
}

impl TryFrom<VariantCandidatePileup> for PileupMetrics {
    type Error = color_eyre::eyre::Report;

    fn try_from(pileup: VariantCandidatePileup) -> Result<Self, Self::Error> {
        let by_allele = pileup.by_allele();
        let [reference, alts_reads @ ..] = by_allele.as_slice() else {
            // todo: do we reach this point also when we have no reads covering
            // this region but there is a cpg in the ref?
            bail!("No alleles found in pileup");
        };
        trace!(
            pos = pileup.pos,
            ref_base = ?reference.base,
            alt_bases = ?alts_reads.iter().map(|b| b.base).collect::<Vec<_>>(),
            "Found candidate"
        );

        let ref_base = reference.base;

        #[derive(Debug, Error)]
        #[error("Failed to compute allele metric for {0}")]
        struct AlleleMetricError(Base);

        let pos_metrics = PositionMetrics::from_pileup(&pileup);
        let ref_metrics = if reference.is_empty() {
            // this can happen at canonical cpg sites with no evidence
            AlleleMetrics { base: ref_base, ..AlleleMetrics::default() }
        } else {
            AlleleMetrics::from_bases(reference, &pileup).wrap_err(AlleleMetricError(ref_base))?
        };

        let alts = alts_reads
            .iter()
            .map(|pile| {
                let base = pile.base;
                AlleleMetrics::from_bases(pile, &pileup)
                    .wrap_err(AlleleMetricError(base))
                    .map(|metrics| Alt { base, metrics, filters: AltFilters::default() })
            })
            .collect::<Result<SmallVec<_, 2>>>()
            .wrap_err("Failed to compute allele metrics for alt alleles")?;

        drop(by_allele);

        let pos_filters = SmallVec::new();

        Ok(PileupMetrics { pileup, pos_metrics, pos_filters, ref_metrics, alts })
    }
}

impl PileupMetrics {
    /// Get reference base
    pub fn ref_base(&self) -> Base {
        self.pileup.reference_base
    }

    pub fn contig(&self) -> SmolStr {
        self.pileup.segment.contig.clone()
    }

    pub fn pos(&self) -> u32 {
        self.pileup.pos
    }

    pub fn alt(&self, alt: Base) -> Option<&AlleleMetrics> {
        self.alts.iter().find(|a| a.base == alt).map(|a| &a.metrics)
    }

    pub fn alt_metrics(&self, alt: Base) -> Option<MetricsForAlt<'_>> {
        let alt = self.alts.iter().find(|a| a.base == alt);
        alt.map(|alt| MetricsForAlt { metrics: self, alt: &alt.metrics })
    }

    pub fn alt_filters_mut(&mut self, alt: Base) -> Option<&mut AltFilters> {
        self.alts.iter_mut().find(|a| a.base == alt).map(|a| &mut a.filters)
    }

    /// Get all alternative bases in the pileup (for lookup with mutation)
    pub fn alts(&self) -> SmallVec<Base, 4> {
        self.alts.iter().map(|a| a.base).collect()
    }

    pub fn alts_metrics(&self) -> impl Iterator<Item = &AlleleMetrics> {
        self.alts.iter().map(|a| &a.metrics)
    }

    pub fn ref_alts_metrics(&self) -> impl Iterator<Item = &AlleleMetrics> {
        std::iter::once(&self.ref_metrics).chain(self.alts.iter().map(|a| &a.metrics))
    }
}

pub struct PositionMetrics {
    pub read_depth: usize,
    /// Base quality
    pub baseq: RootMeanSquare,
    /// Mapping quality
    pub mapq: RootMeanSquare,
    /// Number of reads with mapping quality 0
    pub mapq0: usize,
    /// Entropy of the region around the position
    pub region_entropy: f64,
    /// Sequence context around the position in the reference
    pub sequence_context: SequenceContext,
    pub cpg: InCpG,
    /// Is this position a de-novo cpg candidate?
    pub de_novo_cpg_candidate: DeNovoCpGCandidate,
    pub genotype: Option<EstimatedGenotype>,
    pub methylated: Methylated,
}

impl PositionMetrics {
    pub fn from_pileup(pileup: &VariantCandidatePileup) -> Self {
        PositionMetrics {
            read_depth: pileup.reads.len(),
            baseq: pileup.reads.iter().map(|x| x.qual).collect::<RootMeanSquare>(),
            mapq: pileup.reads.iter().map(|x| x.mapq).collect::<RootMeanSquare>(),
            mapq0: pileup.reads.iter().filter(|x| x.mapq == 0).count(),
            region_entropy: pileup.entropy(),
            sequence_context: SequenceContext::from(pileup),
            cpg: InCpG::from(pileup),
            de_novo_cpg_candidate: DeNovoCpGCandidate::from(pileup),
            // this is set later in `call`
            genotype: None,
            methylated: Methylated::Unknown,
        }
    }
}

#[derive(Debug, Default)]
pub struct AlleleMetrics {
    pub base: Base,
    pub depth: u32,
    /// base quality for the allele
    pub baseq: RootMeanSquare,
    /// mapping quality for the allele
    pub mapq: RootMeanSquare,
    pub strand_count: ByStrand<u32>,
    pub baseq_s: ByStrand<RootMeanSquare>,
    pub mapq_s: ByStrand<RootMeanSquare>,
    /// number of aligned bases in read
    pub num_aligned_bases: RootMeanSquare,
    /// number of indels in read
    pub num_indels: RootMeanSquare,
    /// relative position in read
    pub position_in_read: RootMeanSquare,
    /// does this alt form a de-novo cpg?
    pub denovo: DeNovoCpGCandidate,
}

#[derive(Default)]
pub struct AltFilters {
    /// ML prediction: probability this is a true variant
    pub ml: Option<Probability>,

    // TODO: Turn this into a `struct` that can be const-constructed
    pub filters: SmallVec<SmolStr, 5>,
}

impl AlleleMetrics {
    pub fn from_bases(reads: &[&SimpleRead], pileup: &VariantCandidatePileup) -> Result<Self> {
        use Base::*;
        use Strand::*;

        let allele_depth = u32::try_from(reads.len()).wrap_err("read count fits into u32")?;
        if allele_depth == 0 {
            trace!(
                %pileup.pos,
                base=%pileup.reference_base,
                pileup_reads=?pileup.reads.len(),
                allele_reads=?reads.len(),
                "Why are we here? No reads for allele metrics calculation"
            );

            return Ok(AlleleMetrics { base: pileup.reference_base, ..Default::default() });
        }

        let base = reads[0].base;

        let i = || reads.iter();
        let ot = || i().filter(|x| x.strand == OT);
        let ob = || i().filter(|x| x.strand == OB);

        let denovo = {
            if base == pileup.reference_base {
                DeNovoCpGCandidate::NotCandidate
            } else if pileup.ref_before() == C && base == G {
                DeNovoCpGCandidate::Candidate {
                    ref_base: base,
                    alt_base: G,
                    alt_index: pileup.alts().iter().position(|b| *b == base).unwrap_or(42),
                }
            } else if pileup.ref_after() == G && base == C {
                DeNovoCpGCandidate::Candidate {
                    ref_base: base,
                    alt_base: C,
                    alt_index: pileup.alts().iter().position(|b| *b == base).unwrap_or(42),
                }
            } else {
                DeNovoCpGCandidate::NotCandidate
            }
        };

        Ok(AlleleMetrics {
            base,
            depth: allele_depth,
            baseq: i().map(|x| x.qual).rms(),
            mapq: i().map(|x| x.mapq).rms(),
            strand_count: ByStrand {
                base,
                ot: u32::try_from(ot().count()).wrap_err("read count fits into u32")?,
                ob: u32::try_from(ob().count()).wrap_err("read count fits into u32")?,
            },
            baseq_s: ByStrand {
                base,
                ot: ot().map(|x| x.qual).collect::<RootMeanSquare>(),
                ob: ob().map(|x| x.qual).collect::<RootMeanSquare>(),
            },
            mapq_s: ByStrand {
                base,
                ot: ot().map(|x| x.mapq).collect::<RootMeanSquare>(),
                ob: ob().map(|x| x.mapq).collect::<RootMeanSquare>(),
            },
            num_aligned_bases: i().map(|x| x.matching_bases).rms(),
            num_indels: i().map(|x| x.indels).rms(),
            position_in_read: i()
                .map(|b| f64::from(b.position.pos) / f64::from(b.position.read_length))
                .rms(),
            denovo,
        })
    }
}

pub struct MetricsForAlt<'p> {
    pub metrics: &'p PileupMetrics,
    pub alt: &'p AlleleMetrics,
}
