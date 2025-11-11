use crate::{
    call::{
        pileup::{Pileup, SimpleRead},
        variant_calling::EstimatedGenotype,
    },
    utils::{ByStrand, IntoF64, default, logging::ThisIsABug},
    vcf::{InCpG, Methylated},
};
use better_default::Default;
use color_eyre::{
    Result,
    eyre::{Context, bail},
};
use rastair_types::{Base, Probability, RootMeanSquare, Strand, rms::RootMeanSquareExt};
use rastair_vcf::VcfFilter;
use smallvec::SmallVec;
use smol_str::SmolStr;
use std::ops::Deref;
use thiserror::Error;
use tracing::trace;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PileupMetrics {
    /// The underlying pileup
    pub pileup: Pileup,
    /// Metrics about the position itself
    pub pos_metrics: PositionMetrics,
    /// Filters that apply to the entire pileup
    pub pos_filters: Filters,
    /// Metrics for the reference allele
    pub ref_metrics: AlleleMetrics,
    /// Metrics and filters for each alternative allele
    pub alts: SmallVec<Alt, 2>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Alt {
    pub base: Base,
    pub metrics: AlleleMetrics,
    pub filters: AltFilters,
}

impl PileupMetrics {
    /// Create new metrics from pileup
    ///
    /// NOTE: The extended metrics in `PositionMetrics` are not set here and
    /// need to be set later using `set_extended_metrics`.
    pub fn new(pileup: Pileup) -> Result<Self> {
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

        // Compute initial metrics but keep extended empty to set later
        let pos_metrics = PositionMetrics::from_pileup(&pileup, PositionMetricsExt::default());

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
            .collect::<Result<_>>()
            .wrap_err("Failed to compute allele metrics for alt alleles")?;

        drop(by_allele);

        let pos_filters = Filters::default();

        let metrics = PileupMetrics { pileup, pos_metrics, pos_filters, ref_metrics, alts };

        Ok(metrics)
    }

    pub fn set_extended_metrics(&mut self, ext: PositionMetricsExt) {
        self.pos_metrics.extended = ext;
    }

    /// Get reference base
    pub fn ref_base(&self) -> Base {
        self.pileup.reference_base
    }

    pub fn contig(&self) -> SmolStr {
        self.pileup.region.contig.clone()
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

    pub fn forms_denovo(&self) -> bool {
        *self.pos_metrics.denovo_adj || self.alts.iter().any(|a| *a.metrics.denovo)
    }

    pub fn pass(&self, ml_threshold: Option<Probability>) -> bool {
        if self.pos_filters.other_pos_in_cpg_passes {
            return true;
        }
        self.pos_filters.pass() && self.alts.iter().all(|a| a.filters.pass(ml_threshold))
    }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[cfg_attr(test, derive(Default))]
pub struct PositionMetrics {
    /// Read depth, i.e., number of reads covering this position
    pub depth: usize,
    /// Base quality
    pub baseq: RootMeanSquare,
    /// Mapping quality
    pub mapq: RootMeanSquare,
    /// Number of reads with mapping quality 0
    pub mapq0: usize,
    /// Is this position in a CpG context in the reference?
    pub cpg: InCpG,

    /// Extended metrics
    // set by `call` later since they depend on more context
    // todo: explore using type-state for this
    #[serde(flatten)]
    extended: PositionMetricsExt,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct PositionMetricsExt {
    /// Entropy of the surrounding region
    pub region_entropy: f64,
    /// Estimated genotype
    pub genotype: Option<EstimatedGenotype>,
    /// Methylation beta
    pub methylated: Methylated,
    /// Is this position a de-novo cpg candidate?
    pub denovo_adj: DenovoAdjecent,
}

#[derive(Debug, Clone, Copy, Default, serde::Serialize, serde::Deserialize)]
pub enum DenovoAdjecent {
    #[default]
    No,
    ThisIsTheMatchingC,
    ThisIsTheMatchingG,
}

impl Deref for DenovoAdjecent {
    type Target = bool;

    fn deref(&self) -> &Self::Target {
        match self {
            DenovoAdjecent::No => &false,
            _ => &true,
        }
    }
}

impl PositionMetrics {
    pub fn from_pileup(pileup: &Pileup, extended: PositionMetricsExt) -> Self {
        PositionMetrics {
            depth: pileup.reads.len(),
            baseq: pileup.reads.iter().map(|x| x.qual).collect(),
            mapq: pileup.reads.iter().map(|x| x.mapq).collect(),
            mapq0: pileup.reads.iter().filter(|x| x.mapq == 0).count(),
            cpg: InCpG::from(pileup),

            // These fields are given by all
            extended,
        }
    }
}

impl Deref for PositionMetrics {
    type Target = PositionMetricsExt;

    fn deref(&self) -> &Self::Target {
        &self.extended
    }
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct AlleleMetrics {
    pub base: Base,
    /// Read depth, i.e. number of reads supporting this allele
    pub depth: u32,
    /// base quality for the allele
    pub baseq: RootMeanSquare,
    /// mapping quality for the allele
    pub mapq: RootMeanSquare,
    /// count of reads by strand, also known as strand bias
    pub strand_count: ByStrand<u32>,
    /// base quality by strand
    pub baseq_s: ByStrand<RootMeanSquare>,
    /// mapping quality by strand
    pub mapq_s: ByStrand<RootMeanSquare>,
    /// number of aligned bases in read
    pub num_aligned_bases: RootMeanSquare,
    /// number of indels in read
    pub num_indels: RootMeanSquare,
    /// relative position in read
    pub position_in_read: RootMeanSquare,
    /// Allele frequency
    pub allele_frequency: Probability,
    /// does this alt form a de-novo cpg?
    pub denovo: FormsDenovo,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
pub enum FormsDenovo {
    #[default]
    No,
    ThisBecomesC,
    ThisBecomesG,
}

impl Deref for FormsDenovo {
    type Target = bool;

    fn deref(&self) -> &Self::Target {
        match self {
            FormsDenovo::No => &false,
            _ => &true,
        }
    }
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct AltFilters {
    /// ML prediction: probability this is a true variant
    pub ml: Option<Probability>,
    pub filters: Filters,
}

impl AltFilters {
    pub fn pass(&self, ml_threshold: Option<Probability>) -> bool {
        if self.filters.other_pos_in_cpg_passes {
            return true;
        }
        if ml_threshold.is_some() { self.ml > ml_threshold } else { self.filters.is_empty() }
    }
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct Filters {
    pub other_pos_in_cpg_passes: bool,
    filters: SmallVec<SmolStr, 6>,
}

impl Filters {
    pub fn add(&mut self, filter: impl VcfFilter, condition: impl FnOnce() -> bool) {
        if condition() {
            self.filters.push(filter.filter());
        }
    }

    pub fn merge(&mut self, other: Filters) {
        for filter in other.filters {
            if !self.filters.contains(&filter) {
                self.filters.push(filter);
            }
        }
    }

    pub fn pass(&self) -> bool {
        self.other_pos_in_cpg_passes || self.filters.is_empty()
    }
}

impl Deref for Filters {
    type Target = SmallVec<SmolStr, 6>;

    fn deref(&self) -> &Self::Target {
        &self.filters
    }
}

impl AlleleMetrics {
    pub fn from_bases(reads: &[&SimpleRead], pileup: &Pileup) -> Result<Self> {
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

            return Ok(AlleleMetrics { base: pileup.reference_base, ..default() });
        }

        let base = reads[0].base;

        let i = || reads.iter();
        let ot = || i().filter(|x| x.strand == OT);
        let ob = || i().filter(|x| x.strand == OB);

        let denovo = {
            if base == pileup.reference_base {
                // we are looking at the ref allele
                FormsDenovo::No
            } else if pileup.ref_before() == C && base == G {
                FormsDenovo::ThisBecomesG
            } else if pileup.ref_after() == G && base == C {
                FormsDenovo::ThisBecomesC
            } else {
                FormsDenovo::No
            }
        };

        Ok(AlleleMetrics {
            base,
            depth: allele_depth,
            baseq: i().map(|x| x.qual).rms(),
            mapq: i().map(|x| x.mapq).rms(),
            strand_count: ByStrand {
                base,
                ot: u32::try_from(ot().count()).wrap_err("read count should fit into u32")?,
                ob: u32::try_from(ob().count()).wrap_err("read count should fit into u32")?,
            },
            baseq_s: ByStrand {
                base,
                ot: ot().map(|x| x.qual).collect(),
                ob: ob().map(|x| x.qual).collect(),
            },
            mapq_s: ByStrand {
                base,
                ot: ot().map(|x| x.mapq).collect(),
                ob: ob().map(|x| x.mapq).collect(),
            },
            num_aligned_bases: i().map(|x| x.matching_bases).rms(),
            num_indels: i().map(|x| x.indels).rms(),
            position_in_read: i()
                .map(|b| f64::from(b.position.pos) / f64::from(b.position.read_length))
                .rms(),
            allele_frequency: Probability::new(allele_depth.f() / pileup.reads.len().f())
                .wrap_err("allele frequency not in in [0,1]")
                .this_is_a_bug()?,
            denovo,
        })
    }
}

pub struct MetricsForAlt<'p> {
    pub metrics: &'p PileupMetrics,
    pub alt: &'p AlleleMetrics,
}

impl MetricsForAlt<'_> {
    pub fn is_evidence_for_methylation(&self) -> bool {
        (self.metrics.pos_metrics.cpg == InCpG::C && self.alt.base == Base::T)
            || (self.metrics.pos_metrics.cpg == InCpG::G && self.alt.base == Base::A)
    }
}
