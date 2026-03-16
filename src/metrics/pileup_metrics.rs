use crate::{
    call::{
        pileup::{Pileup, SimpleRead},
        variant_calling::EstimatedGenotype,
    },
    metrics::MethylationEvidenceStrandInfo,
    utils::{ByStrand, IntoF64, default, logging::ThisIsABug},
    vcf::{InCpG, Methylated},
};
use better_default::Default;
use color_eyre::{
    Result,
    eyre::{Context, bail},
};
use rastair_types::SmallVec;
use rastair_types::SmolStr;
use rastair_types::{Base, Probability, RmsAccumulator, RootMeanSquare, Strand};
use rastair_vcf::VcfFilter;
use std::ops::Deref;
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
    /// "Tags" for this positions, which will become calls
    pub tags: RecordTags,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct RecordTags {
    // The tags have been calculated. Mainly to debug :)
    pub set: bool,
    /// This position has coverage
    pub covered: bool,
    /// This is in a CpG site
    pub cpg: bool,
    /// This is a de-novo CpG site (not the partner)
    pub denovo_cpg: bool,
    /// This is the partner position of a de-novo CpG site
    pub denovo_cpg_partner: bool,
    /// This position is a variant (but not a denovo CpG)
    pub variant: bool,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Alt {
    pub base: Base,
    pub metrics: AlleleMetrics,
    pub filters: AltFilters,
    pub call: AltCall,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum AltCall {
    #[default]
    Uncalled,
    RealVariant,
    MethylationEvidenceOnly {
        for_base: Base,
    },
    ReadError,
}

impl PileupMetrics {
    /// Create new metrics from pileup
    ///
    /// NOTE: The extended metrics in `PositionMetrics` are not set here and
    /// need to be set later using `set_extended_metrics`.
    pub fn new(pileup: Pileup) -> Result<Self> {
        let ref_base = pileup.reference_base;
        let total_reads = pileup.reads.len();

        // Single pass: accumulate per-base metrics, position-level RMS, and discover alleles
        let mut accumulators = PerBaseAccumulators::default();
        let mut pos_baseq = RmsAccumulator::new();
        let mut pos_mapq = RmsAccumulator::new();
        let mut mapq0: u32 = 0;
        let mut alt_bases: SmallVec<Base, 4> = SmallVec::new();
        for read in pileup.reads.iter() {
            accumulators.accumulate(read);
            pos_baseq.add(f64::from(read.qual));
            pos_mapq.add(f64::from(read.mapq));
            if read.mapq == 0 {
                mapq0 += 1;
            }
            if read.base.known_index().is_some()
                && read.base != ref_base
                && !alt_bases.contains(&read.base)
            {
                alt_bases.push(read.base);
            }
        }

        trace!(pos = pileup.pos, ?ref_base, ?alt_bases, "New pileup");

        let pos_metrics = PositionMetrics::from_pileup(
            &pileup,
            PositionMetricsExt::default(),
            pos_baseq.finish(),
            pos_mapq.finish(),
            mapq0,
        );

        // Reference base can be Unknown (N in FASTA) — no accumulator slot exists for it,
        // and no reads will ever match it, so just use default metrics.
        let ref_metrics = if let Some(acc) = accumulators.take(ref_base) {
            acc.finish(ref_base, total_reads, &pileup)
                .wrap_err("Failed to compute allele metrics for reference")?
        } else {
            AlleleMetrics { base: ref_base, ..default() }
        };

        let alts = alt_bases
            .iter()
            .map(|&base| {
                // alt_bases only contains known bases (filtered above), so take always succeeds
                let acc = accumulators
                    .take(base)
                    .ok_or_else(|| color_eyre::eyre::eyre!("unknown base {base} in alt_bases"))?;
                let metrics = acc
                    .finish(base, total_reads, &pileup)
                    .wrap_err("Failed to compute allele metrics for alt")?;
                Ok(Alt { base, metrics, filters: AltFilters::default(), call: default() })
            })
            .collect::<Result<_>>()?;

        Ok(PileupMetrics {
            pileup,
            pos_metrics,
            pos_filters: Filters::default(),
            ref_metrics,
            alts,
            tags: RecordTags::default(),
        })
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

    pub fn contig_pos(&self) -> SmolStr {
        use std::fmt::Write as _;
        let mut res = rastair_types::smol_str::SmolStrBuilder::new();

        write!(&mut res, "{}:{}", self.contig(), self.pos()).expect("works");

        res.finish()
    }

    pub fn alt(&self, alt: Base) -> Option<&AlleleMetrics> {
        self.alts.iter().find(|a| a.base == alt).map(|a| &a.metrics)
    }

    pub fn allele(&self, base: Base) -> Option<&AlleleMetrics> {
        if base == self.ref_base() {
            Some(&self.ref_metrics)
        } else {
            self.alts.iter().find(|a| a.base == base).map(|a| &a.metrics)
        }
    }

    pub fn alt_metrics(&self, alt: Base) -> Option<MetricsForAlt<'_>> {
        let alt = self.alts.iter().find(|a| a.base == alt);
        alt.map(|alt| MetricsForAlt { metrics: self, alt: &alt.metrics })
    }

    pub fn alt_filters(&self, alt: Base) -> Option<&AltFilters> {
        self.alts.iter().find(|a| a.base == alt).map(|a| &a.filters)
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
        if self.pos_filters.other_pos_in_denovo_passes {
            return true;
        }
        self.pos_filters.pass() && self.alts.iter().any(|a| a.filters.pass(ml_threshold))
    }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[cfg_attr(test, derive(Default))]
pub struct PositionMetrics {
    /// Read depth, i.e., number of reads covering this position
    pub depth: u32,
    /// Base quality
    pub baseq: RootMeanSquare,
    /// Mapping quality
    pub mapq: RootMeanSquare,
    /// Number of reads with mapping quality 0
    pub mapq0: u32,
    /// Is this position in a CpG context in the reference?
    pub cpg: InCpG,

    /// Extended metrics
    // set by `call` later since they depend on more context
    #[serde(flatten)]
    pub extended: PositionMetricsExt,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct PositionMetricsExt {
    /// Entropy of the surrounding region
    pub region_entropy: f64,
    /// Estimated genotype
    pub genotype: Option<EstimatedGenotype>,
    /// Methylation strand info
    pub methylation_strand_info: MethylationEvidenceStrandInfo,
    /// Methylation beta
    pub methylated: Methylated,
    /// Is this position a de-novo cpg candidate?
    pub denovo_adj: DenovoAdjecent,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
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
    pub fn from_pileup(
        pileup: &Pileup,
        extended: PositionMetricsExt,
        baseq: RootMeanSquare,
        mapq: RootMeanSquare,
        mapq0: u32,
    ) -> Self {
        PositionMetrics {
            depth: u32::try_from(pileup.reads.len()).expect("depth fits into u32"),
            baseq,
            mapq,
            mapq0,
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

impl FormsDenovo {
    pub fn some(&self) -> Option<Self> {
        match self {
            FormsDenovo::No => None,
            _ => Some(*self),
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
        if self.filters.other_pos_in_denovo_passes {
            return true;
        }
        if let Some(ml_threshold) = ml_threshold
            && let Some(ml) = self.ml
        {
            ml >= ml_threshold
        } else {
            self.filters.is_empty()
        }
    }
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct Filters {
    pub other_pos_in_denovo_passes: bool,
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
        self.other_pos_in_denovo_passes || self.filters.is_empty()
    }
}

impl Deref for Filters {
    type Target = SmallVec<SmolStr, 6>;

    fn deref(&self) -> &Self::Target {
        &self.filters
    }
}

#[derive(Debug, Clone, Default)]
struct AlleleAccumulator {
    depth: u32,
    baseq: RmsAccumulator,
    mapq: RmsAccumulator,
    baseq_ot: RmsAccumulator,
    baseq_ob: RmsAccumulator,
    mapq_ot: RmsAccumulator,
    mapq_ob: RmsAccumulator,
    aligned: RmsAccumulator,
    indels: RmsAccumulator,
    pos_in_read: RmsAccumulator,
    ot_count: u32,
    ob_count: u32,
}

impl AlleleAccumulator {
    fn add(&mut self, read: &SimpleRead) {
        self.depth += 1;
        let q = f64::from(read.qual);
        let m = f64::from(read.mapq);
        self.baseq.add(q);
        self.mapq.add(m);
        match read.strand {
            Strand::OT => {
                self.ot_count += 1;
                self.baseq_ot.add(q);
                self.mapq_ot.add(m);
            }
            Strand::OB => {
                self.ob_count += 1;
                self.baseq_ob.add(q);
                self.mapq_ob.add(m);
            }
            Strand::Unknown => {}
        }
        self.aligned.add(f64::from(read.matching_bases));
        self.indels.add(f64::from(read.indels));
        self.pos_in_read.add(f64::from(read.position.pos) / f64::from(read.position.read_length));
    }

    fn finish(self, base: Base, total_reads: usize, pileup: &Pileup) -> Result<AlleleMetrics> {
        use Base::*;

        if self.depth == 0 {
            // Can happen at canonical CpG sites with no evidence for a particular allele
            trace!(
                pos = pileup.pos,
                ref_base = ?pileup.reference_base,
                ?base,
                pileup_reads = pileup.reads.len(),
                "No reads for allele"
            );
            return Ok(AlleleMetrics { base, ..default() });
        }

        // Should be impossible: depth > 0 implies total_reads > 0 since reads were counted
        // from the same pileup. Guard defensively since a division by zero here would produce
        // NaN/Inf which Probability::new rejects anyway, but better to fail with a clear message.
        if total_reads == 0 {
            bail!("allele has depth {} but pileup has 0 total reads — this is a bug", self.depth);
        }

        let denovo = if base == pileup.reference_base {
            FormsDenovo::No
        } else if pileup.ref_before() == C && base == G {
            FormsDenovo::ThisBecomesG
        } else if pileup.ref_after() == G && base == C {
            FormsDenovo::ThisBecomesC
        } else {
            FormsDenovo::No
        };

        Ok(AlleleMetrics {
            base,
            depth: self.depth,
            baseq: self.baseq.finish(),
            mapq: self.mapq.finish(),
            strand_count: ByStrand { base, ot: self.ot_count, ob: self.ob_count },
            baseq_s: ByStrand { base, ot: self.baseq_ot.finish(), ob: self.baseq_ob.finish() },
            mapq_s: ByStrand { base, ot: self.mapq_ot.finish(), ob: self.mapq_ob.finish() },
            num_aligned_bases: self.aligned.finish(),
            num_indels: self.indels.finish(),
            position_in_read: self.pos_in_read.finish(),
            allele_frequency: Probability::new(self.depth.f() / total_reads.f())
                .wrap_err("allele frequency not in [0,1]")
                .this_is_a_bug()?,
            denovo,
        })
    }
}

/// Per-base accumulators indexed by [`Base::known_index`], one slot per `Base::KNOWN`.
#[derive(Debug, Default)]
struct PerBaseAccumulators([AlleleAccumulator; 4]);

impl PerBaseAccumulators {
    fn accumulate(&mut self, read: &SimpleRead) {
        let Some(idx) = read.base.known_index() else { return };
        self.0[idx].add(read);
    }

    fn take(&mut self, base: Base) -> Option<AlleleAccumulator> {
        let idx = base.known_index()?;
        Some(std::mem::take(&mut self.0[idx]))
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
