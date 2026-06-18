use crate::{
    call::{
        pileup::{Pileup, SimpleRead, indels},
        variant_calling::EstimatedGenotype,
    },
    metrics::{MethylationEvidenceStrandInfo, PairedCounts, ReadKey},
    sequence::ChunkRegion,
    utils::{ByStrand, IntoF64, SequenceContext, default, logging::ThisIsABug},
    vcf::{InCpG, Methylated, RastairFilter},
};
use better_default::Default;
use color_eyre::{
    Result,
    eyre::{Context, bail},
};
use seqair_types::SmallVec;
use seqair_types::SmolStr;
use seqair_types::{Base, Probability, RmsAccumulator, RootMeanSquare, Strand};
use std::ops::Deref;
use tracing::{trace, warn};

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PileupMetrics {
    pub region: ChunkRegion,
    pub pos: u32,
    pub reference_base: Base,
    pub context: SequenceContext,
    pub pos_metrics: PositionMetrics,
    pub pos_filters: Filters,
    pub ref_metrics: AlleleMetrics,
    pub alts: SmallVec<Alt, 2>,
    /// Counts of (`my_base`, `before_base`) pairs by strand
    pub before_counts: PairedCounts,
    /// Counts of (`my_base`, `after_base`) pairs by strand
    pub after_counts: PairedCounts,
    /// "Tags" for this positions, which will become calls
    pub tags: RecordTags,
    #[serde(default)]
    pub indel_data: Option<Box<indels::IndelData>>,
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
        let Pileup {
            region,
            pos,
            reference_base,
            context,
            indel_observations,
            homopolymer_run,
            dinucleotide_run,
            soft_clip_count,
            reads,
            noisy_ref_count,
            indel_ref_window,
            indel_ref_anchor,
        } = pileup;
        let mut accumulators = PerBaseAccumulators::default();
        let mut pos_baseq = RmsAccumulator::new();
        let mut pos_mapq = RmsAccumulator::new();
        let mut mapq0: u32 = 0;
        let mut alt_bases: SmallVec<Base, 4> = SmallVec::new();
        let mut total_reads: usize = 0;
        for read in reads.iter() {
            total_reads += 1;
            let qual_sq = f64::from(read.qual).powi(2);
            let mapq_sq = f64::from(read.mapq).powi(2);
            accumulators.accumulate(read, qual_sq, mapq_sq);
            pos_baseq.add_squared(qual_sq);
            pos_mapq.add_squared(mapq_sq);
            if read.mapq == 0 {
                mapq0 += 1;
            }
            if read.base.known_index().is_some()
                && read.base != reference_base
                && !alt_bases.contains(&read.base)
            {
                alt_bases.push(read.base);
            }
        }

        trace!(pos, ?reference_base, ?alt_bases, "New pileup");

        let pos_metrics = PositionMetrics::new(
            total_reads,
            reference_base,
            context.before_1,
            context.after_1,
            pos_baseq.finish(),
            pos_mapq.finish(),
            mapq0,
        );

        let ref_metrics = if let Some(acc) = accumulators.take(reference_base) {
            acc.finish(reference_base, total_reads, pos, reference_base, &context)
                .wrap_err("Failed to compute allele metrics for reference")?
        } else {
            AlleleMetrics { base: reference_base, ..default() }
        };

        let alts = alt_bases
            .iter()
            .map(|&base| {
                let acc = accumulators
                    .take(base)
                    .ok_or_else(|| color_eyre::eyre::eyre!("unknown base {base} in alt_bases"))?;
                let metrics = acc
                    .finish(base, total_reads, pos, reference_base, &context)
                    .wrap_err("Failed to compute allele metrics for alt")?;
                Ok(Alt { base, metrics, filters: AltFilters::default(), call: default() })
            })
            .collect::<Result<_>>()?;

        let indel_data = if indel_observations.is_empty() {
            None
        } else {
            let counts =
                aggregate_indels(&indel_observations, total_reads, noisy_ref_count, pos);
            Some(Box::new(indels::IndelData {
                observations: indel_observations,
                ref_window: indel_ref_window,
                ref_anchor: indel_ref_anchor,
                homopolymer_run,
                dinucleotide_run,
                soft_clip_count,
                counts,
                calls: Vec::new(),
            }))
        };

        let mut before_counts = PairedCounts::default();
        let mut after_counts = PairedCounts::default();
        for read in reads.iter() {
            if read.strand == Strand::Unknown {
                continue;
            }
            if let Some(before) = read.before_base {
                before_counts.increment(ReadKey {
                    strand: read.strand,
                    current: read.base,
                    adj: before,
                });
            }
            if let Some(after) = read.after_base {
                after_counts.increment(ReadKey {
                    strand: read.strand,
                    current: read.base,
                    adj: after,
                });
            }
        }

        Ok(PileupMetrics {
            region,
            pos,
            reference_base,
            context,
            pos_metrics,
            pos_filters: Filters::default(),
            ref_metrics,
            alts,
            before_counts,
            after_counts,
            tags: RecordTags::default(),
            indel_data,
        })
    }

    pub fn ref_base(&self) -> Base {
        self.reference_base
    }

    pub fn contig(&self) -> SmolStr {
        self.region.contig.clone()
    }

    pub fn contig_name(&self) -> &str {
        &self.region.contig
    }

    pub fn pos(&self) -> u32 {
        self.pos
    }

    pub fn idx(&self) -> usize {
        self.region.pos_to_idx(self.pos).expect("valid position")
    }

    pub fn contig_pos(&self) -> SmolStr {
        use std::fmt::Write as _;
        let mut res = seqair_types::smol_str::SmolStrBuilder::new();

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
    pub fn new(
        total_reads: usize,
        reference_base: Base,
        before_1: Option<Base>,
        after_1: Option<Base>,
        baseq: RootMeanSquare,
        mapq: RootMeanSquare,
        mapq0: u32,
    ) -> Self {
        PositionMetrics {
            depth: u32::try_from(total_reads).expect("depth fits into u32"),
            baseq,
            mapq,
            mapq0,
            cpg: InCpG::new(reference_base, before_1, after_1),
            extended: PositionMetricsExt::default(),
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
    filters: SmallVec<RastairFilter, 6>,
}

impl Filters {
    pub fn add(&mut self, filter: RastairFilter, condition: impl FnOnce() -> bool) {
        if condition() && !self.filters.contains(&filter) {
            self.filters.push(filter);
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
    type Target = SmallVec<RastairFilter, 6>;

    fn deref(&self) -> &Self::Target {
        &self.filters
    }
}

#[derive(Debug, Clone, Default)]
pub(crate) struct AlleleAccumulator {
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
    pub(crate) fn add(&mut self, read: &SimpleRead, qual_sq: f64, mapq_sq: f64) {
        self.add_fields(
            qual_sq,
            mapq_sq,
            read.strand,
            read.matching_bases,
            read.indels,
            read.position.pos,
            read.position.read_length,
        );
    }

    pub(crate) fn add_fields(
        &mut self,
        qual_sq: f64,
        mapq_sq: f64,
        strand: Strand,
        matching_bases: u32,
        indels: u32,
        pos_in_read: u32,
        read_length: u32,
    ) {
        self.depth += 1;
        self.baseq.add_squared(qual_sq);
        self.mapq.add_squared(mapq_sq);
        match strand {
            Strand::OT => {
                self.ot_count += 1;
                self.baseq_ot.add_squared(qual_sq);
                self.mapq_ot.add_squared(mapq_sq);
            }
            Strand::OB => {
                self.ob_count += 1;
                self.baseq_ob.add_squared(qual_sq);
                self.mapq_ob.add_squared(mapq_sq);
            }
            Strand::Unknown => {}
        }
        self.aligned.add(f64::from(matching_bases));
        self.indels.add(f64::from(indels));
        self.pos_in_read.add(f64::from(pos_in_read) / f64::from(read_length));
    }

    pub(crate) fn finish(
        self,
        base: Base,
        total_reads: usize,
        pos: u32,
        ref_base: Base,
        context: &SequenceContext,
    ) -> Result<AlleleMetrics> {
        use Base::*;

        if self.depth == 0 {
            trace!(pos, ref_base = ?ref_base, ?base, pileup_reads = total_reads, "No reads for allele");
            return Ok(AlleleMetrics { base, ..default() });
        }

        if total_reads == 0 {
            bail!("allele has depth {} but pileup has 0 total reads — this is a bug", self.depth);
        }

        let denovo = if base == ref_base {
            FormsDenovo::No
        } else if context.before_1 == Some(C) && base == G {
            FormsDenovo::ThisBecomesG
        } else if context.after_1 == Some(G) && base == C {
            FormsDenovo::ThisBecomesC
        } else {
            FormsDenovo::No
        };

        Ok(AlleleMetrics {
            base,
            depth: self.depth,
            baseq: self.baseq.finish(),
            mapq: self.mapq.finish(),
            strand_count: ByStrand { ot: self.ot_count, ob: self.ob_count },
            baseq_s: ByStrand { ot: self.baseq_ot.finish(), ob: self.baseq_ob.finish() },
            mapq_s: ByStrand { ot: self.mapq_ot.finish(), ob: self.mapq_ob.finish() },
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
pub(crate) struct PerBaseAccumulators([AlleleAccumulator; 4]);

impl PerBaseAccumulators {
    pub(crate) fn accumulate(&mut self, read: &SimpleRead, qual_sq: f64, mapq_sq: f64) {
        let Some(idx) = read.base.known_index() else { return };
        self.0[idx].add(read, qual_sq, mapq_sq);
    }

    #[cfg(feature = "experimental-seqair")]
    pub(crate) fn accumulate_fields(
        &mut self,
        base: Base,
        qual_sq: f64,
        mapq_sq: f64,
        strand: Strand,
        matching_bases: u32,
        indels: u32,
        pos_in_read: u32,
        read_length: u32,
    ) {
        let Some(idx) = base.known_index() else { return };
        self.0[idx].add_fields(
            qual_sq,
            mapq_sq,
            strand,
            matching_bases,
            indels,
            pos_in_read,
            read_length,
        );
    }

    pub(crate) fn take(&mut self, base: Base) -> Option<AlleleAccumulator> {
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

pub struct MetricsForIndel<'p> {
    pub metrics: &'p PileupMetrics,
    pub indel: &'p crate::call::variant_calling::indel_calling::IndelCall,
}

pub(crate) fn aggregate_indels(
    indel_observations: &SmallVec<indels::IndelObservation, 0>,
    total_reads: usize,
    noisy_ref_count: u32,
    pos: u32,
) -> indels::IndelCounts {
    if indel_observations.is_empty() {
        return indels::IndelCounts {
            ref_count: total_reads as u32,
            noisy_ref_count,
            ..Default::default()
        };
    }

    let mut alleles: SmallVec<indels::IndelAlleleCounts, 2> = SmallVec::new();

    for obs in indel_observations {
        let entry = match alleles.iter_mut().find(|e| e.allele == obs.allele) {
            Some(entry) => entry,
            None => {
                alleles.push(indels::IndelAlleleCounts {
                    allele: obs.allele.clone(),
                    ot: 0,
                    ob: 0,
                    unknown_strand: 0,
                    noisy: 0,
                });
                alleles.last_mut().expect("just pushed")
            }
        };
        match obs.strand {
            Strand::OT => entry.ot += 1,
            Strand::OB => entry.ob += 1,
            Strand::Unknown => entry.unknown_strand += 1,
        }
        if obs.noisy {
            entry.noisy += 1;
        }
    }

    let total_indel_reads: u32 = alleles.iter().map(|a| a.total()).sum();
    let depth = total_reads as u32;
    // Both counts are drawn from one pass over the same alignments, so every
    // indel-carrying fragment is also part of the depth. If that stops holding,
    // `ref_count` floors to zero and every VAF here silently reads 1.0.
    if total_indel_reads > depth {
        warn!(
            pos,
            total_indel_reads,
            depth,
            "More indel-supporting fragments than reads at this position; the VAF \
             denominator is wrong. This is a bug in rastair, please report it."
        );
    }
    let ref_count = depth.saturating_sub(total_indel_reads);

    indels::IndelCounts { alleles, ref_count, noisy_ref_count }
}

#[cfg(test)]
mod size_tests {
    use super::*;

    // The budget covers the two inline `PairedCounts` tables (128 bytes each);
    // without them the struct fits in the original 800.
    #[test]
    fn pileup_metrics_size() {
        let size = std::mem::size_of::<PileupMetrics>();
        assert!(size < 1024, "PileupMetrics grew to {size} bytes");
    }
}
