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
use enumset::EnumSet;
use seqair_types::SmallVec;
use seqair_types::SmolStr;
use seqair_types::{Base, Probability, RmsAccumulator, RootMeanSquare, Strand};
use std::ops::Deref;
use std::sync::Arc;
use tracing::{trace, warn};

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PileupMetrics {
    /// The region this position came from, shared by every `PileupMetrics` in
    /// it. Inline it is 64 bytes — a `SmolStr` contig plus five `u64` — copied
    /// once per covered base and identical every time; behind an `Arc` it is 8.
    pub region: Arc<ChunkRegion>,
    pub pos: u32,
    pub reference_base: Base,
    pub context: SequenceContext,
    pub pos_metrics: PositionMetrics,
    pub pos_filters: Filters,
    pub ref_metrics: AlleleMetrics,
    /// Alternate alleles at this position.
    ///
    /// Inline capacity **one**, not two. An `Alt` is 152 bytes, so each inline
    /// slot is charged to every position in the genome, and almost no position
    /// uses the second one. Measured over 10,035,867 positions of chr12
    /// (NA12878, ~26x):
    ///
    /// | alts | positions | share |
    /// | ---: | ---: | ---: |
    /// | 0 | 9,573,352 | 95.39 % |
    /// | 1 | 454,755 | 4.53 % |
    /// | 2 | 7,564 | 0.08 % |
    /// | 3+ | 196 | 0.002 % |
    ///
    /// So the second inline slot cost 152 bytes at every position to save
    /// 7,760 heap allocations per 10 Mb — about one allocation per 1,300
    /// positions, against 1.5 GB of extra memory traffic over the same span.
    ///
    /// To re-measure after any change to alt calling, drop this into
    /// `get_pileups` in `src/call/process/pileups.rs`, just before
    /// `SlidingEntropy::new`:
    ///
    /// ```ignore
    /// let mut hist = [0u64; 8];
    /// for pm in &pileup_metrics {
    ///     hist[pm.alts.len().min(7)] += 1;
    /// }
    /// eprintln!("ALTSTAT {hist:?}");
    /// ```
    ///
    /// then sum the arrays over a run:
    ///
    /// ```text
    /// rastair call --gpu -f hg38.fa.gz in.bam -@ 8 -l chr12:20000000-30000000 --vcf /dev/null \
    ///   2>&1 | grep ALTSTAT | ...
    /// ```
    pub alts: SmallVec<Alt, 1>,
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
            let counts = aggregate_indels(&indel_observations, total_reads, noisy_ref_count, pos);
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

/// The set of FILTER codes a position or an alt allele has earned.
///
/// FILTER is a set, so this is a bitset: 13 variants fit in the `u16` that
/// `RastairFilter`'s `#[enumset(repr)]` pins down. Iteration is therefore in
/// discriminant order, which is also header registration order.
#[derive(Debug, Clone, Copy, Default, serde::Serialize, serde::Deserialize)]
pub struct Filters {
    pub other_pos_in_denovo_passes: bool,
    filters: EnumSet<RastairFilter>,
}

impl Filters {
    pub fn add(&mut self, filter: RastairFilter, condition: impl FnOnce() -> bool) {
        if condition() {
            self.filters.insert(filter);
        }
    }

    pub fn merge(&mut self, other: Filters) {
        self.filters |= other.filters;
    }

    pub fn pass(&self) -> bool {
        self.other_pos_in_denovo_passes || self.filters.is_empty()
    }

    pub fn is_empty(&self) -> bool {
        self.filters.is_empty()
    }

    /// The codes themselves, for a caller that has to union or emit them.
    /// Iterating an [`EnumSet`] yields discriminant order, which is the order
    /// the filters are registered in the VCF header.
    pub fn as_set(&self) -> EnumSet<RastairFilter> {
        self.filters
    }
}

/// Sum of squared values whose count is kept by the owner.
///
/// [`RmsAccumulator`] carries its own `count`, and [`AlleleAccumulator`] holds
/// nine of them whose counts are all `depth`, `ot_count` or `ob_count` — three
/// numbers it already maintains three lines away. It is updated once per read
/// per column, so those nine redundant increments were the hottest single line
/// in `call` (5.2 % of worker CPU), and the padding they carry made the struct
/// 160 bytes where 88 does — which matters again because
/// `PerBaseAccumulators::default()` re-zeroes four of them at every column.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct SumOfSquares(f64);

impl SumOfSquares {
    #[inline]
    fn add_squared(&mut self, x_sq: f64) {
        self.0 = self.0.algebraic_add(x_sq);
    }

    #[inline]
    fn add(&mut self, x: f64) {
        self.add_squared(x.algebraic_mul(x));
    }

    /// The RMS of the `count` values added.
    ///
    /// Deliberately routed through a one-element [`RmsAccumulator`] rather than
    /// taking the square root here: that accumulator divides by its own count of
    /// 1, which is exact, so `finish` evaluates the very same
    /// `sum.algebraic_div(count).sqrt()` the nine accumulators used to, bit for
    /// bit — including `RootMeanSquare(0.0)` for an empty one.
    fn finish(self, count: u32) -> RootMeanSquare {
        if count == 0 {
            return RootMeanSquare::default();
        }
        let mut acc = RmsAccumulator::new();
        acc.add_squared(self.0.algebraic_div(f64::from(count)));
        acc.finish()
    }
}

#[derive(Debug, Clone, Default)]
pub(crate) struct AlleleAccumulator {
    depth: u32,
    ot_count: u32,
    ob_count: u32,
    baseq: SumOfSquares,
    mapq: SumOfSquares,
    baseq_ot: SumOfSquares,
    baseq_ob: SumOfSquares,
    mapq_ot: SumOfSquares,
    mapq_ob: SumOfSquares,
    aligned: SumOfSquares,
    indels: SumOfSquares,
    pos_in_read: SumOfSquares,
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
            baseq: self.baseq.finish(self.depth),
            mapq: self.mapq.finish(self.depth),
            strand_count: ByStrand { ot: self.ot_count, ob: self.ob_count },
            baseq_s: ByStrand {
                ot: self.baseq_ot.finish(self.ot_count),
                ob: self.baseq_ob.finish(self.ob_count),
            },
            mapq_s: ByStrand {
                ot: self.mapq_ot.finish(self.ot_count),
                ob: self.mapq_ob.finish(self.ob_count),
            },
            num_aligned_bases: self.aligned.finish(self.depth),
            num_indels: self.indels.finish(self.depth),
            position_in_read: self.pos_in_read.finish(self.depth),
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
    indel_observations: &[indels::IndelObservation],
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

    /// A region holds one of these per covered base — 100,401 at the default
    /// `--segment-max-length` — and the pipeline walks that vec six or seven
    /// times, so this number is the memory traffic of the whole back half of
    /// `call`. It came down from 928 by sizing three fields for the common
    /// case rather than the tail (see `alts`, `PairedCounts`, and `region`),
    /// and from 568 by making `Filters` a bitset instead of a list; pinning it
    /// exactly means growing it again is a decision someone makes on purpose,
    /// with a measurement, rather than a field that slipped in.
    #[test]
    fn pileup_metrics_stays_small() {
        let size = std::mem::size_of::<PileupMetrics>();
        assert_eq!(
            size, 520,
            "PileupMetrics is {size} bytes. If that is deliberate, measure what it costs \
             (chr12:20–30 Mb, `--gpu -@ 8`, user CPU and peak RSS) and update this number."
        );
    }

    /// The nine per-allele sums gave up their own counts, so what says this was
    /// a refactor and not a numerical change is that `finish` returns *exactly*
    /// what `RmsAccumulator` returns for the same values — the same bits, not
    /// "close enough", because these reach the VCF as printed floats.
    #[test]
    fn sum_of_squares_is_bit_identical_to_rms_accumulator() {
        let cases: &[&[f64]] = &[
            &[],
            &[0.0],
            &[37.0],
            &[37.0, 41.0, 12.0],
            &[60.0; 1000],
            &[1e-8, 1e8, 3.5, 0.25],
            &[0.0, 0.0, 40.0],
        ];
        for values in cases {
            let mut reference = RmsAccumulator::new();
            let mut ours = SumOfSquares::default();
            for &v in *values {
                reference.add(v);
                ours.add(v);
            }
            let count = u32::try_from(values.len()).expect("test case fits in u32");
            assert_eq!(
                reference.finish().to_bits(),
                ours.finish(count).to_bits(),
                "diverged on {values:?}"
            );
        }
    }

    /// Four of these are built and zeroed at *every* column, and one is updated
    /// once per read per column, so this size is per-position memory traffic in
    /// the hottest loop `call` has. It was 160 when each of the nine sums
    /// carried a count that `depth`/`ot_count`/`ob_count` already held.
    #[test]
    fn allele_accumulator_stays_small() {
        assert_eq!(size_of::<AlleleAccumulator>(), 88);
        assert_eq!(size_of::<PerBaseAccumulators>(), 352);
    }

    // A `Filters` is copied into every `Alt`, so it is one of the few places
    // where a byte or two is worth pinning down.
    #[test]
    fn filters_are_a_bitset() {
        assert_eq!(size_of::<Filters>(), 4);
    }
}

#[cfg(test)]
mod filter_set_tests {
    use super::*;
    use crate::vcf::RastairFilter::{DnCpgBq, LowDp, LowMlScore, MVaf};

    #[test]
    fn a_filter_added_twice_is_present_once() {
        let mut filters = Filters::default();
        filters.add(LowDp, || true);
        filters.add(LowDp, || true);

        assert_eq!(filters.as_set().iter().collect::<Vec<_>>(), [LowDp]);
    }

    #[test]
    fn a_filter_is_only_added_when_its_condition_holds() {
        let mut filters = Filters::default();
        filters.add(LowDp, || false);

        assert!(filters.is_empty());
        assert!(filters.pass());
    }

    #[test]
    fn merging_unions_the_two_sets() {
        let mut filters = Filters::default();
        filters.add(MVaf, || true);
        filters.add(LowDp, || true);

        let mut other = Filters::default();
        other.add(DnCpgBq, || true);
        other.add(LowDp, || true);
        filters.merge(other);

        assert_eq!(filters.as_set().iter().collect::<Vec<_>>(), [LowDp, DnCpgBq, MVaf]);
    }

    // The FILTER column prints in this order, and `Schema::filter` indexes its
    // `FilterId` table by the same discriminant, so both follow the enum.
    #[test]
    fn iteration_follows_header_registration_order() {
        let mut filters = Filters::default();
        for filter in [LowMlScore, MVaf, DnCpgBq, LowDp] {
            filters.add(filter, || true);
        }

        assert_eq!(filters.as_set().iter().collect::<Vec<_>>(), [LowDp, DnCpgBq, MVaf, LowMlScore]);
    }

    // `other_pos_in_denovo_passes` overrides the set rather than living in it:
    // it is not a VCF FILTER code.
    #[test]
    fn the_denovo_override_passes_a_non_empty_set() {
        let mut filters = Filters::default();
        filters.add(LowDp, || true);
        assert!(!filters.pass());

        filters.other_pos_in_denovo_passes = true;
        assert!(filters.pass());
        assert!(!filters.is_empty());
    }
}
