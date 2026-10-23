use crate::{
    call::{phasing::RecordPhases, pileup::indels, variant_calling::EstimatedGenotype},
    metrics::{MethylationEvidenceStrandInfo, PairedCounts},
    sequence::ChunkRegion,
    utils::{ByStrand, IntoF64, SequenceContext, default, logging::ThisIsABug},
    vcf::{InCpG, Methylated, RastairFilter, TapsContext},
};
use better_default::Default;
use color_eyre::{
    Result,
    eyre::{Context, bail},
};
use enumset::EnumSet;
use seqair_types::SmallVec;
use seqair_types::SmolStr;
use seqair_types::{Base, Probability, RootMeanSquare, Strand, SumOfSquares};
use std::ops::Deref;
use std::sync::Arc;
use tracing::trace;

/// One covered base's metrics.
///
/// Pinned in size (`pileup_metrics_stays_small`): a region holds one per
/// covered base and the pipeline walks that vec many times. Per-feature data a
/// run may not use lives behind one `Option<Box<_>>` (`phase`, `indel_data`),
/// so an unphased run pays one pointer.
#[derive(Debug, Clone)]
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
    pub indel_data: Option<Box<indels::IndelData>>,
    /// Each record's place in its phase block, set by `phase_segment` on
    /// the columns it phased.
    pub phase: Option<Box<RecordPhases>>,
    /// The pair-HMM's scores of the cluster holding this column's
    /// substitution events, if it has any.
    pub snv_scores: Option<Arc<crate::call::pair_hmm::HaplotypeScores>>,
}

#[derive(Debug, Clone, Default)]
pub struct RecordTags {
    // The tags have been calculated. Mainly to debug :)
    pub set: bool,
    /// This position has coverage
    pub covered: bool,
    /// This is a de-novo CpG site (not the partner)
    pub denovo_cpg: bool,
    /// This is the partner position of a de-novo CpG site
    pub denovo_cpg_partner: bool,
    /// This position is a variant (but not a denovo CpG)
    pub variant: bool,
}

#[derive(Debug, Clone)]
pub struct Alt {
    pub base: Base,
    pub metrics: AlleleMetrics,
    pub filters: AltFilters,
    pub call: AltCall,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
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

    /// Whether the SNV engine scored `alt`'s substitution at this column.
    pub fn engine_scored(&self, alt: Base) -> bool {
        self.snv_scores.as_deref().is_some_and(|s| s.base_event_at(self.pos, alt).is_some())
    }

    /// Get all alternative bases in the pileup (for lookup with mutation)
    pub fn alts(&self) -> SmallVec<Base, 4> {
        self.alts.iter().map(|a| a.base).collect()
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

/// Alternate alleles are listed by support, deepest first; equal support
/// keeps the order they were observed in.
///
/// Observation order is the pileup engine's read order at equal start
/// positions, an accident of the engine, while the ALT column, `AD`, `ML` and
/// the genotype indices all follow this order. Ties are left alone on purpose:
/// genotyping and methylation break their own ties by this order.
pub fn order_alts(alts: &mut [Alt]) {
    alts.sort_by_key(|alt| std::cmp::Reverse(alt.metrics.depth));
}

#[derive(Debug, Clone)]
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
    pub extended: PositionMetricsExt,
}

#[derive(Debug, Clone, Default)]
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

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
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
            cpg: TapsContext::new(reference_base, before_1, after_1).in_cpg(),
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

#[derive(Debug, Clone, Default)]
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

/// Whether an alt forms a CpG; [`TapsContext::forms_cpg`] decides it from the
/// reference, `call::process::mark_variant_pair_cpgs` from a called neighbour.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
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

#[derive(Debug, Clone, Default)]
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
/// FILTER is a set, so this is a bitset: 14 variants fit in the `u16` that
/// `RastairFilter`'s `#[enumset(repr)]` pins down. Iteration is therefore in
/// discriminant order, which is also header registration order.
#[derive(Debug, Clone, Copy, Default)]
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

    /// Insert `filter` when `on`, remove it otherwise.
    pub fn set(&mut self, filter: RastairFilter, on: bool) {
        if on {
            self.filters.insert(filter);
        } else {
            self.filters.remove(filter);
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

#[derive(Debug, Clone, Default)]
pub(crate) struct AlleleAccumulator {
    depth: u32,
    /// `[OT, OB]`; kept out of `by_strand` so neither array needs padding.
    strand_depth: [u32; 2],
    // Field order is load-bearing. rustc keeps declaration order among equally
    // aligned fields, and LLVM pairs two sums into one 128-bit `fadd.2d` only
    // when they are adjacent *and* reached at a constant offset. `baseq`/`mapq`
    // and `aligned`/`indels` satisfy both and vectorise; the strand-split pair
    // is `[[f64; 2]; 2]` rather than four named fields for the same reason —
    // see `add_fields`.
    baseq: SumOfSquares,
    mapq: SumOfSquares,
    /// `[OT, OB]`, each `[baseq, mapq]`.
    by_strand: [[f64; 2]; 2],
    aligned: SumOfSquares,
    indels: SumOfSquares,
    pos_in_read: SumOfSquares,
}

/// One read's contribution to its allele's accumulator.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ReadFields {
    pub qual_sq: f64,
    pub mapq_sq: f64,
    pub strand: Strand,
    pub matching_bases: u32,
    pub indels: u32,
    pub pos_in_read: u32,
    pub read_length: u32,
}

/// The `[OT, OB]` slot a known strand accumulates into.
const fn strand_slot(strand: Strand) -> Option<usize> {
    match strand {
        Strand::OT => Some(0),
        Strand::OB => Some(1),
        Strand::Unknown => None,
    }
}

impl AlleleAccumulator {
    pub(crate) fn add_fields(&mut self, read: ReadFields) {
        let ReadFields {
            qual_sq,
            mapq_sq,
            strand,
            matching_bases,
            indels,
            pos_in_read,
            read_length,
        } = read;
        self.depth += 1;
        self.baseq.add_squared(qual_sq);
        self.mapq.add_squared(mapq_sq);
        // Indexed, not matched. Whether consecutive reads are OT or OB is close
        // to a coin flip, so the `match` this replaces mispredicted once per
        // read per column; and LLVM merged its two arms into a common tail
        // addressed by a register, which is what stopped it pairing the two
        // sums into one `fadd.2d`. One in-range slot fixes both.
        if let Some(slot) = strand_slot(strand)
            && let Some(strand_depth) = self.strand_depth.get_mut(slot)
            && let Some([baseq, mapq]) = self.by_strand.get_mut(slot)
        {
            *strand_depth += 1;
            *baseq = baseq.algebraic_add(qual_sq);
            *mapq = mapq.algebraic_add(mapq_sq);
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
        if self.depth == 0 {
            trace!(pos, ref_base = ?ref_base, ?base, pileup_reads = total_reads, "No reads for allele");
            return Ok(AlleleMetrics { base, ..default() });
        }

        if total_reads == 0 {
            bail!("allele has depth {} but pileup has 0 total reads — this is a bug", self.depth);
        }

        let denovo = TapsContext::new(ref_base, context.before_1, context.after_1).forms_cpg(base);

        let [ot_depth, ob_depth] = self.strand_depth;
        let [[ot_baseq, ot_mapq], [ob_baseq, ob_mapq]] = self.by_strand;

        Ok(AlleleMetrics {
            base,
            depth: self.depth,
            baseq: self.baseq.finish(self.depth),
            mapq: self.mapq.finish(self.depth),
            strand_count: ByStrand { ot: ot_depth, ob: ob_depth },
            baseq_s: ByStrand {
                ot: SumOfSquares::from(ot_baseq).finish(ot_depth),
                ob: SumOfSquares::from(ob_baseq).finish(ob_depth),
            },
            mapq_s: ByStrand {
                ot: SumOfSquares::from(ot_mapq).finish(ot_depth),
                ob: SumOfSquares::from(ob_mapq).finish(ob_depth),
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
    pub(crate) fn accumulate_fields(&mut self, base: Base, read: ReadFields) {
        if let Some(allele) = base.known_index().and_then(|idx| self.0.get_mut(idx)) {
            allele.add_fields(read);
        }
    }

    pub(crate) fn take(&mut self, base: Base) -> Option<AlleleAccumulator> {
        self.0.get_mut(base.known_index()?).map(std::mem::take)
    }
}

pub struct MetricsForAlt<'p> {
    pub metrics: &'p PileupMetrics,
    pub alt: &'p AlleleMetrics,
}

impl MetricsForAlt<'_> {
    /// See [`TapsContext::is_conversion`].
    pub fn is_conversion(&self) -> bool {
        self.metrics.pos_metrics.cpg.mod_base() == Some(self.alt.base)
    }
}

#[cfg(test)]
mod size_tests {
    use super::*;
    use seqair_types::RmsAccumulator;

    /// A region holds one of these per covered base — 100,401 at the default
    /// `--segment-max-length` — and the pipeline walks that vec six or seven
    /// times, so this number is the memory traffic of the whole back half of
    /// `call`. It came down from 928 by sizing three fields for the common
    /// case rather than the tail (see `alts`, `PairedCounts`, and `region`),
    /// and from 568 by making `Filters` a bitset instead of a list; pinning it
    /// exactly means growing it again is a decision someone makes on purpose,
    /// with a measurement, rather than a field that slipped in. Phasing adds
    /// one pointer, behind which all of its per-column data lives.
    #[test]
    fn pileup_metrics_stays_small() {
        let size = std::mem::size_of::<PileupMetrics>();
        assert_eq!(
            size, 528,
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
