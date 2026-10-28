use crate::{
    call::{
        phasing::RecordPhases,
        pileup::indels,
        variant_calling::{EstimatedGenotype, GtAllele},
    },
    metrics::PairedCounts,
    sequence::ChunkRegion,
    utils::{ByStrand, IntoF64, SequenceContext, default, logging::ThisIsABug},
    vcf::{CpgSide, Methylated, RastairFilter, TapsContext},
};
use better_default::Default;
use color_eyre::{
    Result,
    eyre::{Context, bail},
};
use enumset::EnumSet;
use seqair_types::SmallVec;
use seqair_types::SmolStr;
use seqair_types::{Base, Pos0, Probability, RootMeanSquare, Strand, SumOfSquares};
use std::{num::NonZeroU8, sync::Arc};
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
    pub pos: Pos0,
    pub reference_base: Base,
    pub context: SequenceContext,
    pub pos_metrics: PositionMetrics,
    /// The CpG sides a base here takes beside a called variant next door: an
    /// alt, or the base a called insertion puts beside it
    /// ([`CpgSide::beside`], set by `process::mark_called_cpgs`).
    pub beside_called_alts: EnumSet<CpgSide>,
    /// The CpG sides an alt here takes with a reference neighbour that no
    /// chromosome carrying the alt keeps: the neighbour's called alt sits on
    /// each (set by `process::mark_called_cpgs`).
    pub lost_to_called_alts: EnumSet<CpgSide>,
    pub ref_metrics: AlleleMetrics,
    /// Alternate alleles at this position.
    ///
    /// Inline capacity **one**, not two. An `Alt` is 120 bytes, so each inline
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
    /// So a second inline slot would cost 120 bytes at every position to save
    /// about one allocation per 1,300 positions.
    pub alts: SmallVec<Alt, 1>,
    /// Counts of (`my_base`, `before_base`) pairs by strand
    pub before_counts: PairedCounts,
    /// Counts of (`my_base`, `after_base`) pairs by strand
    pub after_counts: PairedCounts,
    pub indel_data: Option<Box<indels::IndelData>>,
    /// Each record's place in its phase block, set by `phase_segment` on
    /// the columns it phased.
    pub phase: Option<Box<RecordPhases>>,
    /// The pair-HMM's scores of the cluster holding this column's
    /// substitution events, if it has any.
    pub snv_scores: Option<Arc<crate::call::pair_hmm::HaplotypeScores>>,
}

#[derive(Debug, Clone)]
pub struct Alt {
    pub base: Base,
    pub metrics: AlleleMetrics,
    pub verdict: AltVerdict,
}

impl Alt {
    /// A real variant; otherwise methylation evidence or a read error.
    pub fn called(&self) -> bool {
        matches!(self.verdict, AltVerdict::Called(_))
    }
}

/// An alt's stage in calling: the forests score it, and `process::set_alt_calls`
/// alone calls or rejects it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AltVerdict {
    /// Not scored, or its score cleared for a rescoring.
    Unscored,
    /// The forests' probability that it is a true variant.
    Scored(Probability),
    /// A true variant, at the forests' probability.
    Called(Probability),
    Rejected(Rejection),
}

impl AltVerdict {
    /// The forests' probability that the alt is a true variant, where they
    /// scored it.
    pub fn ml(self) -> Option<Probability> {
        match self {
            Self::Scored(ml)
            | Self::Called(ml)
            | Self::Rejected(Rejection::LowScore(ml) | Rejection::NoGenotypeCopy(ml)) => Some(ml),
            Self::Rejected(Rejection::TooDeep(ml)) => ml,
            Self::Unscored | Self::Rejected(Rejection::PreMl | Rejection::NoScore) => None,
        }
    }
}

/// Why an alt is not called; each reason is one FILTER code.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Rejection {
    /// Too little evidence to be worth the forests.
    PreMl,
    /// A methylation candidate more reads show than `--m-max-coverage`:
    /// decided after [`Self::PreMl`] and before the score, which it keeps
    /// where the forests gave one.
    TooDeep(Option<Probability>),
    /// Worth scoring, but its features gave the forests no row.
    NoScore,
    /// Scored below the threshold.
    LowScore(Probability),
    /// Scored at the threshold, but the SNV engine's genotype carries no copy
    /// of it.
    NoGenotypeCopy(Probability),
}

impl Rejection {
    pub fn filter(self) -> RastairFilter {
        match self {
            Self::PreMl => RastairFilter::PreMl,
            Self::TooDeep(_) => RastairFilter::MHighDp,
            Self::NoScore | Self::LowScore(_) => RastairFilter::LowMlScore,
            Self::NoGenotypeCopy(_) => RastairFilter::NoGenotypeCopy,
        }
    }
}

impl PileupMetrics {
    pub fn ref_base(&self) -> Base {
        self.reference_base
    }

    pub fn contig_name(&self) -> &str {
        &self.region.contig
    }

    pub fn pos(&self) -> Pos0 {
        self.pos
    }

    pub fn contig_pos(&self) -> SmolStr {
        seqair_types::smol_str::format_smolstr!("{}:{}", self.contig_name(), self.pos())
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

    /// The `n`-th alt, as genotypes number it.
    pub fn alt_at(&self, n: NonZeroU8) -> Option<&Alt> {
        self.alts.get(usize::from(n.get() - 1))
    }

    /// The base a genotype allele reads as.
    pub fn allele_base(&self, allele: GtAllele) -> Option<Base> {
        match allele {
            GtAllele::Ref => Some(self.reference_base),
            GtAllele::Alt(n) => self.alt_at(n).map(|alt| alt.base),
        }
    }

    pub fn find_alt(&self, alt: Base) -> Option<&Alt> {
        self.alts.iter().find(|a| a.base == alt)
    }

    pub fn find_alt_mut(&mut self, alt: Base) -> Option<&mut Alt> {
        self.alts.iter_mut().find(|a| a.base == alt)
    }

    /// Whether the SNV engine scored `alt`'s substitution at this column.
    pub fn engine_scored(&self, alt: Base) -> bool {
        self.snv_scores.as_deref().is_some_and(|s| s.base_event_at(self.pos, alt).is_some())
    }

    /// The side of a reference CpG this column is, if any.
    pub fn reference_cpg(&self) -> Option<CpgSide> {
        TapsContext::of(self).in_cpg()
    }

    /// A read reaches this column.
    pub fn covered(&self) -> bool {
        self.pos_metrics.depth > 0
    }

    /// The CpG side `base` takes here: beside its reference neighbours, or
    /// beside a called variant next door.
    pub fn cpg_side_of(&self, base: Base) -> Option<CpgSide> {
        let reference = TapsContext::of(self).sides() - self.lost_to_called_alts;
        CpgSide::taken_by(base, reference | self.beside_called_alts)
    }

    /// The side this reference base takes in a CpG a called variant next
    /// door makes, and the reference does not: a reference CpG's `C` stays
    /// one beside an inserted `G`.
    pub fn denovo_partner(&self) -> Option<CpgSide> {
        let side = CpgSide::taken_by(self.reference_base, self.beside_called_alts);
        side.filter(|_| self.reference_cpg().is_none())
    }

    /// A called alt that makes no CpG, or a called indel allele.
    pub fn is_variant(&self) -> bool {
        self.alts.iter().any(|a| a.called() && self.cpg_side_of(a.base).is_none())
            || self.called_indels().next().is_some()
    }

    /// The indel alleles anchored here whose call carries them.
    pub fn called_indels(&self) -> impl Iterator<Item = &indels::IndelAllele> {
        let scored = self.indel_data.as_ref().and_then(|d| d.scored.as_ref());
        let calls = scored.into_iter().flat_map(|s| s.records.scored_calls());
        calls.filter(|(call, _)| call.genotype.carries_an_alt()).map(|(call, _)| &call.allele)
    }

    /// In a CpG, the reference's or one a called alt makes, here or next
    /// door: the columns `M5mC` is written for.
    pub fn in_cpg(&self) -> bool {
        !super::methylation::origins(self).is_empty()
    }

    /// In a CpG a called alt makes, here or next door.
    pub fn in_denovo_cpg(&self) -> bool {
        self.denovo_partner().is_some()
            || super::methylation::origins(self).contains(&crate::vcf::CpgOrigin::DeNovo)
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
    /// Entropy of the surrounding region
    pub region_entropy: f64,
    /// Estimated genotype, set once the alts are called
    pub genotype: Option<EstimatedGenotype>,
    /// Methylation beta, set once the alts are called
    pub methylated: Methylated,
}

impl PositionMetrics {
    pub fn new(
        depth: u32,
        baseq: RootMeanSquare,
        mapq: RootMeanSquare,
        mapq0: u32,
        region_entropy: f64,
    ) -> Self {
        PositionMetrics {
            depth,
            baseq,
            mapq,
            mapq0,
            region_entropy,
            genotype: None,
            methylated: Methylated::default(),
        }
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
    pub strand: Option<Strand>,
    pub matching_bases: u32,
    pub indels: u32,
    pub pos_in_read: u32,
    pub read_length: u32,
}

/// The `[OT, OB]` slot a strand accumulates into.
const fn strand_slot(strand: Strand) -> usize {
    match strand {
        Strand::OT => 0,
        Strand::OB => 1,
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
        if let Some(slot) = strand.map(strand_slot)
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
        pos: Pos0,
        ref_base: Base,
    ) -> Result<AlleleMetrics> {
        if self.depth == 0 {
            trace!(%pos, ref_base = ?ref_base, ?base, pileup_reads = total_reads, "No reads for allele");
            return Ok(AlleleMetrics { base, ..default() });
        }

        if total_reads == 0 {
            bail!("allele has depth {} but pileup has 0 total reads — this is a bug", self.depth);
        }

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
    /// The side of a de-novo CpG the alt makes with a reference neighbour:
    /// what the forests see, before any alt is called.
    pub fn denovo_side(&self) -> Option<CpgSide> {
        TapsContext::of(self.metrics).forms_cpg(self.alt.base)
    }

    /// Is the alt what a methylated CpG reference base reads as? Only where
    /// the reference base takes the CpG side, in a reference CpG, where the
    /// CpG forest is trained, or beside a called indel's base (alts are
    /// called after the forests score), and on either strand: it picks the
    /// forest for the alt, not a weight for its reads.
    pub fn is_conversion(&self) -> bool {
        let metrics = self.metrics;
        let conversion = CpgSide::converting(metrics.reference_base, self.alt.base);
        conversion.is_some_and(|side| metrics.cpg_side_of(metrics.reference_base) == Some(side))
    }
}

#[cfg(test)]
mod size_tests {
    use super::*;
    use seqair_types::RmsAccumulator;

    /// A region holds one of these per covered base — 100,401 at the default
    /// `--segment-max-length` — and the pipeline walks that vec six or seven
    /// times, so this number is the memory traffic of the whole back half of
    /// `call`. Pinned exactly, so growing it is a decision made on purpose,
    /// with a measurement.
    #[test]
    fn pileup_metrics_stays_small() {
        let size = std::mem::size_of::<PileupMetrics>();
        assert_eq!(
            size, 616,
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
    /// the hottest loop `call` has.
    #[test]
    fn allele_accumulator_stays_small() {
        assert_eq!(size_of::<AlleleAccumulator>(), 88);
        assert_eq!(size_of::<PerBaseAccumulators>(), 352);
    }

    #[test]
    fn alt_stays_small() {
        assert_eq!(size_of::<AltVerdict>(), 16);
        assert_eq!(size_of::<Alt>(), 120, "the `alts` doc reasons with this size");
    }
}
