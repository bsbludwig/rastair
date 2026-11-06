//! One column's reads, accumulated into a [`ColumnDraft`] and finished into
//! [`PileupMetrics`].
//!
//! [`ColumnDraft::accumulate`] counts a column's reads; [`ColumnDraft::finish`]
//! turns a draft with [`ColumnEvidence`] into metrics. Nothing skips a column
//! with zero reads or an `N` reference (`Base::Unknown`) before it is
//! finished, so both are real paths: return default metrics, not an error.
//!
//! **Driving the engine.** `PileupEngine::new` takes a `PileupInput`, and only
//! `store.prepare_for_pileup()` makes one, sorted by position with mates
//! linked; without either, columns go missing or overlap dedup silently does
//! nothing. `PileupColumn::position_of` and `alignment_at` find a view's
//! linked mate in the column, all `counted_mate` needs; pairwise rules stay on
//! this side of the boundary.
//! `PileupEngine::reclaim_allocation` returns an *empty* store that keeps its
//! capacity for the next region.
//!
//! **Indel observations.** Each rule guards a class of mistake:
//!
//! - A deletion's REF comes from the deletion site: `view.qpos()` is a read
//!   offset, not the column's genomic `pos` (seqair's `QPos` keeps them apart).
//! - Overlap dedup keeps one read per fragment and can keep the mate that does
//!   not span the indel, so the dropped mate stands in:
//!   `build_indel_observation(view).or_else(mate)`.
//! - The CIGAR counts only pick candidates: `IndelCounts` holds, per allele,
//!   the fragments whose alignment carries it, and nothing reads it for a
//!   call. There are no read-level indel filters; the pair-HMM weighs each
//!   read's alignment itself.

use crate::utils::MappingQualityExt;
use crate::{
    call::{
        ColumnEvidence,
        pair_hmm::HaplotypeScratch,
        pileup::indels::{IndelAllele, IndelCounts, IndelData},
        process::{ColumnKeep, Engine, PileupMappingParams},
    },
    metrics::{
        AlleleMetrics, Alt, AltVerdict, PairedCounts, PerBaseAccumulators, PileupMetrics,
        PositionMetrics, ReadFields, ReadKey, order_alts,
    },
    sequence::{RastairReadExtras, Segment},
    utils::SequenceContext,
    vcf::{CpgSide, TapsContext},
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use enumset::EnumSet;
use seqair::bam::{
    RecordIdx,
    pileup::{AlignmentView, Indel, PileupColumn, PileupOp},
};
use seqair_types::{Base, Offset, Pos0, QPos, SmallVec, Strand, SumOfSquares};

/// The distance from a column to its neighbour.
const ONE_BASE: Offset = Offset::new(1);
use std::rc::Rc;
use tracing::{instrument, warn};

/// One column's reads, accumulated but not yet reduced to a [`PileupMetrics`].
///
/// The split exists so a column can be rejected before the expensive half runs:
/// [`ColumnDraft::finish`] pays nine divisions and square roots per allele, a
/// `PositionMetrics`, and a ~530-byte struct write, and roughly seven columns in
/// eight are then dropped for want of [`ColumnEvidence`]. Everything that
/// filter reads is already known here — see [`ColumnDraft::evidence`].
pub(crate) struct ColumnDraft {
    segment: Rc<Segment>,
    pos: Pos0,
    idx: usize,
    reference_base: Base,
    context: SequenceContext,
    accumulators: PerBaseAccumulators,
    /// Over `total_depth` reads, which is their count.
    pos_baseq: SumOfSquares,
    pos_mapq: SumOfSquares,
    mapq0: u32,
    total_depth: usize,
    alt_bases: SmallVec<Base, 4>,
    indel_observations: SmallVec<IndelAllele, 3>,
    before_counts: PairedCounts,
    after_counts: PairedCounts,
}

/// What a column contributes to its neighbours' de-novo adjacency
/// ([`CpgSide::beside`]): the CpG sides of its alts, and of the last base
/// its reads insert after it, the next column's neighbour on their
/// chromosome (it may keep a `C` before it needlessly).
#[derive(Debug, Clone, Copy)]
pub(crate) struct DenovoNeighbour {
    pos: Pos0,
    alt_sides: EnumSet<CpgSide>,
}

impl DenovoNeighbour {
    /// Its alts' sides if it sits at `pos`, none otherwise.
    fn sides_at(neighbour: Option<Self>, pos: Option<Pos0>) -> EnumSet<CpgSide> {
        neighbour.filter(|n| Some(n.pos) == pos).map(|n| n.alt_sides).unwrap_or_default()
    }
}

/// The de-novo neighbour of a column that [`ColumnDraft::accumulate`] would
/// reject without a trace, or `None` when the column needs the full walk.
///
/// That is a column where every read shows the reference base with nothing
/// inserted or deleted after it, so it has no alt and no indel, and whose
/// reference is an `A`, a `T`, or a `G` that neither the reference (a `C`
/// before it) nor `previous` (an adjacent alt `C`) makes a CpG. A `C` always
/// takes the walk: the column after it may still rescue it. Every read, not only the ones the
/// filters keep, so this never claims a column the walk would keep.
pub(crate) fn plain_reference_column(
    column: &PileupColumn<'_, RastairReadExtras>,
    segment: &Segment,
    params: &PileupMappingParams,
    previous: Option<DenovoNeighbour>,
) -> Option<DenovoNeighbour> {
    // Every column is kept, so none is plain.
    if params.columns == ColumnKeep::Every {
        return None;
    }
    let pos = column.pos();
    let idx = segment.pos_to_idx(pos).ok()?;
    let reference = Base::from(*segment.sequence.get(idx)?);
    let settled = match reference {
        Base::A | Base::T => true,
        Base::G => {
            let after_c = idx
                .checked_sub(1)
                .and_then(|i| segment.sequence.get(i))
                .is_some_and(|&b| Base::from(b) == Base::C);
            let before = DenovoNeighbour::sides_at(previous, pos.checked_sub_offset(ONE_BASE));
            let partners = CpgSide::beside(before, EnumSet::empty());
            !after_c && CpgSide::taken_by(reference, partners).is_none()
        }
        _ => false,
    };
    if !settled {
        return None;
    }
    let plain = column.raw_alignments().all(|aln| {
        matches!(aln.op, PileupOp::Match { base, .. } if base == reference)
            && matches!(aln.indel_after(), Indel::None)
    });
    plain.then_some(DenovoNeighbour { pos, alt_sides: EnumSet::empty() })
}

impl ColumnDraft {
    /// `scratch` is reset here; it lives across columns only so a deep pileup
    /// does not allocate once per position.
    #[instrument(level = "trace", skip_all)]
    pub(crate) fn accumulate(
        column: &PileupColumn<'_, RastairReadExtras>,
        segment: Rc<Segment>,
        params: &PileupMappingParams,
        scratch: &mut ColumnScratch,
    ) -> Result<ColumnDraft> {
        let pos = column.pos();
        let idx = segment.pos_to_idx(pos)?;
        let depth = column.depth();

        let reference_base: Base =
            segment.sequence.get(idx).wrap_err("failed to get reference base")?.into();

        let context =
            SequenceContext::new(idx, &segment).wrap_err("failed to get sequence context")?;

        let dedup_overlaps = !params.keep_overlapping_reads;
        scratch.begin(depth);

        let mut draft = Self::empty(segment, pos, idx, reference_base, context);

        for (slot, view) in column.alignments().enumerate() {
            // Already resolved if the overlap rule reached forward for this
            // read while deciding its mate; otherwise this is its first ask.
            let resolved = match scratch.mates.resolved(slot) {
                Some(known) => known.observed,
                None => observed(&view, params, reference_base, &draft.context),
            };
            let Some(Observed { base, baseq, qpos }) = resolved else {
                continue;
            };
            let strand = view.extra().strand;
            let overlap = if dedup_overlaps {
                drops_overlapping_mate(
                    column,
                    &view,
                    base,
                    params,
                    reference_base,
                    &draft.context,
                    scratch,
                )
            } else {
                Overlap::Keep
            };
            if overlap == Overlap::Drop {
                continue;
            }
            draft.total_depth += 1;

            let qual_sq = f64::from(baseq).powi(2);
            let mapq_sq = f64::from(view.mapq.or_unique()).powi(2);
            draft.accumulators.accumulate_fields(
                base,
                ReadFields {
                    qual_sq,
                    mapq_sq,
                    strand,
                    matching_bases: view.matching_bases,
                    indels: view.indel_bases,
                    pos_in_read: qpos.get(),
                    read_length: view.seq_len,
                },
            );
            draft.pos_baseq.add_squared(qual_sq);
            draft.pos_mapq.add_squared(mapq_sq);
            if view.mapq.get() == Some(0) {
                draft.mapq0 += 1;
            }
            if base.known_index().is_some()
                && base != reference_base
                && !draft.alt_bases.contains(&base)
            {
                draft.alt_bases.push(base);
            }

            if let Some(strand) = strand {
                let seq = view.seq();
                if let Some(&adj) = qpos.checked_sub(1).and_then(|i| seq.get(i.as_usize())) {
                    draft.before_counts.increment(ReadKey { strand, current: base, adj });
                }
                // The read's neighbour, across an indel too: an inserted `G`
                // after a `C` makes a CpG of it.
                if let Some(&adj) = qpos.checked_add(1).and_then(|i| seq.get(i.as_usize())) {
                    draft.after_counts.increment(ReadKey { strand, current: base, adj });
                }
            }

            if params.engine == Engine::PairHmm {
                let mate = dedup_overlaps
                    .then(|| {
                        counted_mate(
                            column,
                            &view,
                            params,
                            reference_base,
                            &draft.context,
                            &mut scratch.mates,
                        )
                    })
                    .flatten();
                // A deduplicated column keeps one read per fragment, and that
                // read is not necessarily the one carrying the fragment's
                // indel: when the mates disagree the rule can drop exactly the
                // mate that spans it, and the fragment's only indel observation
                // would go with it. The mate then speaks for the fragment.
                let evidence =
                    build_indel_observation(&view, idx, draft.segment.as_ref()).or_else(|| {
                        build_indel_observation(mate.as_ref()?, idx, draft.segment.as_ref())
                    });
                if let Some(obs) = evidence {
                    draft.indel_observations.push(obs);
                }
            }
        }

        Ok(draft)
    }

    /// The draft of a position no read reaches, or `None` unless the
    /// reference has a CpG there: the only such positions any output reports.
    pub(crate) fn uncovered_cpg(pos: Pos0, segment: &Rc<Segment>) -> Result<Option<ColumnDraft>> {
        let idx = segment.pos_to_idx(pos)?;
        let base_at =
            |idx: Option<usize>| idx.and_then(|i| segment.sequence.get(i)).map(|&b| Base::from(b));
        let reference_base = base_at(Some(idx)).wrap_err("failed to get reference base")?;
        // The CpG test is all most gap positions need, so it comes before the
        // full context.
        let cpg = TapsContext::new(
            reference_base,
            base_at(idx.checked_sub(1)),
            base_at(idx.checked_add(1)),
        )
        .in_cpg();
        if cpg.is_none() {
            return Ok(None);
        }
        let context =
            SequenceContext::new(idx, segment).wrap_err("failed to get sequence context")?;
        Ok(Some(Self::empty(segment.clone(), pos, idx, reference_base, context)))
    }

    fn empty(
        segment: Rc<Segment>,
        pos: Pos0,
        idx: usize,
        reference_base: Base,
        context: SequenceContext,
    ) -> ColumnDraft {
        ColumnDraft {
            segment,
            pos,
            idx,
            reference_base,
            context,
            accumulators: PerBaseAccumulators::default(),
            pos_baseq: SumOfSquares::new(),
            pos_mapq: SumOfSquares::new(),
            mapq0: 0,
            total_depth: 0,
            alt_bases: SmallVec::new(),
            indel_observations: SmallVec::new(),
            before_counts: PairedCounts::default(),
            after_counts: PairedCounts::default(),
        }
    }

    fn taps_context(&self) -> TapsContext {
        TapsContext::new(self.reference_base, self.context.before_1, self.context.after_1)
    }

    /// The sequence index of this column, for the caller's sliding entropy.
    pub(crate) fn idx(&self) -> usize {
        self.idx
    }

    pub(crate) fn denovo_neighbour(&self) -> DenovoNeighbour {
        let inserted = self.indel_observations.iter().filter_map(IndelAllele::last_inserted);
        let bases = self.alt_bases.iter().copied().chain(inserted);
        DenovoNeighbour { pos: self.pos, alt_sides: bases.filter_map(CpgSide::of_base).collect() }
    }

    /// Only a reference `C` can still be rescued by the column after it, so
    /// every other draft's fate is settled the moment it is accumulated.
    pub(crate) fn awaits_successor(&self) -> bool {
        self.reference_base == Base::C
    }

    /// The draft's [`ColumnEvidence`].
    ///
    /// `before`/`after` are the neighbouring columns the engine yielded, kept
    /// or not; each contributes only when it is genomically adjacent.
    pub(crate) fn evidence(
        &self,
        before: Option<DenovoNeighbour>,
        after: Option<DenovoNeighbour>,
    ) -> EnumSet<ColumnEvidence> {
        let before = DenovoNeighbour::sides_at(before, self.pos.checked_sub_offset(ONE_BASE));
        let after = DenovoNeighbour::sides_at(after, self.pos.checked_add_offset(ONE_BASE));
        let mut evidence = EnumSet::empty();
        if !self.alt_bases.is_empty() {
            evidence |= ColumnEvidence::Alt;
        }
        if CpgSide::taken_by(
            self.reference_base,
            self.taps_context().sides() | CpgSide::beside(before, after),
        )
        .is_some()
        {
            evidence |= ColumnEvidence::Cpg;
        }
        if !self.indel_observations.is_empty() {
            evidence |= ColumnEvidence::Indel;
        }
        evidence
    }

    /// Reduce the accumulated column to the metrics the rest of the pipeline
    /// consumes. The entropy of the reference around it is the caller's: it
    /// slides one window over the columns it keeps.
    pub(crate) fn finish(self, region_entropy: f64) -> Result<PileupMetrics> {
        let ColumnDraft {
            segment,
            pos,
            idx: _,
            reference_base,
            context,
            mut accumulators,
            pos_baseq,
            pos_mapq,
            mapq0,
            total_depth,
            alt_bases,
            indel_observations,
            before_counts,
            after_counts,
        } = self;

        let depth = u32::try_from(total_depth).wrap_err("column depth exceeds u32")?;
        let pos_metrics = PositionMetrics::new(
            depth,
            pos_baseq.finish(depth),
            pos_mapq.finish(depth),
            mapq0,
            region_entropy,
        );

        let ref_metrics = if let Some(acc) = accumulators.take(reference_base) {
            acc.finish(reference_base, total_depth, pos, reference_base)?
        } else {
            AlleleMetrics { base: reference_base, ..Default::default() }
        };

        let mut alts: SmallVec<Alt, 1> = alt_bases
            .iter()
            .map(|&base| {
                let acc = accumulators
                    .take(base)
                    .ok_or_else(|| color_eyre::eyre::eyre!("unknown base {base} in alt_bases"))?;
                let metrics = acc.finish(base, total_depth, pos, reference_base)?;
                Ok(Alt { base, metrics, verdict: AltVerdict::Unscored })
            })
            .collect::<Result<_>>()?;
        order_alts(&mut alts);

        let indel_data = if indel_observations.is_empty() {
            None
        } else {
            let counts = IndelCounts::aggregate(&indel_observations);
            Some(Box::new(IndelData { counts, scored: None }))
        };

        Ok(PileupMetrics {
            region: segment.range.clone(),
            pos,
            reference_base,
            context,
            pos_metrics,
            beside_called_alts: EnumSet::empty(),
            lost_to_called_alts: EnumSet::empty(),
            unconfounded: EnumSet::empty(),
            ref_metrics,
            alts,
            before_counts,
            after_counts,
            indel_data,
            phase: None,
            snv_scores: None,
        })
    }
}

/// What the overlap rule makes of a read in this column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Overlap {
    Drop,
    Keep,
}

/// Decide the fate of a read that shares a mate overlap with another read in
/// this column.
///
/// Called once per read inside an overlap, in column order. The left mate
/// (lower record index) is always seen first, so it is the one that decides:
/// it looks its mate up in the column, applies the overlap rule — drop the
/// second-in-template read, or the later one when both show the same base —
/// and, when it wins, records the mate for the `mate_drops` check the right
/// mate then hits. Deciding at the left mate keeps every kept read
/// accumulating in column order.
fn drops_overlapping_mate(
    column: &PileupColumn<'_, RastairReadExtras>,
    view: &AlignmentView<'_, '_, RastairReadExtras>,
    base: Base,
    params: &PileupMappingParams,
    reference_base: Base,
    context: &SequenceContext,
    scratch: &mut ColumnScratch,
) -> Overlap {
    if !view.in_mate_overlap() {
        return Overlap::Keep;
    }
    let this_idx = view.alignment().record_idx();
    let Some(mate_idx) = view.alignment().mate_idx() else { return Overlap::Keep };

    if this_idx > mate_idx {
        return if scratch.mate_drops.binary_search(&this_idx).is_ok() {
            Overlap::Drop
        } else {
            Overlap::Keep
        };
    }

    // A mate that is absent from this column, or that fails a filter here,
    // never formed a pair — this read stands on its own.
    let Some(slot) = column.position_of(mate_idx) else { return Overlap::Keep };
    let Some(mate) = scratch.mates.resolve(column, slot, params, reference_base, context) else {
        return Overlap::Keep;
    };
    let Some(Observed { base: mate_base, .. }) = mate.observed else {
        return Overlap::Keep;
    };

    // A pair counts once: by read 1's base when the bases disagree, by the
    // earlier read when they agree.
    let disagree = base != mate_base;
    if disagree && !mate.second_in_template {
        return Overlap::Drop;
    }
    if let Err(slot) = scratch.mate_drops.binary_search(&mate_idx) {
        scratch.mate_drops.insert(slot, mate_idx);
    }
    Overlap::Keep
}

/// What a column requires of a read before it counts: a called base, a base
/// quality, a query position, and the base-quality and read-end mask filters.
/// Mapping quality is the loader's (`RastairRecordFilter::filter_raw`): a read
/// below it never reaches a column.
///
/// One definition, applied to a read the column iterates *and* to a mate
/// consulted through it, so a mate can never contribute evidence that a read in
/// its own right would have been denied.
fn observed(
    view: &AlignmentView<'_, '_, RastairReadExtras>,
    params: &PileupMappingParams,
    reference_base: Base,
    context: &SequenceContext,
) -> Option<Observed> {
    let qual = view.qual()?;
    let baseq = qual.get()?;
    let base = view.base()?;
    let qpos = view.qpos()?;
    if !params.quality.keeps_base(qual) {
        return None;
    }
    if !passes_read_masking(view, reference_base, context) {
        return None;
    }
    Some(Observed { base, baseq, qpos })
}

/// What a read contributes at one column, once it has passed [`observed`].
#[derive(Clone, Copy)]
struct Observed {
    base: Base,
    baseq: u8,
    qpos: QPos,
}

/// Reusable buffers for one column, owned by the caller of
/// [`ColumnDraft::accumulate`] and reset there.
#[derive(Default)]
pub(crate) struct ColumnScratch {
    /// Right mates whose left mate won the overlap at this column, in ascending
    /// record index so membership is a binary search. Only the reads inside a
    /// mate overlap ever touch it; at high coverage there can be hundreds,
    /// which is why the buffer is the caller's and not an inline `SmallVec`
    /// that would spill to the heap once per column.
    mate_drops: Vec<RecordIdx>,
    mates: MateObservations,
    pub(super) haplotypes: HaplotypeScratch,
    /// `None` until the first cluster, and after it when the pair-HMM is
    /// scored on the CPU.
    pub(super) gpu: Option<super::gpu::Batch>,
}

impl ColumnScratch {
    fn begin(&mut self, depth: usize) {
        self.mate_drops.clear();
        self.mates.begin(depth);
    }
}

/// [`observed`] verdicts for reads the column was asked about before its walk
/// reached them, so no read is asked twice.
///
/// Three sites want the same answer for the same read: the read itself, the
/// *left* mate of an overlapping pair — which has to know what the column will
/// make of a read it has not reached yet — and, on the indel path, whichever
/// mate `counted_mate` reaches for. Measured over chr12, that is 1.19 asks per
/// read without indel calling and 1.38 with it, for ~90 instructions of
/// quality filtering, read-end masking and, for a rescued clip, a CpG-partner
/// test each time.
///
/// One slot per alignment, addressed by the read's position in the column
/// rather than by its record index, so both storing and claiming a verdict are
/// a single indexed access. `generation` is what makes that affordable: slots
/// are stale rather than cleared between columns, which is why the buffer is
/// only ever grown.
#[derive(Default)]
struct MateObservations {
    slots: Vec<MateObservation>,
    generation: u64,
}

#[derive(Clone, Copy, Default)]
struct MateObservation {
    /// The column this slot was written for; anything else means "empty".
    generation: u64,
    second_in_template: bool,
    observed: Option<Observed>,
}

impl MateObservations {
    fn begin(&mut self, depth: usize) {
        self.generation = self.generation.wrapping_add(1);
        if self.slots.len() < depth {
            self.slots.resize(depth, MateObservation::default());
        }
    }

    /// The verdict for the read at `slot`, resolving it against the column and
    /// remembering it when this is the first ask.
    ///
    /// `None` means the column has no entry at that position; a read that
    /// merely failed a filter is `Some` with `observed: None`.
    fn resolve(
        &mut self,
        column: &PileupColumn<'_, RastairReadExtras>,
        slot: usize,
        params: &PileupMappingParams,
        reference_base: Base,
        context: &SequenceContext,
    ) -> Option<MateObservation> {
        if let Some(known) = self.slots.get(slot).filter(|s| s.generation == self.generation) {
            return Some(*known);
        }
        let mate = column.alignment_at(slot)?;
        let entry = MateObservation {
            generation: self.generation,
            second_in_template: mate.flags.is_second_in_template(),
            observed: observed(&mate, params, reference_base, context),
        };
        *self.slots.get_mut(slot)? = entry;
        Some(entry)
    }

    /// The verdict already resolved for the read at `slot`, for the column walk.
    ///
    /// Slots are left in place: `counted_mate` can reach back to a read the
    /// walk has passed.
    fn resolved(&self, slot: usize) -> Option<&MateObservation> {
        self.slots.get(slot).filter(|s| s.generation == self.generation)
    }
}

/// The linked mate of `view` in this column, when it would count in its own
/// right.
///
/// Deduplication keeps one read per fragment, but the indel observation is a
/// property of the *fragment*. Reading it off the surviving read alone loses
/// whatever only the dropped mate showed.
///
/// The mate must clear [`observed`], the same bar the kept read cleared, so it
/// cannot contribute evidence a read in its own right would have been denied.
fn counted_mate<'a, 'eng>(
    column: &'a PileupColumn<'eng, RastairReadExtras>,
    view: &AlignmentView<'a, 'eng, RastairReadExtras>,
    params: &PileupMappingParams,
    reference_base: Base,
    context: &SequenceContext,
    mates: &mut MateObservations,
) -> Option<AlignmentView<'a, 'eng, RastairReadExtras>> {
    // Same reject `mate_of` applies before its search: outside the pair's
    // overlap the mate cannot be in this column.
    if !view.in_mate_overlap() {
        return None;
    }
    let slot = column.position_of(view.mate_idx()?)?;
    let entry = mates.resolve(column, slot, params, reference_base, context)?;
    entry.observed.is_some().then(|| column.alignment_at(slot)).flatten()
}

pub(super) fn passes_read_masking(
    view: &AlignmentView<'_, '_, RastairReadExtras>,
    reference_base: Base,
    context: &SequenceContext,
) -> bool {
    if view.is_soft_clip() {
        // A rescued fringe base is a read-end base by construction, so the
        // read-end mask would always reject it; the CpG-partner check is its
        // filter instead.
        let Some(observed) = view.base() else { return false };
        view.extra()
            .strand
            .is_some_and(|strand| soft_clip_cpg_partner(reference_base, observed, context, strand))
    } else {
        view.qpos().is_some_and(|qpos| view.extra().keeps(qpos.get()))
    }
}

/// Is the observed soft clipped base a CpG position?
fn soft_clip_cpg_partner(
    reference_base: Base,
    observed: Base,
    context: &SequenceContext,
    strand: Strand,
) -> bool {
    match (reference_base, strand) {
        (Base::C, Strand::OT) => {
            context.after_1 == Some(Base::G) && matches!(observed, Base::C | Base::T)
        }
        (Base::G, Strand::OB) => {
            context.before_1 == Some(Base::C) && matches!(observed, Base::G | Base::A)
        }
        _ => false,
    }
}

fn build_indel_observation(
    view: &AlignmentView<'_, '_, RastairReadExtras>,
    idx: usize,
    segment: &Segment,
) -> Option<IndelAllele> {
    match view.alignment().indel_after() {
        Indel::None => None,
        Indel::Insertion(_) => {
            let bases: SmallVec<Base, 4> = view.inserted_bases().iter().copied().collect();
            (!bases.is_empty()).then_some(IndelAllele::Insertion(bases))
        }
        Indel::Deletion(del_len) => {
            let ref_start = idx.checked_add(1)?;
            let ref_end = ref_start.checked_add(usize::try_from(del_len).ok()?)?;
            let bases: SmallVec<Base, 4> = segment
                .sequence
                .get(ref_start..ref_end)
                .map(|s| s.iter().map(|&b| Base::from(b)).collect())
                .unwrap_or_default();
            (!bases.is_empty()).then_some(IndelAllele::Deletion(bases))
        }
    }
}
