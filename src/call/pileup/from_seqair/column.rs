//! One column's reads, accumulated into a [`ColumnDraft`] and finished into
//! [`PileupMetrics`].
//!
//! [`ColumnDraft::accumulate`] counts a column's reads; [`ColumnDraft::finish`]
//! turns a draft the record filters keep into metrics. Nothing skips a column
//! with zero reads or an `N` reference (`Base::Unknown`) before it is
//! finished, so both are real paths: return default metrics, not an error.
//!
//! **Driving the engine.** `PileupEngine::new` takes a `PileupInput`, and only
//! `store.prepare_for_pileup()` makes one, sorted by position with mates
//! linked; without either, columns go missing or overlap dedup silently does
//! nothing. `PileupColumn::mate_of` gives a view's linked mate in the column,
//! all `counted_mate` needs; pairwise rules stay on this side of the boundary.
//! `PileupEngine::reclaim_allocation` returns an *empty* store that keeps its
//! capacity for the next region.
//!
//! **Indel observations.** Each rule guards a class of mistake:
//!
//! - A deletion's REF comes from the deletion site: `view.qpos()` is a read
//!   offset, not the column's genomic `pos` (seqair's `QPos` keeps them apart).
//! - Overlap dedup keeps one read per fragment and can keep the mate that does
//!   not span the indel, so the dropped mate stands in:
//!   `build_indel_observation(view).or_else(mate)`. seqair's `pair_indel`
//!   ("own indel wins, the mate is never consulted") is not used.
//! - The CIGAR counts only pick candidates: `IndelCounts` holds, per allele,
//!   the fragments whose alignment carries it, and nothing reads it for a
//!   call. There are no read-level indel filters; the pair-HMM weighs each
//!   read's alignment itself.

use crate::call::pileup::indels::{IndelAllele, IndelCounts};
use crate::{
    call::{PreFilterInputs, pair_hmm::HaplotypeScratch, process::PileupMappingParams},
    metrics::{
        Alt, AltFilters, Filters, FormsDenovo, PairedCounts, PerBaseAccumulators, PileupMetrics,
        ReadKey, RecordTags, order_alts,
    },
    sequence::{RastairReadExtras, Segment},
    utils::SequenceContext,
    vcf::TapsContext,
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use seqair::bam::{
    RecordIdx,
    pileup::{AlignmentView, Indel, PileupColumn, PileupOp},
};
use seqair_types::{Base, QPos, SmallVec, Strand, SumOfSquares};
use std::rc::Rc;
use tracing::{instrument, warn};

/// One column's reads, accumulated but not yet reduced to a [`PileupMetrics`].
///
/// The split exists so a column can be rejected before the expensive half runs:
/// [`ColumnDraft::finish`] pays nine divisions and square roots per allele, a
/// `PositionMetrics`, and a ~900-byte struct write, and roughly seven columns in
/// eight are then dropped by [`RecordFilters::pre_filter`]. Everything that
/// filter reads is already known here — see [`ColumnDraft::pre_filter_inputs`].
pub(crate) struct ColumnDraft {
    segment: Rc<Segment>,
    pos_u32: u32,
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

/// What a column contributes to its neighbours' de-novo adjacency, and nothing
/// else: [`crate::call::process::set_denovo_adj`] asks only whether the column
/// on one side carries an alt that would create the other half of a CpG.
#[derive(Debug, Clone, Copy)]
pub(crate) struct DenovoNeighbour {
    pos: u32,
    becomes_c: bool,
    becomes_g: bool,
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
    // Without an early rejection every column is kept, so none is plain.
    if params.early_reject.is_none() {
        return None;
    }
    let pos = u32::try_from(column.pos().as_u64()).ok()?;
    let idx = segment.pos_to_idx(pos).ok()?;
    let reference = Base::from(*segment.sequence.get(idx)?);
    let settled = match reference {
        Base::A | Base::T => true,
        Base::G => {
            let after_c = idx
                .checked_sub(1)
                .and_then(|i| segment.sequence.get(i))
                .is_some_and(|&b| Base::from(b) == Base::C);
            let rescued =
                previous.is_some_and(|p| p.becomes_c && p.pos.checked_add(1) == Some(pos));
            !after_c && !rescued
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
    plain.then_some(DenovoNeighbour { pos, becomes_c: false, becomes_g: false })
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
        let pos = column.pos().as_u64();
        let pos_u32 = u32::try_from(pos).wrap_err("pileup position exceeds u32")?;
        let idx = segment.pos_to_idx(pos_u32)?;
        let depth = column.depth();

        let reference_base: Base =
            segment.sequence.get(idx).wrap_err("failed to get reference base")?.into();

        let context =
            SequenceContext::new(idx, &segment).wrap_err("failed to get sequence context")?;

        scratch.begin(depth);

        let mut accumulators = PerBaseAccumulators::default();
        let mut pos_baseq = SumOfSquares::new();
        let mut pos_mapq = SumOfSquares::new();
        let mut mapq0: u32 = 0;
        let mut total_depth: usize = 0;
        let mut alt_bases: SmallVec<Base, 4> = SmallVec::new();
        let mut indel_observations = SmallVec::new();
        let mut before_counts = PairedCounts::default();
        let mut after_counts = PairedCounts::default();

        for (slot, view) in column.alignments().enumerate() {
            // Already resolved if the overlap rule reached forward for this
            // read while deciding its mate; otherwise this is its first ask.
            let resolved = match scratch.mates.resolved(slot) {
                Some(known) => known.observed,
                None => observed(&view, params),
            };
            let Some(Observed { base, baseq, qpos }) = resolved else {
                continue;
            };
            let strand = view.extra().strand;
            let overlap = drops_overlapping_mate(column, &view, base, params, scratch);
            if overlap == Overlap::Drop {
                continue;
            }
            total_depth += 1;

            let qual_sq = f64::from(baseq).powi(2);
            let mapq_sq = f64::from(view.mapq).powi(2);
            accumulators.accumulate_fields(
                base,
                crate::metrics::ReadFields {
                    qual_sq,
                    mapq_sq,
                    strand,
                    matching_bases: view.matching_bases,
                    indels: view.indel_bases,
                    pos_in_read: qpos.get(),
                    read_length: view.seq_len,
                },
            );
            pos_baseq.add_squared(qual_sq);
            pos_mapq.add_squared(mapq_sq);
            if view.mapq == 0 {
                mapq0 += 1;
            }
            if base.known_index().is_some() && base != reference_base && !alt_bases.contains(&base)
            {
                alt_bases.push(base);
            }

            if strand != Strand::Unknown {
                let seq = view.seq();
                if let Some(&adj) = qpos.checked_sub(1).and_then(|i| seq.get(i.as_usize())) {
                    before_counts.increment(ReadKey { strand, current: base, adj });
                }
                // An adjacent base is only adjacent in the alignment when nothing
                // is inserted or deleted between them.
                if matches!(view.alignment().indel_after(), Indel::None)
                    && let Some(&adj) = qpos.checked_add(1).and_then(|i| seq.get(i.as_usize()))
                {
                    after_counts.increment(ReadKey { strand, current: base, adj });
                }
            }

            if params.call_indels {
                let mate = counted_mate(column, &view, params, &mut scratch.mates);
                // A deduplicated column keeps one read per fragment, and that
                // read is not necessarily the one carrying the fragment's
                // indel: when the mates disagree the rule can drop exactly the
                // mate that spans it, and the fragment's only indel observation
                // would go with it. The mate then speaks for the fragment.
                let evidence = build_indel_observation(&view, pos, segment.as_ref())
                    .or_else(|| build_indel_observation(mate.as_ref()?, pos, segment.as_ref()));
                if let Some(obs) = evidence {
                    indel_observations.push(obs);
                }
            }
        }

        Ok(ColumnDraft {
            segment,
            pos_u32,
            idx,
            reference_base,
            context,
            accumulators,
            pos_baseq,
            pos_mapq,
            mapq0,
            total_depth,
            alt_bases,
            indel_observations,
            before_counts,
            after_counts,
        })
    }

    /// The draft of a position no read reaches, or `None` unless the
    /// reference has a CpG there: the only such positions any output reports.
    pub(crate) fn uncovered_cpg(pos: u64, segment: &Rc<Segment>) -> Result<Option<ColumnDraft>> {
        let pos_u32 = u32::try_from(pos).wrap_err("position exceeds u32")?;
        let idx = segment.pos_to_idx(pos_u32)?;
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
        if !*cpg {
            return Ok(None);
        }
        let context =
            SequenceContext::new(idx, segment).wrap_err("failed to get sequence context")?;
        Ok(Some(ColumnDraft {
            segment: segment.clone(),
            pos_u32,
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
        }))
    }

    fn taps_context(&self) -> TapsContext {
        TapsContext::new(self.reference_base, self.context.before_1, self.context.after_1)
    }

    /// The sequence index of this column, for the caller's sliding entropy.
    pub(crate) fn idx(&self) -> usize {
        self.idx
    }

    pub(crate) fn denovo_neighbour(&self) -> DenovoNeighbour {
        let mut becomes_c = false;
        let mut becomes_g = false;
        for &base in self.alt_bases.iter() {
            match self.taps_context().forms_cpg(base) {
                FormsDenovo::ThisBecomesC => becomes_c = true,
                FormsDenovo::ThisBecomesG => becomes_g = true,
                FormsDenovo::No => {}
            }
        }
        DenovoNeighbour { pos: self.pos_u32, becomes_c, becomes_g }
    }

    /// Only a reference `C` can still be rescued by the column after it, so
    /// every other draft's fate is settled the moment it is accumulated.
    pub(crate) fn awaits_successor(&self) -> bool {
        self.reference_base == Base::C
    }

    /// The three questions [`RecordFilters::pre_filter`] asks of a finished
    /// `PileupMetrics`, answered from the draft.
    ///
    /// `before`/`after` are the neighbouring *emitted* columns, exactly what
    /// `map_surrounding` would hand `set_denovo_adj`; each contributes only
    /// when it is genomically adjacent.
    pub(crate) fn pre_filter_inputs(
        &self,
        before: Option<DenovoNeighbour>,
        after: Option<DenovoNeighbour>,
    ) -> PreFilterInputs {
        let denovo_adj = (self.reference_base == Base::G
            && before.is_some_and(|b| b.becomes_c && b.pos.checked_add(1) == Some(self.pos_u32)))
            || (self.reference_base == Base::C
                && after
                    .is_some_and(|a| self.pos_u32.checked_add(1) == Some(a.pos) && a.becomes_g));
        let neighbour = self.denovo_neighbour();
        PreFilterInputs {
            has_alts: !self.alt_bases.is_empty(),
            cpg: *self.taps_context().in_cpg()
                || denovo_adj
                || neighbour.becomes_c
                || neighbour.becomes_g,
            has_indels: !self.indel_observations.is_empty(),
        }
    }

    /// Reduce the accumulated column to the metrics the rest of the pipeline
    /// consumes.
    pub(crate) fn finish(self) -> Result<PileupMetrics> {
        let ColumnDraft {
            segment,
            pos_u32,
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
        let pos_metrics = crate::metrics::PositionMetrics::new(
            total_depth,
            reference_base,
            context.before_1,
            context.after_1,
            pos_baseq.finish(depth),
            pos_mapq.finish(depth),
            mapq0,
        );

        let ref_metrics = if let Some(acc) = accumulators.take(reference_base) {
            acc.finish(reference_base, total_depth, pos_u32, reference_base, &context)?
        } else {
            crate::metrics::AlleleMetrics { base: reference_base, ..Default::default() }
        };

        let mut alts: SmallVec<Alt, 1> = alt_bases
            .iter()
            .map(|&base| {
                let acc = accumulators
                    .take(base)
                    .ok_or_else(|| color_eyre::eyre::eyre!("unknown base {base} in alt_bases"))?;
                let metrics = acc.finish(base, total_depth, pos_u32, reference_base, &context)?;
                Ok(Alt { base, metrics, filters: AltFilters::default(), call: Default::default() })
            })
            .collect::<Result<_>>()?;
        order_alts(&mut alts);

        let indel_data = if indel_observations.is_empty() {
            None
        } else {
            let counts = IndelCounts::aggregate(&indel_observations);
            Some(Box::new(crate::call::pileup::indels::IndelData { counts, scored: None }))
        };

        Ok(PileupMetrics {
            region: segment.range.clone(),
            pos: pos_u32,
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
/// it looks its mate up in the column, applies the same rule the name-based
/// collector used — drop the second-in-template read, or the later one when
/// both show the same base — and, when it wins, records the mate for the
/// `mate_drops` check the right mate then hits. Deciding at the left mate
/// keeps every kept read accumulating in column order.
fn drops_overlapping_mate(
    column: &PileupColumn<'_, RastairReadExtras>,
    view: &AlignmentView<'_, '_, RastairReadExtras>,
    base: Base,
    params: &PileupMappingParams,
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
    let Some(mate) = scratch.mates.resolve(column, slot, params) else {
        return Overlap::Keep;
    };
    let Some(Observed { base: mate_base, .. }) = mate.observed else {
        return Overlap::Keep;
    };

    // The name-based collector resolved the pair when it reached the *later*
    // read: it dropped that read if the bases agreed or it was read 2, and
    // otherwise dropped the earlier one.
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
) -> Option<Observed> {
    let baseq = view.qual()?.get()?;
    let base = view.base()?;
    let qpos = view.qpos()?;
    if baseq < params.quality.min_baseq {
        return None;
    }
    if !passes_read_masking(view) {
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
/// quality filtering and read-end masking each time.
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
    ) -> Option<MateObservation> {
        if let Some(known) = self.slots.get(slot).filter(|s| s.generation == self.generation) {
            return Some(*known);
        }
        let mate = column.alignment_at(slot)?;
        let entry = MateObservation {
            generation: self.generation,
            second_in_template: mate.flags.is_second_in_template(),
            observed: observed(&mate, params),
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
    mates: &mut MateObservations,
) -> Option<AlignmentView<'a, 'eng, RastairReadExtras>> {
    // Same reject `mate_of` applies before its search: outside the pair's
    // overlap the mate cannot be in this column.
    if !view.in_mate_overlap() {
        return None;
    }
    let slot = column.position_of(view.mate_idx()?)?;
    let entry = mates.resolve(column, slot, params)?;
    entry.observed.is_some().then(|| column.alignment_at(slot)).flatten()
}

pub(super) fn passes_read_masking(view: &AlignmentView<'_, '_, RastairReadExtras>) -> bool {
    view.qpos().is_some_and(|qpos| view.extra().mask.contains(&qpos.get()))
}

fn build_indel_observation(
    view: &AlignmentView<'_, '_, RastairReadExtras>,
    genomic_pos: u64,
    segment: &Segment,
) -> Option<IndelAllele> {
    match view.alignment().indel_after() {
        Indel::None => None,
        Indel::Insertion(_) => {
            let bases: SmallVec<Base, 4> = view.inserted_bases().iter().copied().collect();
            (!bases.is_empty()).then_some(IndelAllele::Insertion(bases))
        }
        Indel::Deletion(del_len) => {
            let ref_start = genomic_pos.checked_add(1)?.saturating_sub(segment.range.region.start);
            let ref_start = usize::try_from(ref_start).ok()?;
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
