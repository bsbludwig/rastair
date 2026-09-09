#![cfg(feature = "experimental-seqair")]

use super::{
    indels::{IndelAllele, IndelObservation},
    ref_features::{indel_ref_window_at, indel_tract_runs_at},
};
use crate::{
    call::{PreFilterInputs, process::PileupMappingParams},
    metrics::{
        Alt, AltFilters, Filters, FormsDenovo, PairedCounts, PerBaseAccumulators, PileupMetrics,
        ReadKey, RecordTags, aggregate_indels, alt_forms_denovo,
    },
    sequence::{RastairReadExtras, Segment},
    utils::SequenceContext,
    vcf::InCpG,
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use seqair::bam::pileup::{AlignmentView, Indel, PileupColumn};
use seqair_types::{Base, QPos, RmsAccumulator, SmallVec, Strand};
use std::rc::Rc;
use tracing::{debug, instrument, trace};

/// One column's reads, accumulated but not yet reduced to a [`PileupMetrics`].
///
/// The split exists so a column can be rejected before the expensive half runs:
/// [`ColumnDraft::finish`] pays nine divisions and square roots per allele, a
/// `PositionMetrics`, and a ~900-byte struct write, and roughly seven columns in
/// eight are then dropped by [`RecordFilters::pre_filter`]. Everything that
/// filter reads is already known here — see [`ColumnDraft::pre_filter_inputs`].
pub(crate) struct ColumnDraft {
    segment: Rc<Segment>,
    pos: u64,
    pos_u32: u32,
    idx: usize,
    reference_base: Base,
    context: SequenceContext,
    accumulators: PerBaseAccumulators,
    pos_baseq: RmsAccumulator,
    pos_mapq: RmsAccumulator,
    mapq0: u32,
    total_depth: usize,
    alt_bases: SmallVec<Base, 4>,
    indel_observations: SmallVec<IndelObservation, 3>,
    depth_offset: u32,
    soft_clip_count: u32,
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

impl ColumnDraft {
    /// `mate_drops` is scratch: a reusable buffer for the right mates this
    /// column drops. It is cleared here, and lives across columns only so a
    /// deep pileup does not allocate one per position.
    #[instrument(level = "trace", skip_all)]
    pub(crate) fn accumulate(
        column: &PileupColumn<'_, RastairReadExtras>,
        segment: Rc<Segment>,
        params: &PileupMappingParams,
        mate_drops: &mut Vec<u32>,
    ) -> Result<ColumnDraft> {
        let pos = column.pos().as_u64();
        let pos_u32 = u32::try_from(pos).wrap_err("pileup position exceeds u32")?;
        let idx = segment.pos_to_idx(pos_u32)?;
        let depth = column.depth();
        let max_reads = params.max_coverage.clamp(depth);
        if depth > max_reads {
            debug!(pos, depth, "Capping number of reads in pileup to {max_reads}");
        }

        let reference_base: Base =
            segment.sequence.get(idx).wrap_err("failed to get reference base")?.into();

        let context =
            SequenceContext::new(idx, &segment).wrap_err("failed to get sequence context")?;

        let dedup_overlaps = !params.keep_overlapping_reads;
        // Right mates whose left mate won the overlap at this column, in
        // ascending record index so membership is a binary search. Only the
        // reads inside a mate overlap ever touch it; at high coverage there can
        // be hundreds, which is why the buffer is the caller's and not an
        // inline `SmallVec` that would spill to the heap once per column.
        mate_drops.clear();

        let mut accumulators = PerBaseAccumulators::default();
        let mut pos_baseq = RmsAccumulator::new();
        let mut pos_mapq = RmsAccumulator::new();
        let mut mapq0: u32 = 0;
        let mut total_depth: usize = 0;
        let mut alt_bases: SmallVec<Base, 4> = SmallVec::new();
        let mut indel_observations = SmallVec::new();
        let mut depth_offset: u32 = 0;
        let mut soft_clip_count: u32 = 0;
        let mut before_counts = PairedCounts::default();
        let mut after_counts = PairedCounts::default();

        for view in column.alignments() {
            if total_depth >= max_reads {
                break;
            }
            let Some(Observed { base, baseq, qpos }) =
                observed(&view, params, reference_base, &context)
            else {
                continue;
            };
            let strand = view.extra().strand;
            if dedup_overlaps
                && drops_overlapping_mate(
                    column,
                    &view,
                    base,
                    params,
                    reference_base,
                    &context,
                    mate_drops,
                )
            {
                continue;
            }
            total_depth += 1;

            let qual_sq = f64::from(baseq).powi(2);
            let mapq_sq = f64::from(view.mapq).powi(2);
            accumulators.accumulate_fields(
                base,
                qual_sq,
                mapq_sq,
                strand,
                view.matching_bases,
                view.indel_bases,
                qpos.get(),
                view.seq_len,
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
                let aln = view.alignment();
                let extras = view.extra();
                // Every count below describes the *fragment*, not the read that
                // happened to survive dedup, so the dropped mate's shape has to
                // be OR'd in — `from_hts` accumulates these over both mates
                // before it deduplicates.
                let mate = dedup_overlaps
                    .then(|| counted_mate(column, &view, params, reference_base, &context))
                    .flatten();
                let mate_extras = mate.as_ref().map(|mate| mate.extra());
                let has_soft_clip =
                    extras.has_soft_clip || mate_extras.is_some_and(|extras| extras.has_soft_clip);
                let has_repeat =
                    extras.has_repeat || mate_extras.is_some_and(|extras| extras.has_repeat);

                if has_soft_clip {
                    soft_clip_count += 1;
                }

                // A deduplicated column keeps one read per fragment, and that
                // read is not necessarily the one carrying the fragment's
                // indel: when the mates disagree the rule can drop exactly the
                // mate that spans it, and the fragment's only indel observation
                // would go with it. `from_hts` never had this problem because
                // it votes per fragment *before* deduplicating — every mate
                // gets a turn, and the first surviving observation is the
                // fragment's.
                //
                // So the fallback is on the *observation*, not on the indel:
                // this read's own indel may exist and still be rejected by
                // `build_indel_observation`'s end-of-read and mismatch filters,
                // and the mate then still speaks for the fragment. Testing the
                // indel instead (seqair's `PileupColumn::pair_indel`, whose
                // rule is "own indel wins, mate never consulted") silently
                // drops those fragments — it cannot know about filters the
                // caller owns.
                let evidence = build_indel_observation(&view, pos, segment.as_ref(), params)
                    .or_else(|| {
                        let mate = mate.as_ref()?;
                        build_indel_observation(mate, pos, segment.as_ref(), params)
                    });
                // The noisy-reference count and the alternate side are two
                // sides of one split, and `IndelCounts::clean_depth` subtracts
                // a fragment counted on both twice. So a fragment is
                // noisy-reference only when *this* read shows no indel (the
                // alignment shape it slipped from is this read's) and the
                // fragment contributed no observation at all — the pair's
                // verdict, not the kept read's. `from_hts` reaches the same
                // rule through per-fragment votes.
                if matches!(aln.indel_after(), Indel::None)
                    && (has_repeat || has_soft_clip)
                    && evidence.is_none()
                {
                    depth_offset += 1;
                }

                if let Some(obs) = evidence {
                    indel_observations.push(obs);
                }
            }
        }

        Ok(ColumnDraft {
            segment,
            pos,
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
            depth_offset,
            soft_clip_count,
            before_counts,
            after_counts,
        })
    }

    /// The sequence index of this column, for the caller's sliding entropy.
    pub(crate) fn idx(&self) -> usize {
        self.idx
    }

    pub(crate) fn denovo_neighbour(&self) -> DenovoNeighbour {
        let mut becomes_c = false;
        let mut becomes_g = false;
        for &base in self.alt_bases.iter() {
            match alt_forms_denovo(base, self.reference_base, &self.context) {
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
            cpg: *InCpG::new(self.reference_base, self.context.before_1, self.context.after_1)
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
            pos,
            pos_u32,
            idx,
            reference_base,
            context,
            mut accumulators,
            pos_baseq,
            pos_mapq,
            mapq0,
            total_depth,
            alt_bases,
            indel_observations,
            depth_offset,
            soft_clip_count,
            before_counts,
            after_counts,
        } = self;

        let pos_metrics = crate::metrics::PositionMetrics::new(
            total_depth,
            reference_base,
            context.before_1,
            context.after_1,
            pos_baseq.finish(),
            pos_mapq.finish(),
            mapq0,
        );

        let ref_metrics = if let Some(acc) = accumulators.take(reference_base) {
            acc.finish(reference_base, total_depth, pos_u32, reference_base, &context)?
        } else {
            crate::metrics::AlleleMetrics { base: reference_base, ..Default::default() }
        };

        let alts = alt_bases
            .iter()
            .map(|&base| {
                let acc = accumulators
                    .take(base)
                    .ok_or_else(|| color_eyre::eyre::eyre!("unknown base {base} in alt_bases"))?;
                let metrics = acc.finish(base, total_depth, pos_u32, reference_base, &context)?;
                Ok(Alt { base, metrics, filters: AltFilters::default(), call: Default::default() })
            })
            .collect::<Result<_>>()?;

        let indel_data = if indel_observations.is_empty() {
            None
        } else {
            let counts = aggregate_indels(&indel_observations, total_depth, depth_offset, pos_u32);
            let (indel_ref_window, indel_ref_anchor) = indel_ref_window_at(idx, &segment);
            let runs = indel_tract_runs_at(pos, &segment);
            Some(Box::new(crate::call::pileup::indels::IndelData {
                observations: indel_observations,
                ref_window: indel_ref_window,
                ref_anchor: indel_ref_anchor,
                homopolymer_run: runs.homopolymer,
                dinucleotide_run: runs.dinucleotide,
                soft_clip_count,
                counts,
                calls: Vec::new(),
            }))
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
        })
    }
}

impl PileupMetrics {
    /// Accumulate and finish one column in one go.
    ///
    /// Production splits these two halves so a column can be rejected between
    /// them; the tests want every column, and this is the form they ask for.
    #[cfg(test)]
    pub(crate) fn from_seqair(
        column: &PileupColumn<'_, RastairReadExtras>,
        segment: Rc<Segment>,
        params: &PileupMappingParams,
        mate_drops: &mut Vec<u32>,
    ) -> Result<PileupMetrics> {
        ColumnDraft::accumulate(column, segment, params, mate_drops)?.finish()
    }
}

/// Decide the fate of a read that shares a mate overlap with another read in
/// this column, returning `true` when *this* read is the one to drop.
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
    reference_base: Base,
    context: &SequenceContext,
    mate_drops: &mut Vec<u32>,
) -> bool {
    if !view.in_mate_overlap() {
        return false;
    }
    let this_idx = view.alignment().record_idx();
    let Some(mate_idx) = view.alignment().mate_idx() else { return false };

    if this_idx > mate_idx {
        return mate_drops.binary_search(&this_idx).is_ok();
    }

    // A mate that is absent from this column, or that fails a filter here,
    // never formed a pair — this read stands on its own.
    let Some(mate) = column.find_record(mate_idx) else { return false };
    let Some(Observed { base: mate_base, .. }) = observed(&mate, params, reference_base, context)
    else {
        return false;
    };

    // The name-based collector resolved the pair when it reached the *later*
    // read: it dropped that read if the bases agreed or it was read 2, and
    // otherwise dropped the earlier one.
    if base != mate_base && !mate.flags.is_second_in_template() {
        return true;
    }
    if let Err(slot) = mate_drops.binary_search(&mate_idx) {
        mate_drops.insert(slot, mate_idx);
    }
    false
}

/// What a column requires of a read before it counts: a called base, a base
/// quality, a query position, and the mapping/base-quality and read-end mask
/// filters.
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
    let baseq = view.qual()?.get()?;
    let base = view.base()?;
    let qpos = view.qpos()?;
    if !params.quality.filter_fields(view.mapq, baseq) {
        return None;
    }
    if !passes_read_masking(view, reference_base, context) {
        return None;
    }
    Some(Observed { base, baseq, qpos })
}

/// What a read contributes at one column, once it has passed [`observed`].
struct Observed {
    base: Base,
    baseq: u8,
    qpos: QPos,
}

/// The linked mate of `view` in this column, when it would count in its own
/// right.
///
/// Deduplication keeps one read per fragment, but the soft-clip and
/// noisy-reference counts are properties of the *fragment*: `from_hts` votes
/// per fragment before it deduplicates, so whatever either mate shows counts
/// once. Reading them off the surviving read alone loses whatever only the
/// dropped mate showed.
///
/// The mate must clear [`observed`], the same bar the kept read cleared, so it
/// cannot contribute evidence a read in its own right would have been denied.
fn counted_mate<'a, 'eng>(
    column: &'a PileupColumn<'eng, RastairReadExtras>,
    view: &AlignmentView<'a, 'eng, RastairReadExtras>,
    params: &PileupMappingParams,
    reference_base: Base,
    context: &SequenceContext,
) -> Option<AlignmentView<'a, 'eng, RastairReadExtras>> {
    let mate = column.mate_of(view)?;
    observed(&mate, params, reference_base, context).is_some().then_some(mate)
}

fn passes_read_masking(
    view: &AlignmentView<'_, '_, RastairReadExtras>,
    reference_base: Base,
    context: &SequenceContext,
) -> bool {
    if view.is_soft_clip() {
        // A rescued fringe base is a read-end base by construction, so the
        // read-end mask would always reject it; the CpG-partner check is its
        // filter instead.
        let Some(observed) = view.base() else { return false };
        soft_clip_cpg_partner(reference_base, observed, context, view.extra().strand)
    } else {
        let Some(qpos) = view.qpos() else {
            return false;
        };
        view.extra().mask.contains(&qpos.get())
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
    genomic_pos: u64,
    segment: &Segment,
    params: &PileupMappingParams,
) -> Option<IndelObservation> {
    let aln = view.alignment();
    let extras = view.extra();

    let indel = aln.indel_after();
    if matches!(indel, Indel::None) {
        return None;
    }

    let qpos = aln.qpos()?;
    let qpos_usize = qpos.as_usize();
    let read_len = aln.seq_len as usize;
    let indel_cutoff = params.indel_end_of_read_cutoff;

    if qpos_usize < indel_cutoff || qpos_usize >= read_len.saturating_sub(indel_cutoff) {
        trace!(qpos = qpos.get(), read_len, "Indel skipped: too close to read end");
        return None;
    }
    if extras.taps_aware_mismatches > params.indel_max_mismatches {
        trace!(
            mismatches = extras.taps_aware_mismatches,
            max = params.indel_max_mismatches,
            "Indel skipped: too many non-TAPS mismatches"
        );
        return None;
    }

    let segment_start = segment.range.region.start as usize;
    let (allele, insertion_base_quals, post_del_base_qual) = match indel {
        Indel::Insertion(_) => {
            let bases: SmallVec<Base, 4> = view.inserted_bases().iter().copied().collect();
            if bases.is_empty() {
                return None;
            }
            let quals: SmallVec<u8, 4> =
                view.inserted_quals().iter().filter_map(|q| q.get()).collect();
            (IndelAllele::Insertion(bases), quals, 0)
        }
        Indel::Deletion(del_len) => {
            let ref_start = (genomic_pos as usize + 1).saturating_sub(segment_start);
            let ref_end = ref_start + del_len as usize;
            let bases: SmallVec<Base, 4> = segment
                .sequence
                .get(ref_start..ref_end)
                .map(|s| s.iter().map(|&b| Base::from(b)).collect())
                .unwrap_or_default();
            if bases.is_empty() {
                return None;
            }
            let post_del = qpos
                .checked_add(1)
                .and_then(|i| view.qualities().get(i.as_usize()))
                .and_then(|q| q.get())
                .unwrap_or(0);
            (IndelAllele::Deletion(bases), SmallVec::new(), post_del)
        }
        Indel::None => unreachable!("matched a non-None indel above"),
    };

    let base_qual = view.qualities().get(qpos.as_usize()).and_then(|q| q.get()).unwrap_or(0);

    Some(IndelObservation {
        allele,
        strand: extras.strand,
        reverse: aln.flags.is_reverse(),
        pos_in_read: qpos.get(),
        read_length: aln.seq_len,
        mapq: aln.mapq,
        base_qual,
        matching_bases: aln.matching_bases,
        num_indels_in_read: aln.indel_bases,
        insertion_base_quals,
        post_del_base_qual,
        has_repeat: extras.has_repeat,
        noisy: extras.has_repeat || extras.has_soft_clip,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::process::PileupMappingParams;
    use crate::call::variant_calling::{ReadMaskParams, ReadMaskSetting};
    use crate::sequence::{ChunkRegion, Region, Segment};
    use crate::utils::default;
    use seqair::bam::cigar::{CigarOp, CigarOpType};
    use seqair::bam::pileup::PileupEngine;
    use seqair::bam::record_store::{CustomizeRecordStore, RecordStore, SlimRecord};
    use seqair_types::{BamFlags, Pos0, Strand};

    /// The parts of `RastairRecordFilter::compute` these tests need: strand
    /// from the flags, soft-clip detection, and the read-end mask window.
    #[derive(Default, Clone)]
    struct TestExtras(ReadMaskParams);

    impl CustomizeRecordStore for TestExtras {
        type Extra = RastairReadExtras;

        fn compute(
            &mut self,
            rec: &SlimRecord,
            store: &RecordStore<RastairReadExtras>,
        ) -> RastairReadExtras {
            let has_soft_clip = rec
                .cigar(store)
                .map(|ops| ops.iter().any(|op| op.op_type() == CigarOpType::SoftClip))
                .unwrap_or(false);
            let strand = Strand::from(rec.flags);
            RastairReadExtras {
                strand,
                has_soft_clip,
                has_repeat: false,
                taps_aware_mismatches: 0,
                mask: self
                    .0
                    .keep_range(strand, rec.flags.is_reverse(), rec.seq_len)
                    .unwrap_or(0..0),
            }
        }
    }

    /// A read to push into a test store: everything the overlap dedup looks at.
    #[derive(Clone, Debug)]
    struct TestRead {
        qname: Vec<u8>,
        pos: u32,
        flags: u16,
        bases: Vec<Base>,
        quals: Vec<u8>,
        mapq: u8,
        cigar: Vec<CigarOp>,
        mate_pos: i32,
    }

    impl TestRead {
        /// A plain `<len>M` read whose every base is `base`.
        fn matching(qname: &[u8], pos: u32, len: usize, base: Base, flags: u16) -> Self {
            Self {
                qname: qname.to_vec(),
                pos,
                flags,
                bases: vec![base; len],
                quals: vec![40; len],
                mapq: 60,
                cigar: vec![CigarOp::new(CigarOpType::Match, len as u32)],
                mate_pos: -1,
            }
        }

        fn with_base_at(mut self, offset: usize, base: Base) -> Self {
            if let Some(slot) = self.bases.get_mut(offset) {
                *slot = base;
            }
            self
        }

        fn with_qual(mut self, qual: u8) -> Self {
            self.quals = vec![qual; self.bases.len()];
            self
        }

        fn with_mapq(mut self, mapq: u8) -> Self {
            self.mapq = mapq;
            self
        }

        fn ref_span(&self) -> u32 {
            self.cigar
                .iter()
                .filter(|op| matches!(op.op_type(), CigarOpType::Match | CigarOpType::Deletion))
                .map(|op| op.len())
                .sum()
        }
    }

    /// Push `reads` as one mate pair per qname, wire up the mate fields, and
    /// link — what `Readers::pileup` does for a real fetch.
    fn store_of(reads: &[TestRead], masking: &ReadMaskParams) -> RecordStore<RastairReadExtras> {
        let mut extras = TestExtras(masking.clone());
        let mut store = RecordStore::<RastairReadExtras>::new();
        for read in reads {
            let mate_pos = reads
                .iter()
                .find(|other| other.qname == read.qname && other.pos != read.pos)
                .map_or(read.mate_pos, |other| other.pos as i32);
            let end = read.pos + read.ref_span().max(1) - 1;
            store
                .push_fields(
                    Pos0::new(read.pos).unwrap(),
                    Pos0::new(end).unwrap(),
                    BamFlags::from(read.flags),
                    read.mapq,
                    read.bases.len() as u32,
                    0,
                    &read.qname,
                    &read.cigar,
                    &read.bases,
                    &read.quals,
                    &[],
                    0,
                    0,
                    mate_pos,
                    0,
                    &mut extras,
                )
                .unwrap();
        }
        // Linking (and sorting) is `prepare_for_pileup`'s job at the point the
        // engine is built, so this fixture no longer has to remember either —
        // reads may be listed in any order.
        store
    }

    /// The overlapping-pair rule as it was before mate links: group the
    /// column's surviving alignments by qname and, on the second one, drop
    /// whichever `resolve_pair` chose. An independent oracle — it reaches the
    /// answer by matching names, which is exactly what the new code does not do.
    fn kept_by_name_collector(
        column: &PileupColumn<'_, RastairReadExtras>,
        params: &PileupMappingParams,
        reference_base: Base,
        context: &SequenceContext,
    ) -> Vec<usize> {
        let passing: Vec<usize> = column
            .alignments()
            .enumerate()
            .filter_map(|(idx, view)| {
                let baseq = view.qual()?.get()?;
                let _base = view.base()?;
                view.qpos()?;
                if !params.quality.filter_fields(view.mapq, baseq) {
                    return None;
                }
                if !passes_read_masking(&view, reference_base, context) {
                    return None;
                }
                Some(idx)
            })
            .collect();

        let alignments: Vec<AlignmentView<'_, '_, RastairReadExtras>> =
            column.alignments().collect();
        let mut first_by_name: Vec<(&[u8], usize)> = Vec::new();
        let mut removed: Vec<usize> = Vec::new();
        for &idx in &passing {
            let view = &alignments[idx];
            let name = view.qname();
            match first_by_name.iter().find(|(seen, _)| *seen == name) {
                None => first_by_name.push((name, idx)),
                Some(&(_, other)) => {
                    let this_base = view.base();
                    let other_base = alignments[other].base();
                    if this_base == other_base || view.flags.is_second_in_template() {
                        removed.push(idx);
                    } else {
                        removed.push(other);
                    }
                }
            }
        }
        passing.into_iter().filter(|idx| !removed.contains(idx)).collect()
    }

    /// Run both implementations over every column of `reads` and require the
    /// same surviving observations: total depth, and per allele the depth and
    /// the OT/OB split.
    fn assert_same_as_name_collector(reads: &[TestRead], seq: &[u8], params: &PileupMappingParams) {
        let seg = segment(seq);
        let store = store_of(reads, &params.read_masking);
        let last = seq.len().saturating_sub(1) as u32;
        let mut engine = PileupEngine::new(
            store.prepare_for_pileup().input,
            Pos0::new(0).unwrap(),
            Pos0::new(last).unwrap(),
        );
        if params.rescue_soft_clip_cpg {
            engine.set_soft_clip_overhang(1);
        }
        if let Some(cap) = params.max_coverage.per_column() {
            engine.set_max_depth(cap);
        }

        let mut scratch = Vec::new();
        let mut columns = 0;
        while let Some(col) = engine.pileups() {
            let pos = col.pos().as_u64() as usize;
            let reference_base = Base::from(*seq.get(pos).unwrap());
            let context = SequenceContext::new(pos, &seg).unwrap();
            let expected = kept_by_name_collector(&col, params, reference_base, &context);

            let alignments: Vec<AlignmentView<'_, '_, RastairReadExtras>> =
                col.alignments().collect();
            let mut expected_depth = 0u32;
            let mut expected_by_base: Vec<(Base, u32, u32)> = Vec::new();
            for &idx in &expected {
                let view = &alignments[idx];
                let Some(base) = view.base() else { continue };
                expected_depth += 1;
                let strand = view.extra().strand;
                let slot = match expected_by_base.iter_mut().find(|(b, _, _)| *b == base) {
                    Some(slot) => slot,
                    None => {
                        expected_by_base.push((base, 0, 0));
                        expected_by_base.last_mut().unwrap()
                    }
                };
                match strand {
                    Strand::OT => slot.1 += 1,
                    Strand::OB => slot.2 += 1,
                    Strand::Unknown => {}
                }
            }

            let pm = PileupMetrics::from_seqair(&col, seg.clone(), params, &mut scratch).unwrap();
            assert_eq!(
                pm.pos_metrics.depth as u32, expected_depth,
                "pos {pos}: depth differs from the name-collector rule"
            );
            for (base, ot, ob) in expected_by_base {
                let metrics = if base == pm.reference_base {
                    &pm.ref_metrics
                } else {
                    pm.alts
                        .iter()
                        .find(|alt| alt.base == base)
                        .map(|alt| &alt.metrics)
                        .unwrap_or_else(|| panic!("pos {pos}: no alt for {base}"))
                };
                assert_eq!(metrics.depth, ot + ob, "pos {pos}, {base}: allele depth");
                assert_eq!(metrics.strand_count.ot, ot, "pos {pos}, {base}: OT count");
                assert_eq!(metrics.strand_count.ob, ob, "pos {pos}, {base}: OB count");
            }
            columns += 1;
        }
        assert!(columns > 0, "no columns produced");
    }

    fn segment(seq: &[u8]) -> Rc<Segment> {
        segment_at(0, seq)
    }

    /// A segment placed at `start` on the contig. Everything a segment resolves
    /// against the reference goes through `pos - start`, so a fixture at 0
    /// cannot tell a genomic position from a read offset — see
    /// [`deletion_allele_reads_the_reference_at_the_deletion_site`].
    fn segment_at(start: u64, seq: &[u8]) -> Rc<Segment> {
        let end = start + seq.len().saturating_sub(1) as u64;
        Rc::new(Segment {
            range: std::sync::Arc::new(ChunkRegion {
                region: Region { contig: "chr1".into(), start, end },
                last_position: start + seq.len() as u64,
                overlap_start: 0,
                overlap_end: 0,
            }),
            sequence: seq.to_vec(),
            overlap_start: 0,
            overlap_end: 0,
        })
    }

    /// Overlap dedup keeps one read per fragment, and the rule can keep the
    /// mate that does *not* span the indel. The fragment's indel evidence must
    /// survive that: `from_hts` votes per fragment before deduplicating, so it
    /// never lost one, and the seqair path has to reach the same answer through
    /// the column's `pair_indel`.
    ///
    /// Fixture: mates agree on the anchor base, so the rule drops the right
    /// mate — which is the one carrying the deletion.
    #[test]
    fn a_dropped_mate_does_not_take_the_fragments_indel_with_it() {
        // Segment offsets:  0123456789
        let seq = b"AAAACGTTGG";
        let seg = segment(seq);
        let params =
            PileupMappingParams { call_indels: true, indel_end_of_read_cutoff: 0, ..default() };

        // Left mate, 6M at 0, no indel. Right mate, 2M2D4M at 2: it deletes the
        // segment's C,G at offsets 4-5, anchored at position 3.
        let reads = [
            TestRead {
                qname: b"pair".to_vec(),
                pos: 0,
                flags: 99,
                bases: vec![Base::A, Base::A, Base::A, Base::A, Base::C, Base::G],
                quals: vec![40; 6],
                mapq: 60,
                cigar: vec![CigarOp::new(CigarOpType::Match, 6)],
                mate_pos: 2,
            },
            TestRead {
                qname: b"pair".to_vec(),
                pos: 2,
                flags: 147,
                bases: vec![Base::A, Base::A, Base::T, Base::T, Base::G, Base::G],
                quals: vec![40; 6],
                mapq: 60,
                cigar: vec![
                    CigarOp::new(CigarOpType::Match, 2),
                    CigarOp::new(CigarOpType::Deletion, 2),
                    CigarOp::new(CigarOpType::Match, 4),
                ],
                mate_pos: 0,
            },
        ];

        let alleles_at_anchor = |params: &PileupMappingParams| {
            let store = store_of(&reads, &params.read_masking);
            let mut engine = PileupEngine::new(
                store.prepare_for_pileup().input,
                Pos0::new(0).unwrap(),
                Pos0::new(9).unwrap(),
            );
            let mut scratch = Vec::new();
            let mut out = Vec::new();
            while let Some(col) = engine.pileups() {
                if col.pos().as_u64() != 3 {
                    continue;
                }
                let pm =
                    PileupMetrics::from_seqair(&col, seg.clone(), params, &mut scratch).unwrap();
                out = pm
                    .indel_data
                    .map(|d| d.observations.iter().map(|o| o.allele.clone()).collect())
                    .unwrap_or_default();
            }
            out
        };

        let deletion = IndelAllele::Deletion([Base::C, Base::G].into_iter().collect());
        assert_eq!(
            alleles_at_anchor(&params),
            vec![deletion.clone()],
            "the deduplicated column must still report the dropped mate's deletion"
        );

        // With dedup off both mates are their own fragment, and only one of
        // them spans the deletion — so the count is unchanged, not doubled.
        let mut kept =
            PileupMappingParams { call_indels: true, indel_end_of_read_cutoff: 0, ..default() };
        kept.variant_calling.keep_overlapping_reads = true;
        assert_eq!(
            alleles_at_anchor(&kept),
            vec![deletion],
            "keeping both mates must not turn one fragment's indel into two"
        );
    }

    /// The mate speaks for the fragment when the kept read's *observation* is
    /// missing — not only when its indel is.
    ///
    /// `build_indel_observation` rejects an indel too close to a read end or on
    /// a read with too many non-TAPS mismatches, and those are filters the
    /// caller owns. A rule keyed on the indel alone (seqair's
    /// `PileupColumn::pair_indel`: "own indel wins, the mate is never
    /// consulted") therefore short-circuits on a read whose own indel is about
    /// to be thrown away, and the fragment contributes nothing. `from_hts` gets
    /// this right by giving every mate a turn and taking the first observation
    /// that survives.
    ///
    /// Fixture: the mates disagree on the anchor base and the *right* mate is
    /// first-in-template, so the rule drops the left one — and the surviving
    /// read reaches the deletion 1 bp into itself, inside the cutoff, while the
    /// dropped mate reaches it 3 bp in.
    #[test]
    fn the_mate_speaks_when_the_kept_reads_observation_is_filtered() {
        // Segment offsets:  0123456789
        let seq = b"AAAACGTTGG";
        let seg = segment(seq);
        let params =
            PileupMappingParams { call_indels: true, indel_end_of_read_cutoff: 2, ..default() };

        let reads = [
            // Dropped by the rule; reaches the anchor at qpos 3, outside the cutoff.
            TestRead {
                qname: b"pair".to_vec(),
                pos: 0,
                flags: 147,
                bases: vec![Base::A, Base::A, Base::A, Base::A, Base::T, Base::T, Base::G, Base::G],
                quals: vec![40; 8],
                mapq: 60,
                cigar: vec![
                    CigarOp::new(CigarOpType::Match, 4),
                    CigarOp::new(CigarOpType::Deletion, 2),
                    CigarOp::new(CigarOpType::Match, 4),
                ],
                mate_pos: 2,
            },
            // Kept; reaches the anchor at qpos 1, inside the cutoff.
            TestRead {
                qname: b"pair".to_vec(),
                pos: 2,
                flags: 99,
                bases: vec![Base::A, Base::G, Base::T, Base::T, Base::G, Base::G],
                quals: vec![40; 6],
                mapq: 60,
                cigar: vec![
                    CigarOp::new(CigarOpType::Match, 2),
                    CigarOp::new(CigarOpType::Deletion, 2),
                    CigarOp::new(CigarOpType::Match, 4),
                ],
                mate_pos: 0,
            },
        ];

        let store = store_of(&reads, &params.read_masking);
        let mut engine = PileupEngine::new(
            store.prepare_for_pileup().input,
            Pos0::new(0).unwrap(),
            Pos0::new(9).unwrap(),
        );
        let mut scratch = Vec::new();
        let mut alleles: Vec<IndelAllele> = Vec::new();
        while let Some(col) = engine.pileups() {
            if col.pos().as_u64() != 3 {
                continue;
            }
            let pm = PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut scratch).unwrap();
            alleles = pm
                .indel_data
                .map(|d| d.observations.iter().map(|o| o.allele.clone()).collect())
                .unwrap_or_default();
        }

        assert_eq!(
            alleles,
            vec![IndelAllele::Deletion([Base::C, Base::G].into_iter().collect())],
            "the dropped mate's surviving observation must stand in for the kept read's"
        );
    }

    /// A deletion's REF bases must come from the reference *at the deletion*,
    /// which is `genomic_pos + 1 - segment.start` into the segment's sequence.
    ///
    /// The regression this guards: the loop's genomic `pos` was shadowed by
    /// `view.qpos()`, a read-local offset, so the subtraction saturated to 0 and
    /// the allele became the segment's opening bases — real-looking sequence
    /// from the wrong locus. A segment starting at 0 cannot see this, because
    /// there the clamp and the correct answer coincide for a read at offset 0.
    #[test]
    fn deletion_allele_reads_the_reference_at_the_deletion_site() {
        const START: u64 = 1000;
        // Segment sequence, offset:  0123456789
        //                            AAAACGTTGG
        // A read matches 4 bases from 1000, deletes CG at 1004-1005, matches on.
        // The anchor column is 1003 (the base before the deletion), so the
        // deleted bases are the segment's offsets 4..6 = C,G — and *not* its
        // opening A,A, which is what the read-offset bug produced.
        let seq = b"AAAACGTTGG";
        let seg = segment_at(START, seq);
        let mut params = PileupMappingParams::default();
        params.call_indels = true;
        params.indel_end_of_read_cutoff = 0;

        let mut extras = TestExtras(params.read_masking.clone());
        let mut store = RecordStore::<RastairReadExtras>::new();
        store
            .push_fields(
                Pos0::new(START as u32).unwrap(),
                Pos0::new(START as u32 + 9).unwrap(),
                BamFlags::from(99u16),
                60,
                8,
                0,
                b"deleter",
                &[
                    CigarOp::new(CigarOpType::Match, 4),
                    CigarOp::new(CigarOpType::Deletion, 2),
                    CigarOp::new(CigarOpType::Match, 4),
                ],
                &[Base::A, Base::A, Base::A, Base::A, Base::T, Base::T, Base::G, Base::G],
                &[40u8; 8],
                &[],
                0,
                -1,
                0,
                0,
                &mut extras,
            )
            .unwrap();

        let mut engine = PileupEngine::new(
            store.prepare_for_pileup().input,
            Pos0::new(START as u32).unwrap(),
            Pos0::new(START as u32 + 9).unwrap(),
        );
        let mut allele = None;
        while let Some(col) = engine.pileups() {
            if col.pos().as_u64() != START + 3 {
                continue;
            }
            let metrics =
                PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut Vec::new()).unwrap();
            let data = metrics.indel_data.expect("the anchor column carries the deletion");
            let obs = data.observations.first().expect("one deletion observation").clone();
            allele = Some(obs.allele);
        }

        assert_eq!(
            allele,
            Some(IndelAllele::Deletion([Base::C, Base::G].into_iter().collect())),
            "deleted bases must be read at the deletion site, not at the segment start"
        );
    }

    /// A read aligned to the G of a CpG with its leading base (a methylated T
    /// over the C) soft-clipped is rescued into the pileup at the C: with the
    /// engine overhang on, the T appears as an OT alt at the CpG-C; with it off,
    /// the C position has no column at all.
    #[test]
    fn rescues_soft_clipped_cpg_partner() {
        // Reference: T T C G T T — CpG is C@2 / G@3.
        let seg = segment(b"TTCGTT");
        let params = PileupMappingParams::default();
        let mut extras = TestExtras(params.read_masking.clone());

        let mut build_store = || {
            let mut store = RecordStore::<RastairReadExtras>::new();
            // 1S 3M at pos 3: clip base T over ref C@2, aligned G,T,T over 3,4,5.
            // flags 99 = paired/proper/mate-reverse/first → OT.
            store
                .push_fields(
                    Pos0::new(3).unwrap(),
                    Pos0::new(5).unwrap(),
                    BamFlags::from(99u16),
                    60,
                    3,
                    0,
                    b"clipped",
                    &[CigarOp::new(CigarOpType::SoftClip, 1), CigarOp::new(CigarOpType::Match, 3)],
                    &[Base::T, Base::G, Base::T, Base::T],
                    &[40u8; 4],
                    &[],
                    0,
                    -1,
                    0,
                    0,
                    &mut extras,
                )
                .unwrap();
            store
        };

        let mut metrics_at = |overhang: u32| -> Option<PileupMetrics> {
            let mut engine = PileupEngine::new(
                build_store().prepare_for_pileup().input,
                Pos0::new(0).unwrap(),
                Pos0::new(5).unwrap(),
            );
            engine.set_soft_clip_overhang(overhang);
            let mut out = None;
            while let Some(col) = engine.pileups() {
                if col.pos() == Pos0::new(2).unwrap() {
                    out = Some(
                        PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut Vec::new())
                            .unwrap(),
                    );
                }
            }
            out
        };

        // Overhang off: nothing covers position 2, so no column is produced.
        assert!(metrics_at(0).is_none(), "no rescue without overhang");

        // Overhang on: the clipped T is rescued as an OT alt at the CpG-C.
        let pm = metrics_at(1).expect("CpG-C column emitted via soft-clip rescue");
        assert_eq!(pm.reference_base, Base::C);
        let t = pm.alt(Base::T).expect("rescued T alt present at CpG-C");
        assert_eq!(t.strand_count.ot, 1, "rescued methylation read counted on OT");
        assert_eq!(t.strand_count.ob, 0);
    }

    #[test]
    fn rescues_ob_strand_cpg_partner() {
        // Reference: T T C G T T — CpG is C@2 / G@3.
        let seg = segment(b"TTCGTT");
        let params = PileupMappingParams::default();
        let mut extras = TestExtras(params.read_masking.clone());

        let mut build_store = || {
            let mut store = RecordStore::<RastairReadExtras>::new();
            // 3M 1S at pos 0: aligned T,T,C over ref 0,1,2; clip base A over ref
            // G@3. flag 83 = paired/proper/reverse/first → OB.
            store
                .push_fields(
                    Pos0::new(0).unwrap(),
                    Pos0::new(2).unwrap(),
                    BamFlags::from(83u16),
                    60,
                    3,
                    0,
                    b"clipped",
                    &[CigarOp::new(CigarOpType::Match, 3), CigarOp::new(CigarOpType::SoftClip, 1)],
                    &[Base::T, Base::T, Base::C, Base::A],
                    &[40u8; 4],
                    &[],
                    0,
                    -1,
                    0,
                    0,
                    &mut extras,
                )
                .unwrap();
            store
        };

        let mut metrics_at = |overhang: u32| -> Option<PileupMetrics> {
            let mut engine = PileupEngine::new(
                build_store().prepare_for_pileup().input,
                Pos0::new(0).unwrap(),
                Pos0::new(5).unwrap(),
            );
            engine.set_soft_clip_overhang(overhang);
            let mut out = None;
            while let Some(col) = engine.pileups() {
                if col.pos() == Pos0::new(3).unwrap() {
                    out = Some(
                        PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut Vec::new())
                            .unwrap(),
                    );
                }
            }
            out
        };

        assert!(metrics_at(0).is_none(), "no rescue without overhang");

        let pm = metrics_at(1).expect("CpG-G column emitted via soft-clip rescue");
        assert_eq!(pm.reference_base, Base::G);
        let a = pm.alt(Base::A).expect("rescued A alt present at CpG-G");
        assert_eq!(a.strand_count.ob, 1, "rescued methylation read counted on OB");
        assert_eq!(a.strand_count.ot, 0);
    }

    /// Strand half of the gate: an OT read whose trailing clip lands on the G of
    /// a CpG must *not* be rescued — OB methylation evidence (G→A) cannot come
    /// from an OT read, so this is a plain end-of-read mismatch.
    #[test]
    fn does_not_rescue_ot_clip_over_ref_g() {
        // Reference: T T C G T T — G@3 is the CpG-G, before_1 = C@2.
        let seg = segment(b"TTCGTT");
        let params = PileupMappingParams::default();
        let mut extras = TestExtras(params.read_masking.clone());

        let mut store = RecordStore::<RastairReadExtras>::new();
        // 3M 1S at pos 0, flag 99 → OT. Trailing clip A projects onto ref G@3.
        store
            .push_fields(
                Pos0::new(0).unwrap(),
                Pos0::new(2).unwrap(),
                BamFlags::from(99u16),
                60,
                3,
                0,
                b"clipped",
                &[CigarOp::new(CigarOpType::Match, 3), CigarOp::new(CigarOpType::SoftClip, 1)],
                &[Base::T, Base::T, Base::C, Base::A],
                &[40u8; 4],
                &[],
                0,
                -1,
                0,
                0,
                &mut extras,
            )
            .unwrap();

        let mut engine = PileupEngine::new(
            store.prepare_for_pileup().input,
            Pos0::new(0).unwrap(),
            Pos0::new(5).unwrap(),
        );
        engine.set_soft_clip_overhang(1);
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(3).unwrap() {
                let pm = PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut Vec::new())
                    .unwrap();
                assert!(pm.alt(Base::A).is_none(), "OT clip over ref G must not be rescued");
                assert_eq!(pm.pos_metrics.depth, 0);
            }
        }
    }

    /// The same clipped base over a non-CpG C (no following G) is *not* rescued:
    /// it is a plain end-of-read mismatch, not a methylation partner.
    #[test]
    fn does_not_rescue_outside_cpg_context() {
        // Reference: T T C A T T — C@2 is followed by A, so not a CpG.
        let seg = segment(b"TTCATT");
        let params = PileupMappingParams::default();
        let mut extras = TestExtras(params.read_masking.clone());

        let mut store = RecordStore::<RastairReadExtras>::new();
        store
            .push_fields(
                Pos0::new(3).unwrap(),
                Pos0::new(5).unwrap(),
                BamFlags::from(99u16),
                60,
                3,
                0,
                b"clipped",
                &[CigarOp::new(CigarOpType::SoftClip, 1), CigarOp::new(CigarOpType::Match, 3)],
                &[Base::T, Base::A, Base::T, Base::T],
                &[40u8; 4],
                &[],
                0,
                -1,
                0,
                0,
                &mut extras,
            )
            .unwrap();

        let mut engine = PileupEngine::new(
            store.prepare_for_pileup().input,
            Pos0::new(0).unwrap(),
            Pos0::new(5).unwrap(),
        );
        engine.set_soft_clip_overhang(1);
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(2).unwrap() {
                let pm = PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut Vec::new())
                    .unwrap();
                // The soft-clip view exists but is gated out: no T alt, no depth.
                assert!(pm.alt(Base::T).is_none(), "non-CpG clip must not be rescued");
                assert_eq!(pm.pos_metrics.depth, 0);
            }
        }
    }

    /// A clipped base sitting on the right strand at a CpG partner but whose
    /// *observed* base is not bisulfite-relevant (here a C→G mismatch over the
    /// CpG-C on OT) must not be rescued: it is a fringe SNP/error, not
    /// methylation evidence, and rescuing it would feed noise into variant
    /// calling.
    #[test]
    fn does_not_rescue_non_bisulfite_clip_base() {
        // Reference: T T C G T T — CpG is C@2 / G@3.
        let seg = segment(b"TTCGTT");
        let params = PileupMappingParams::default();
        let mut extras = TestExtras(params.read_masking.clone());

        let mut store = RecordStore::<RastairReadExtras>::new();
        // 1S 3M at pos 3, flag 99 → OT. Clip base G (not T/C) projects onto C@2.
        store
            .push_fields(
                Pos0::new(3).unwrap(),
                Pos0::new(5).unwrap(),
                BamFlags::from(99u16),
                60,
                3,
                0,
                b"clipped",
                &[CigarOp::new(CigarOpType::SoftClip, 1), CigarOp::new(CigarOpType::Match, 3)],
                &[Base::G, Base::G, Base::T, Base::T],
                &[40u8; 4],
                &[],
                0,
                -1,
                0,
                0,
                &mut extras,
            )
            .unwrap();

        let mut engine = PileupEngine::new(
            store.prepare_for_pileup().input,
            Pos0::new(0).unwrap(),
            Pos0::new(5).unwrap(),
        );
        engine.set_soft_clip_overhang(1);
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(2).unwrap() {
                let pm = PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut Vec::new())
                    .unwrap();
                assert!(pm.alt(Base::G).is_none(), "non-bisulfite fringe clip must not be rescued");
                assert_eq!(pm.pos_metrics.depth, 0);
            }
        }
    }

    /// A paired read whose two mates both land on the same CpG-C — one via an
    /// aligned base, the other via a rescued soft-clip partner — must be counted
    /// once. The engine presents both (depth 2); `from_seqair`'s overlapping-pair
    /// dedup collapses them to a single observation, so the rescued fringe base
    /// does not double-count the molecule.
    #[test]
    fn rescued_partner_is_deduped_against_mate() {
        // Reference: T T C G T T — CpG is C@2 / G@3.
        let seg = segment(b"TTCGTT");
        let params = PileupMappingParams::default();
        let mut extras = TestExtras(params.read_masking.clone());

        let mut store = RecordStore::<RastairReadExtras>::new();
        // Mate A: first in template, forward (flag 99 → OT). 3M at pos 2 covers
        // the CpG-C with an aligned C.
        store
            .push_fields(
                Pos0::new(2).unwrap(),
                Pos0::new(4).unwrap(),
                BamFlags::from(99u16),
                60,
                3,
                0,
                b"pair",
                &[CigarOp::new(CigarOpType::Match, 3)],
                &[Base::C, Base::G, Base::T],
                &[40u8; 3],
                &[],
                0,
                0,
                3,
                0,
                &mut extras,
            )
            .unwrap();
        // Mate B: second in template, reverse (flag 147 → OT). 1S 3M at pos 3,
        // its clipped T projecting back onto the same CpG-C.
        store
            .push_fields(
                Pos0::new(3).unwrap(),
                Pos0::new(5).unwrap(),
                BamFlags::from(147u16),
                60,
                3,
                0,
                b"pair",
                &[CigarOp::new(CigarOpType::SoftClip, 1), CigarOp::new(CigarOpType::Match, 3)],
                &[Base::T, Base::G, Base::T, Base::T],
                &[40u8; 4],
                &[],
                0,
                0,
                2,
                0,
                &mut extras,
            )
            .unwrap();

        // `Readers::pileup` links mates after fetching; this test drives the
        // engine directly, so it links by hand.
        let stats = store.link_mates();
        assert_eq!(stats.pairs, 1, "the fixture's mates must link");
        let mut engine = PileupEngine::new(
            store.prepare_for_pileup().input,
            Pos0::new(0).unwrap(),
            Pos0::new(5).unwrap(),
        );
        engine.set_soft_clip_overhang(1);

        let mut checked = false;
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(2).unwrap() {
                // The engine presents both the aligned mate and the rescued clip.
                assert_eq!(col.depth(), 2, "both mates present at the CpG-C before dedup");
                let pm = PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut Vec::new())
                    .unwrap();
                assert_eq!(pm.pos_metrics.depth, 1, "rescued partner deduped against its mate");
                checked = true;
            }
        }
        assert!(checked, "CpG-C column must be produced");
    }

    /// Read-end masking must be applied consistently across both `from_seqair`
    /// passes. A rescued soft-clip CpG partner bypasses masking during counting
    /// (it is a fringe base by construction); if the *dedup* pass still applies
    /// masking to it, the rescued view is dropped from the overlapping-pair
    /// collector while its aligned mate is not — so the pair is never detected
    /// and the molecule is counted twice. With masking that targets only the
    /// clipped mate, the deduped depth must still be 1, not 2.
    #[test]
    fn rescued_partner_masking_is_consistent_across_passes() {
        // Reference: T T C G T T — CpG is C@2 / G@3.
        let seg = segment(b"TTCGTT");
        // Mask one base from the 3' end of OT reverse reads (the clipped mate B);
        // OT forward (mate A) is left untouched.
        let mut params = PileupMappingParams::default();
        params.variant_calling.read_masking =
            ReadMaskParams::new("0,0,0,1".parse().unwrap(), ReadMaskSetting::default());
        let mut extras = TestExtras(params.read_masking.clone());

        let mut store = RecordStore::<RastairReadExtras>::new();
        // Mate A: first in template, forward (flag 99 → OT). 3M at pos 2 covers
        // the CpG-C with an aligned C; OT-forward masking is zero so it survives.
        store
            .push_fields(
                Pos0::new(2).unwrap(),
                Pos0::new(4).unwrap(),
                BamFlags::from(99u16),
                60,
                3,
                0,
                b"pair",
                &[CigarOp::new(CigarOpType::Match, 3)],
                &[Base::C, Base::G, Base::T],
                &[40u8; 3],
                &[],
                0,
                0,
                3,
                0,
                &mut extras,
            )
            .unwrap();
        // Mate B: second in template, reverse (flag 147 → OT). 1S 3M at pos 3,
        // its clipped T (read pos 0) projecting onto the same CpG-C. The OT
        // reverse mask rejects read position 0.
        store
            .push_fields(
                Pos0::new(3).unwrap(),
                Pos0::new(5).unwrap(),
                BamFlags::from(147u16),
                60,
                3,
                0,
                b"pair",
                &[CigarOp::new(CigarOpType::SoftClip, 1), CigarOp::new(CigarOpType::Match, 3)],
                &[Base::T, Base::G, Base::T, Base::T],
                &[40u8; 4],
                &[],
                0,
                0,
                2,
                0,
                &mut extras,
            )
            .unwrap();

        // `Readers::pileup` links mates after fetching; this test drives the
        // engine directly, so it links by hand.
        let stats = store.link_mates();
        assert_eq!(stats.pairs, 1, "the fixture's mates must link");
        let mut engine = PileupEngine::new(
            store.prepare_for_pileup().input,
            Pos0::new(0).unwrap(),
            Pos0::new(5).unwrap(),
        );
        engine.set_soft_clip_overhang(1);

        let mut checked = false;
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(2).unwrap() {
                let pm = PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut Vec::new())
                    .unwrap();
                assert_eq!(
                    pm.pos_metrics.depth, 1,
                    "rescued partner deduped against its mate even with read-end masking active"
                );
                checked = true;
            }
        }
        assert!(checked, "CpG-C column must be produced");
    }

    // ── differential: mate links vs. the old name-collector rule ────────────

    const OT_FIRST: u16 = 99; // paired, proper, mate reverse, first in template
    const OT_SECOND: u16 = 147; // paired, proper, reverse, second in template
    const OB_FIRST: u16 = 83; // paired, proper, reverse, first in template
    const OB_SECOND: u16 = 163; // paired, proper, mate reverse, second in template

    /// `ACGTACGTACGTACGTACGT`, long enough for two 8bp reads to overlap.
    const REF: &[u8] = b"ACGTACGTACGTACGTACGT";

    fn dedup_params() -> PileupMappingParams {
        PileupMappingParams::default()
    }

    /// Mates that agree on every base: the later one goes, whichever way round
    /// the pair is.
    #[test]
    fn dedup_matches_the_old_rule_when_mates_agree() {
        for (a_flags, b_flags) in [(OT_FIRST, OT_SECOND), (OB_SECOND, OB_FIRST)] {
            let reads = vec![
                TestRead::matching(b"pair", 2, 8, Base::A, a_flags),
                TestRead::matching(b"pair", 6, 8, Base::A, b_flags),
            ];
            assert_same_as_name_collector(&reads, REF, &dedup_params());
        }
    }

    /// Mates that disagree at the overlapped base. Which one survives depends
    /// on which is second in the template — the case the single-pass rewrite
    /// had to get right without seeing the pair as a unit.
    #[test]
    fn dedup_matches_the_old_rule_when_mates_disagree() {
        for (a_flags, b_flags) in [(OT_FIRST, OT_SECOND), (OT_SECOND, OT_FIRST)] {
            let reads = vec![
                TestRead::matching(b"pair", 2, 8, Base::A, a_flags),
                TestRead::matching(b"pair", 6, 8, Base::A, b_flags).with_base_at(0, Base::G),
            ];
            assert_same_as_name_collector(&reads, REF, &dedup_params());
        }
    }

    /// A mate that fails base quality, mapping quality, or read-end masking at
    /// the shared column never forms a pair there, so the other one survives
    /// even though both cover the position.
    #[test]
    fn dedup_matches_the_old_rule_when_one_mate_is_filtered() {
        let params = dedup_params();
        let low_baseq = vec![
            TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST),
            TestRead::matching(b"pair", 6, 8, Base::A, OT_SECOND).with_qual(2),
        ];
        assert_same_as_name_collector(&low_baseq, REF, &params);

        let low_mapq = vec![
            TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST).with_mapq(0),
            TestRead::matching(b"pair", 6, 8, Base::A, OT_SECOND),
        ];
        assert_same_as_name_collector(&low_mapq, REF, &params);

        // Mask the first two bases of every read: at the overlap's left edge
        // one mate is masked out while the other is not.
        let mut masked = dedup_params();
        masked.variant_calling.read_masking =
            ReadMaskParams::new("2,0,2,0".parse().unwrap(), "2,0,2,0".parse().unwrap());
        let reads = vec![
            TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST),
            TestRead::matching(b"pair", 6, 8, Base::A, OT_SECOND),
        ];
        assert_same_as_name_collector(&reads, REF, &masked);
    }

    /// A deletion in one mate: at the deleted positions it has no base, so no
    /// pair forms and the other mate stands alone.
    #[test]
    fn dedup_matches_the_old_rule_across_a_deletion() {
        let mut with_del = TestRead::matching(b"pair", 6, 8, Base::A, OT_SECOND);
        with_del.cigar = vec![
            CigarOp::new(CigarOpType::Match, 2),
            CigarOp::new(CigarOpType::Deletion, 2),
            CigarOp::new(CigarOpType::Match, 6),
        ];
        with_del.bases = vec![Base::A; 8];
        with_del.quals = vec![40; 8];
        let reads = vec![TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST), with_del];
        assert_same_as_name_collector(&reads, REF, &dedup_params());
    }

    /// An insertion in one mate — the anchor base is still a normal
    /// observation, and the pair must resolve there like any other column.
    #[test]
    fn dedup_matches_the_old_rule_across_an_insertion() {
        let mut with_ins = TestRead::matching(b"pair", 6, 8, Base::A, OT_SECOND);
        with_ins.cigar = vec![
            CigarOp::new(CigarOpType::Match, 2),
            CigarOp::new(CigarOpType::Insertion, 2),
            CigarOp::new(CigarOpType::Match, 6),
        ];
        let reads = vec![TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST), with_ins];
        assert_same_as_name_collector(&reads, REF, &dedup_params());
    }

    /// The exact edges of the overlap: the pair meets on one base only, at the
    /// last position of the left mate. This is where the half-open/inclusive
    /// mix-up in the overlap interval showed up on real data.
    #[test]
    fn dedup_matches_the_old_rule_at_the_overlap_boundaries() {
        // Left mate covers 2..=9, right mate 9..=16: they share exactly base 9.
        let reads = vec![
            TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST),
            TestRead::matching(b"pair", 9, 8, Base::A, OT_SECOND),
        ];
        assert_same_as_name_collector(&reads, REF, &dedup_params());

        // And one base further apart: no shared position at all.
        let disjoint = vec![
            TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST),
            TestRead::matching(b"pair", 10, 8, Base::A, OT_SECOND),
        ];
        assert_same_as_name_collector(&disjoint, REF, &dedup_params());
    }

    /// Several pairs at one column, interleaved, some agreeing and some not —
    /// the pending-drop bookkeeping has to keep them apart.
    #[test]
    fn dedup_matches_the_old_rule_for_interleaved_pairs() {
        let reads = vec![
            TestRead::matching(b"p1", 2, 8, Base::A, OT_FIRST),
            TestRead::matching(b"p2", 3, 8, Base::A, OB_FIRST),
            TestRead::matching(b"p3", 4, 8, Base::A, OT_SECOND),
            TestRead::matching(b"p1", 6, 8, Base::A, OT_SECOND).with_base_at(0, Base::G),
            TestRead::matching(b"p2", 7, 8, Base::A, OB_SECOND),
            TestRead::matching(b"p3", 8, 8, Base::A, OT_FIRST).with_base_at(0, Base::T),
        ];
        assert_same_as_name_collector(&reads, REF, &dedup_params());
    }

    /// With `--rescue-soft-clip-cpg`, a clipped base is projected onto a column
    /// outside its own alignment. It still belongs to the same molecule as its
    /// mate's aligned base there, so it must dedup against it — which is why
    /// the mate-overlap interval is widened by the overhang.
    #[test]
    fn dedup_matches_the_old_rule_for_a_rescued_soft_clip() {
        // REF has a CpG at 5 (C) / 6 (G). The left mate covers it directly; the
        // right mate's clipped T projects back onto the C.
        let mut clipped = TestRead::matching(b"pair", 6, 3, Base::G, OT_SECOND);
        clipped.cigar =
            vec![CigarOp::new(CigarOpType::SoftClip, 1), CigarOp::new(CigarOpType::Match, 3)];
        clipped.bases = vec![Base::T, Base::G, Base::T, Base::A];
        clipped.quals = vec![40; 4];

        let reads = vec![TestRead::matching(b"pair", 5, 3, Base::C, OT_FIRST), clipped];
        let params = PileupMappingParams { rescue_soft_clip_cpg: true, ..dedup_params() };
        assert_same_as_name_collector(&reads, REF, &params);

        // And the column is genuinely one where dedup has to fire: without the
        // widened interval both mates would be counted at the CpG-C.
        let seg = segment(REF);
        let store = store_of(&reads, &params.read_masking);
        let mut engine = PileupEngine::new(
            store.prepare_for_pileup().input,
            Pos0::new(0).unwrap(),
            Pos0::new(19).unwrap(),
        );
        engine.set_soft_clip_overhang(1);
        let mut scratch = Vec::new();
        let mut depth_at_cpg = None;
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(5).unwrap() {
                assert_eq!(col.depth(), 2, "engine presents both the aligned base and the clip");
                let pm =
                    PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut scratch).unwrap();
                depth_at_cpg = Some(pm.pos_metrics.depth);
            }
        }
        assert_eq!(depth_at_cpg, Some(1), "the rescued clip deduped against its mate");
    }

    /// `--max-coverage` truncates the column before either implementation sees
    /// it, so both must agree on the truncated column too.
    #[test]
    fn dedup_matches_the_old_rule_under_max_depth_truncation() {
        let mut reads = Vec::new();
        for i in 0..6u32 {
            let name = format!("p{i}");
            reads.push(TestRead::matching(name.as_bytes(), 2, 8, Base::A, OT_FIRST));
            reads.push(TestRead::matching(name.as_bytes(), 6, 8, Base::A, OT_SECOND));
        }
        let mut capped = dedup_params();
        capped.variant_calling.max_coverage = crate::call::variant_calling::MaxCoverage::new(3);
        assert_same_as_name_collector(&reads, REF, &capped);
    }

    /// With `--keep-overlapping-reads` no dedup happens at all, so every
    /// filtered observation survives — including both halves of a pair.
    #[test]
    fn keeping_overlapping_reads_keeps_both_mates() {
        let reads = vec![
            TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST),
            TestRead::matching(b"pair", 6, 8, Base::A, OT_SECOND),
        ];
        let mut params = dedup_params();
        params.variant_calling.keep_overlapping_reads = true;
        let seg = segment(REF);
        let store = store_of(&reads, &params.read_masking);
        let mut engine = PileupEngine::new(
            store.prepare_for_pileup().input,
            Pos0::new(0).unwrap(),
            Pos0::new(19).unwrap(),
        );
        let mut scratch = Vec::new();
        let mut overlap_depth = None;
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(7).unwrap() {
                let pm =
                    PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut scratch).unwrap();
                overlap_depth = Some(pm.pos_metrics.depth);
            }
        }
        assert_eq!(overlap_depth, Some(2), "both mates must survive inside the overlap");
    }
}
