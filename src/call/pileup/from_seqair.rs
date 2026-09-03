#![cfg(feature = "experimental-seqair")]

use super::{
    indels::{IndelAllele, IndelObservation},
    ref_features::{dinucleotide_run_at, homopolymer_run_at, indel_ref_window_at},
};
use crate::{
    call::process::PileupMappingParams,
    metrics::{
        Alt, AltFilters, Filters, PairedCounts, PerBaseAccumulators, PileupMetrics, ReadKey,
        RecordTags, aggregate_indels,
    },
    sequence::{RastairReadExtras, Segment},
    utils::SequenceContext,
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use seqair::bam::pileup::{AlignmentView, Indel, PileupColumn};
use seqair_types::{Base, RmsAccumulator, SmallVec, Strand};
use std::rc::Rc;
use tracing::{debug, instrument, trace};

impl PileupMetrics {
    #[instrument(level = "trace", skip_all)]
    /// `mate_drops` is scratch: a reusable buffer for the right mates this
    /// column drops. It is cleared here, and lives across columns only so a
    /// deep pileup does not allocate one per position.
    pub(crate) fn from_seqair(
        column: &PileupColumn<'_, RastairReadExtras>,
        segment: Rc<Segment>,
        params: &PileupMappingParams,
        mate_drops: &mut Vec<u32>,
    ) -> Result<PileupMetrics> {
        let pos = column.pos().as_u64();
        let pos_u32 = u32::try_from(pos).wrap_err("pileup position exceeds u32")?;
        let idx = segment.pos_to_idx(pos_u32)?;
        let depth = column.depth();
        let max_reads = depth.min(params.max_coverage as usize);
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
            let Some(baseq) = view.qual().and_then(|q| q.get()) else {
                continue;
            };
            let Some(base) = view.base() else {
                continue;
            };
            let Some(pos) = view.qpos() else {
                continue;
            };
            let strand = view.extra().strand;
            if !params.quality.filter_fields(view.mapq, baseq) {
                continue;
            }
            if !passes_read_masking(&view, reference_base, &context) {
                continue;
            }
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
                pos as u32,
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
                if let Some(&adj) = pos.checked_sub(1).and_then(|i| seq.get(i)) {
                    before_counts.increment(ReadKey { strand, current: base, adj });
                }
                // An adjacent base is only adjacent in the alignment when nothing
                // is inserted or deleted between them.
                if matches!(view.alignment().indel_after(), Indel::None)
                    && let Some(&adj) = seq.get(pos + 1)
                {
                    after_counts.increment(ReadKey { strand, current: base, adj });
                }
            }

            if params.call_indels {
                let aln = view.alignment();
                let extras = view.extra();

                if extras.has_soft_clip {
                    soft_clip_count += 1;
                }
                if extras.has_repeat && matches!(aln.indel_after(), Indel::None) {
                    depth_offset += 1;
                }

                if let Some(obs) =
                    build_indel_observation(&view, pos as u64, segment.as_ref(), params)
                {
                    indel_observations.push(obs);
                }
            }
        }

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
            let segment_start = segment.range.region.start as usize;
            Some(Box::new(crate::call::pileup::indels::IndelData {
                observations: indel_observations,
                ref_window: indel_ref_window,
                ref_anchor: indel_ref_anchor,
                homopolymer_run: homopolymer_run_at(pos as usize, &segment, segment_start),
                dinucleotide_run: dinucleotide_run_at(pos as usize, &segment, segment_start),
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
    let Some(mate_base) = mate.base() else { return false };
    let Some(mate_baseq) = mate.qual().and_then(|q| q.get()) else { return false };
    if !params.quality.filter_fields(mate.mapq, mate_baseq)
        || !passes_read_masking(&mate, reference_base, context)
    {
        return false;
    }

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
        let Some(qpos) = u32::try_from(view.qpos().unwrap_or(usize::MAX)).ok() else {
            return false;
        };
        view.extra().mask.contains(&qpos)
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
    pos: u64,
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
    let read_len = aln.seq_len as usize;
    let indel_cutoff = params.indel_end_of_read_cutoff;

    if qpos < indel_cutoff || qpos >= read_len.saturating_sub(indel_cutoff) {
        trace!(qpos, read_len, "Indel skipped: too close to read end");
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
            let ref_start = (pos as usize + 1).saturating_sub(segment_start);
            let ref_end = ref_start + del_len as usize;
            let bases: SmallVec<Base, 4> = segment
                .sequence
                .get(ref_start..ref_end)
                .map(|s| s.iter().map(|&b| Base::from(b)).collect())
                .unwrap_or_default();
            if bases.is_empty() {
                return None;
            }
            let post_del = view.qualities().get(qpos + 1).and_then(|q| q.get()).unwrap_or(0);
            (IndelAllele::Deletion(bases), SmallVec::new(), post_del)
        }
        Indel::None => unreachable!("matched a non-None indel above"),
    };

    let base_qual = view.qualities().get(qpos).and_then(|q| q.get()).unwrap_or(0);

    Some(IndelObservation {
        allele,
        strand: extras.strand,
        reverse: aln.flags.is_reverse(),
        pos_in_read: u32::try_from(qpos).ok().unwrap_or(0),
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
        let _stats = store.link_mates();
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
                view.base()?;
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
        let mut engine = PileupEngine::new(store, Pos0::new(0).unwrap(), Pos0::new(last).unwrap());
        if params.rescue_soft_clip_cpg {
            engine.set_soft_clip_overhang(1);
        }
        engine.set_max_depth(params.max_coverage);

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
        let end = seq.len().saturating_sub(1) as u64;
        Rc::new(Segment {
            range: ChunkRegion {
                region: Region { contig: "chr1".into(), start: 0, end },
                last_position: seq.len() as u64,
                overlap_start: 0,
                overlap_end: 0,
            },
            sequence: seq.to_vec(),
            overlap_start: 0,
            overlap_end: 0,
        })
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
            let mut engine =
                PileupEngine::new(build_store(), Pos0::new(0).unwrap(), Pos0::new(5).unwrap());
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
            let mut engine =
                PileupEngine::new(build_store(), Pos0::new(0).unwrap(), Pos0::new(5).unwrap());
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

        let mut engine = PileupEngine::new(store, Pos0::new(0).unwrap(), Pos0::new(5).unwrap());
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

        let mut engine = PileupEngine::new(store, Pos0::new(0).unwrap(), Pos0::new(5).unwrap());
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

        let mut engine = PileupEngine::new(store, Pos0::new(0).unwrap(), Pos0::new(5).unwrap());
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
        store.link_mates();
        let mut engine = PileupEngine::new(store, Pos0::new(0).unwrap(), Pos0::new(5).unwrap());
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
        store.link_mates();
        let mut engine = PileupEngine::new(store, Pos0::new(0).unwrap(), Pos0::new(5).unwrap());
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
        let mut engine = PileupEngine::new(store, Pos0::new(0).unwrap(), Pos0::new(19).unwrap());
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
        capped.variant_calling.max_coverage = 3;
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
        let mut engine = PileupEngine::new(store, Pos0::new(0).unwrap(), Pos0::new(19).unwrap());
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
