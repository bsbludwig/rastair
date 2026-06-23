#![cfg(feature = "experimental-seqair")]

use super::{
    indels::{IndelAllele, IndelObservation},
    overlapping_reads::{DedupInfo, NameCollector, resolve_pair},
    ref_features::{dinucleotide_run_at, homopolymer_run_at, indel_ref_window_at},
};
use crate::{
    call::{process::PileupMappingParams, variant_calling::ReadMaskParams},
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
    pub(crate) fn from_seqair(
        column: &PileupColumn<'_, RastairReadExtras>,
        segment: Rc<Segment>,
        params: &PileupMappingParams,
        collector: &mut NameCollector,
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

        // First pass: determine which overlapping-pair indices to remove.
        let mut to_remove = SmallVec::<usize, 16>::new();
        collector.prepare(max_reads);
        if matches!(collector, NameCollector::Collect(..)) {
            let filtered = column
                .alignments()
                .enumerate()
                .filter_map(|(idx, view)| {
                    let mapq = view.mapq;
                    let baseq = view.qual()?;

                    if !passes_read_masking(&view, reference_base, &context, &params.read_masking) {
                        return None;
                    }
                    if !params.quality.filter_fields(mapq, baseq.get()?) {
                        return None;
                    }

                    let info = DedupInfo {
                        idx,
                        base: view.base()?,
                        second: view.flags.is_second_in_template(),
                    };
                    Some((view.qname(), info))
                })
                .take(max_reads);

            for (name, info) in filtered {
                if let Some(other) = collector.see(name, info) {
                    resolve_pair(&info, info.idx, &other, other.idx, &mut to_remove);
                }
            }
            to_remove.sort_unstable();
        }

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

        for (_idx, view) in column
            .alignments()
            .enumerate()
            .filter(|(idx, _)| to_remove.binary_search(idx).is_err())
            .take(max_reads)
        {
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
            if !passes_read_masking(&view, reference_base, &context, &params.read_masking) {
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

fn passes_read_masking(
    view: &AlignmentView<'_, '_, RastairReadExtras>,
    reference_base: Base,
    context: &SequenceContext,
    read_masking: &ReadMaskParams,
) -> bool {
    let strand = view.extra().strand;
    if view.is_soft_clip() {
        let Some(observed) = view.base() else { return false };
        soft_clip_cpg_partner(reference_base, observed, context, strand)
    } else {
        let Some(pos) = view.qpos() else { return false };
        read_masking.filter_fields(strand, view.flags.is_reverse(), pos as u32, view.seq_len)
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

    #[derive(Default, Clone)]
    struct TestExtras;

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
            RastairReadExtras {
                strand: Strand::from(rec.flags),
                has_soft_clip,
                has_repeat: false,
                taps_aware_mismatches: 0,
            }
        }
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

        let build_store = || {
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
                    &mut TestExtras,
                )
                .unwrap();
            store
        };

        let metrics_at = |overhang: u32| -> Option<PileupMetrics> {
            let mut engine =
                PileupEngine::new(build_store(), Pos0::new(0).unwrap(), Pos0::new(5).unwrap());
            engine.set_soft_clip_overhang(overhang);
            let mut collector = NameCollector::new(&params);
            let mut out = None;
            while let Some(col) = engine.pileups() {
                if col.pos() == Pos0::new(2).unwrap() {
                    out = Some(
                        PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut collector)
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

        let build_store = || {
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
                    &mut TestExtras,
                )
                .unwrap();
            store
        };

        let metrics_at = |overhang: u32| -> Option<PileupMetrics> {
            let mut engine =
                PileupEngine::new(build_store(), Pos0::new(0).unwrap(), Pos0::new(5).unwrap());
            engine.set_soft_clip_overhang(overhang);
            let mut collector = NameCollector::new(&params);
            let mut out = None;
            while let Some(col) = engine.pileups() {
                if col.pos() == Pos0::new(3).unwrap() {
                    out = Some(
                        PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut collector)
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
                &mut TestExtras,
            )
            .unwrap();

        let mut engine = PileupEngine::new(store, Pos0::new(0).unwrap(), Pos0::new(5).unwrap());
        engine.set_soft_clip_overhang(1);
        let mut collector = NameCollector::new(&params);
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(3).unwrap() {
                let pm =
                    PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut collector).unwrap();
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
                &mut TestExtras,
            )
            .unwrap();

        let mut engine = PileupEngine::new(store, Pos0::new(0).unwrap(), Pos0::new(5).unwrap());
        engine.set_soft_clip_overhang(1);
        let mut collector = NameCollector::new(&params);
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(2).unwrap() {
                let pm =
                    PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut collector).unwrap();
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
                &mut TestExtras,
            )
            .unwrap();

        let mut engine = PileupEngine::new(store, Pos0::new(0).unwrap(), Pos0::new(5).unwrap());
        engine.set_soft_clip_overhang(1);
        let mut collector = NameCollector::new(&params);
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(2).unwrap() {
                let pm =
                    PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut collector).unwrap();
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
                -1,
                0,
                0,
                &mut TestExtras,
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
                -1,
                0,
                0,
                &mut TestExtras,
            )
            .unwrap();

        let mut engine = PileupEngine::new(store, Pos0::new(0).unwrap(), Pos0::new(5).unwrap());
        engine.set_soft_clip_overhang(1);
        let mut collector = NameCollector::new(&params);

        let mut checked = false;
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(2).unwrap() {
                // The engine presents both the aligned mate and the rescued clip.
                assert_eq!(col.depth(), 2, "both mates present at the CpG-C before dedup");
                let pm =
                    PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut collector).unwrap();
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
                -1,
                0,
                0,
                &mut TestExtras,
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
                -1,
                0,
                0,
                &mut TestExtras,
            )
            .unwrap();

        let mut engine = PileupEngine::new(store, Pos0::new(0).unwrap(), Pos0::new(5).unwrap());
        engine.set_soft_clip_overhang(1);
        let mut collector = NameCollector::new(&params);

        let mut checked = false;
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(2).unwrap() {
                let pm =
                    PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut collector).unwrap();
                assert_eq!(
                    pm.pos_metrics.depth, 1,
                    "rescued partner deduped against its mate even with read-end masking active"
                );
                checked = true;
            }
        }
        assert!(checked, "CpG-C column must be produced");
    }
}
