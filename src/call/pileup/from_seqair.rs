#![cfg(feature = "experimental-seqair")]

use super::{
    indels::{IndelAllele, IndelObservation},
    overlapping_reads::{DedupInfo, NameCollector, resolve_pair},
    ref_features::{dinucleotide_run_at, homopolymer_run_at, indel_ref_window_at},
};
use crate::{
    call::process::PileupMappingParams,
    metrics::{
        Alt, AltFilters, Filters, PerBaseAccumulators, PileupMetrics, RecordTags, aggregate_indels,
    },
    sequence::Segment,
    utils::SequenceContext,
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use seqair::bam::pileup::{AlignmentView, Indel, PileupColumn};
use seqair_types::{Base, RmsAccumulator, SmallVec};
use std::rc::Rc;
use tracing::{debug, instrument, trace};

use crate::sequence::RastairReadExtras;

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

        // First pass: determine which overlapping-pair indices to remove.
        let mut to_remove = SmallVec::<usize, 16>::new();
        collector.prepare(max_reads);
        if matches!(collector, NameCollector::Collect(..)) {
            let filtered = column
                .alignments()
                .enumerate()
                .filter_map(|(idx, view)| {
                    let strand = view.extra().strand;
                    let reverse = view.flags.is_reverse();
                    let pos = view.qpos()?;
                    let len = view.seq_len;
                    let mapq = view.mapq;
                    let baseq = view.qual()?;

                    if !params.read_masking.filter_fields(strand, reverse, pos as u32, len)
                        || !params.quality.filter_fields(mapq, baseq.get()?)
                    {
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

        let reference_base: Base =
            segment.sequence.get(idx).wrap_err("failed to get reference base")?.into();

        let context =
            SequenceContext::new(idx, &segment).wrap_err("failed to get sequence context")?;

        let mut accumulators = PerBaseAccumulators::default();
        let mut pos_baseq = RmsAccumulator::new();
        let mut pos_mapq = RmsAccumulator::new();
        let mut mapq0: u32 = 0;
        let mut total_depth: usize = 0;
        let mut alt_bases: SmallVec<Base, 4> = SmallVec::new();
        let mut indel_observations = SmallVec::new();
        let mut depth_offset: u32 = 0;
        let mut soft_clip_count: u32 = 0;

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
            let reverse = view.flags.is_reverse();
            if !params.quality.filter_fields(view.mapq, baseq)
                || !params.read_masking.filter_fields(strand, reverse, pos as u32, view.seq_len)
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

        let indels = aggregate_indels(&indel_observations, total_depth, depth_offset, pos_u32);

        let (indel_ref_window, indel_ref_anchor) = if indel_observations.is_empty() {
            (SmallVec::new(), 0)
        } else {
            indel_ref_window_at(idx, &segment)
        };

        let segment_start = segment.range.region.start as usize;

        Ok(PileupMetrics {
            region: segment.range.clone(),
            pos: pos_u32,
            reference_base,
            context,
            indel_observations,
            homopolymer_run: homopolymer_run_at(pos as usize, &segment, segment_start),
            dinucleotide_run: dinucleotide_run_at(pos as usize, &segment, segment_start),
            soft_clip_count,
            indel_ref_window,
            indel_ref_anchor,
            pos_metrics,
            pos_filters: Filters::default(),
            ref_metrics,
            alts,
            tags: RecordTags::default(),
            indels,
            indel_calls: Vec::new(),
        })
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
    })
}
