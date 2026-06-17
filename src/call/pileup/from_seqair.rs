#![cfg(feature = "experimental-seqair")]

use super::{
    indels::{IndelAllele, IndelObservation},
    overlapping_reads::{DedupInfo, NameCollector, resolve_pair},
    ref_features::{dinucleotide_run_at, homopolymer_run_at, indel_ref_window_at},
};
use crate::{
    call::{
        pileup::{Pileup, PositionInRead, SimpleRead, SimpleReads},
        process::PileupMappingParams,
    },
    sequence::Segment,
    utils::SequenceContext,
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use seqair::bam::pileup::{AlignmentView, Indel, PileupColumn};
use seqair_types::{Base, SmallVec};
use std::rc::Rc;
use tracing::{debug, instrument, trace};

use crate::sequence::RastairReadExtras;

impl Pileup {
    #[instrument(level = "trace", skip_all)]
    pub(crate) fn from_seqair(
        column: &PileupColumn<'_, RastairReadExtras>,
        segment: Rc<Segment>,
        params: &PileupMappingParams,
        collector: &mut NameCollector,
    ) -> Result<Pileup> {
        let pos = column.pos().as_u64();
        let pos_u32 = u32::try_from(pos).wrap_err("pileup position exceeds u32")?;
        let idx = segment.pos_to_idx(pos_u32)?;
        let depth = column.depth();
        let max_reads = depth.min(params.max_coverage as usize);
        if depth > max_reads {
            debug!(pos, depth, "Capping number of reads in pileup to {max_reads}");
        }

        // First pass: determine which overlapping-pair indices to remove.
        // The alignments slice is already materialized by seqair, so walking
        // it twice is cheap and lets the second pass aggregate indel stats
        // directly without a per-read buffer.
        let mut to_remove = SmallVec::<usize, 16>::new();
        collector.prepare(max_reads);
        if matches!(collector, NameCollector::Collect(..)) {
            let filtered = column
                .alignments()
                .enumerate()
                .filter_map(|(idx, view)| {
                    let read = extract_read(&view, params)?;
                    Some((idx, view.qname(), read))
                })
                .take(max_reads);

            for (idx, name, read) in filtered {
                let info = DedupInfo { idx, base: read.base, second: read.second };
                if let Some(other) = collector.see(name, info) {
                    resolve_pair(&info, info.idx, &other, other.idx, &mut to_remove);
                }
            }
            to_remove.sort_unstable();
        }

        let mut raw_reads: Vec<SimpleRead> = Vec::with_capacity(max_reads);
        let mut indel_observations = SmallVec::new();
        let mut depth_offset: u32 = 0;
        let mut soft_clip_count: u32 = 0;

        let iter = column
            .alignments()
            .enumerate()
            .filter(|(idx, _)| to_remove.binary_search(idx).is_err())
            .filter_map(|(_, view)| {
                let read = extract_read(&view, params)?;
                Some((read, view))
            })
            .take(max_reads);

        for (read, view) in iter {
            raw_reads.push(read);

            if params.call_indels {
                let aln = view.alignment();
                let extras = view.extra();

                if extras.has_soft_clip {
                    soft_clip_count += 1;
                }
                if extras.has_repeat && matches!(aln.indel_after(), Indel::None) {
                    depth_offset += 1;
                }

                if let Some(obs) = build_indel_observation(&view, pos, segment.as_ref(), params) {
                    indel_observations.push(obs);
                }
            }
        }

        let reads = SimpleReads(raw_reads.into());

        let reference_base: Base =
            segment.sequence.get(idx).wrap_err("failed to get reference base")?.into();

        let context =
            SequenceContext::new(idx, &segment).wrap_err("failed to get sequence context")?;

        let (indel_ref_window, indel_ref_anchor) = if indel_observations.is_empty() {
            (SmallVec::new(), 0)
        } else {
            indel_ref_window_at(idx, &segment)
        };

        let segment_start = segment.range.region.start as usize;

        Ok(Pileup {
            region: segment.range.clone(),
            context,
            pos: pos_u32,
            reads,
            reference_base,
            indel_observations,
            depth_offset,
            homopolymer_run: homopolymer_run_at(pos as usize, &segment, segment_start),
            dinucleotide_run: dinucleotide_run_at(pos as usize, &segment, segment_start),
            soft_clip_count,
            indel_ref_window,
            indel_ref_anchor,
        })
    }
}

fn extract_read(
    view: &AlignmentView<'_, '_, RastairReadExtras>,
    params: &PileupMappingParams,
) -> Option<SimpleRead> {
    let aln = view.alignment();
    let extras = view.extra();

    let qpos = aln.qpos()?;
    let base = aln.base()?;
    let qual_bq = aln.qual()?;
    let qual = qual_bq.get().unwrap_or(0);

    let read = SimpleRead {
        base,
        qual,
        mapq: aln.mapq,
        strand: extras.strand,
        reverse: aln.flags.is_reverse(),
        second: aln.flags.is_second_in_template(),
        position: PositionInRead { pos: u32::try_from(qpos).ok()?, read_length: aln.seq_len },
        matching_bases: aln.matching_bases,
        indels: aln.indel_bases,
    };

    if !params.read_masking.filter(&read) || !params.quality.filter(&read) {
        return None;
    }

    Some(read)
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
