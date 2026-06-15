#![cfg(feature = "experimental-seqair")]

use super::{
    indels::{IndelAllele, IndelObservation},
    overlapping_reads::{NameCollector, resolve_pair},
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
use seqair::bam::pileup::{Indel, PileupColumn};
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

        let mut raw_reads: Vec<SimpleRead> = Vec::with_capacity(max_reads);

        let alignments_with_names = column
            .alignments()
            .filter_map(|view| {
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
                    position: PositionInRead {
                        pos: u32::try_from(qpos).ok()?,
                        read_length: aln.seq_len,
                    },
                    matching_bases: aln.matching_bases,
                    indels: aln.indel_bases,
                };

                Some((view.qname(), read))
            })
            .filter(|(_, read)| params.read_masking.filter(read))
            .filter(|(_, read)| params.quality.filter(read))
            .take(max_reads);

        let reads = match collector {
            NameCollector::Skip => {
                for (_, read) in alignments_with_names {
                    raw_reads.push(read);
                }
                SimpleReads(raw_reads.into())
            }
            NameCollector::Collect(buf) => {
                buf.prepare(max_reads);
                let mut to_remove = SmallVec::<usize, 16>::new();
                for (name, read) in alignments_with_names {
                    let this_idx = raw_reads.len();
                    raw_reads.push(read);
                    if let Some(other_idx) = buf.see(name, this_idx) {
                        resolve_pair(&raw_reads, this_idx, other_idx, &mut to_remove);
                    }
                }
                to_remove.sort_unstable();
                for &idx in to_remove.iter().rev() {
                    raw_reads.swap_remove(idx);
                }
                SimpleReads(raw_reads.into())
            }
        };

        let reference_base: Base =
            segment.sequence.get(idx).wrap_err("failed to get reference base")?.into();

        let context =
            SequenceContext::new(idx, &segment).wrap_err("failed to get sequence context")?;

        // Second pass: collect indels and compute depth_offset / soft_clip_count.
        let indel_cutoff = params.indel_end_of_read_cutoff;
        let mut indel_observations = SmallVec::new();
        let mut depth_offset: u32 = 0;
        let mut soft_clip_count: u32 = 0;
        let segment_start = segment.range.region.start as usize;

        for view in column.alignments() {
            let aln = view.alignment();
            let extras = view.extra();

            if extras.has_soft_clip {
                soft_clip_count += 1;
            }

            // Anchor-addressed indel, htslib-`a.indel()` style. seqair reports
            // it on the anchor (the Match column before a deletion, or the
            // Insertion column), so this mirrors `from_hts`'s `a.indel()` arm
            // directly — no per-read deletion-anchor bookkeeping needed.
            match aln.indel_after() {
                Indel::None => {
                    if extras.has_repeat {
                        depth_offset += 1;
                    }
                }
                indel => {
                    // `qpos` is the anchor's query position (Some for Match and
                    // Insertion ops, which are the only ones with a non-`None`
                    // `indel_after`).
                    let Some(qpos) = aln.qpos() else { continue };
                    let read_len = aln.seq_len as usize;

                    if qpos < indel_cutoff || qpos >= read_len.saturating_sub(indel_cutoff) {
                        trace!(qpos, read_len, "Indel skipped: too close to read end");
                        continue;
                    }
                    if extras.taps_aware_mismatches > params.indel_max_mismatches {
                        trace!(
                            mismatches = extras.taps_aware_mismatches,
                            max = params.indel_max_mismatches,
                            "Indel skipped: too many non-TAPS mismatches"
                        );
                        continue;
                    }

                    let (allele, insertion_base_quals, post_del_base_qual) = match indel {
                        Indel::Insertion(_) => {
                            let bases: SmallVec<Base, 4> =
                                view.inserted_bases().iter().copied().collect();
                            if bases.is_empty() {
                                continue;
                            }
                            let quals: SmallVec<u8, 4> =
                                view.inserted_quals().iter().filter_map(|q| q.get()).collect();
                            (IndelAllele::Insertion(bases), quals, 0)
                        }
                        Indel::Deletion(del_len) => {
                            // Deleted bases come from the reference:
                            // `seq[pos+1 .. pos+1+del_len]` (`pos` is the anchor).
                            let ref_start = (pos as usize + 1).saturating_sub(segment_start);
                            let ref_end = ref_start + del_len as usize;
                            let bases: SmallVec<Base, 4> = segment
                                .sequence
                                .get(ref_start..ref_end)
                                .map(|s| s.iter().map(|&b| Base::from(b)).collect())
                                .unwrap_or_default();
                            if bases.is_empty() {
                                continue;
                            }
                            let post_del =
                                view.qualities().get(qpos + 1).and_then(|q| q.get()).unwrap_or(0);
                            (IndelAllele::Deletion(bases), SmallVec::new(), post_del)
                        }
                        Indel::None => unreachable!("matched a non-None indel above"),
                    };

                    let base_qual = view.qualities().get(qpos).and_then(|q| q.get()).unwrap_or(0);

                    indel_observations.push(IndelObservation {
                        allele,
                        strand: extras.strand,
                        reverse: aln.flags.is_reverse(),
                        pos_in_read: u32::try_from(qpos).wrap_err("qpos overflow")?,
                        read_length: aln.seq_len,
                        mapq: aln.mapq,
                        base_qual,
                        matching_bases: aln.matching_bases,
                        num_indels_in_read: aln.indel_bases,
                        insertion_base_quals,
                        post_del_base_qual,
                        has_repeat: extras.has_repeat,
                    });
                }
            }
        }

        let (indel_ref_window, indel_ref_anchor) = if indel_observations.is_empty() {
            (SmallVec::new(), 0)
        } else {
            indel_ref_window_at(idx, &segment)
        };

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
