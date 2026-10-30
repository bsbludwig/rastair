//! A cluster scored by the pair-HMM, and its scores handed to its columns.
//!
//! [`score_cluster`] scores a cluster once through the sub-segment's
//! [`ClusterContext`], over the reads `records_overlapping(first, last)`
//! returns; the engine's store stays valid until it is reclaimed. Every column
//! of a cluster shares its `Arc<HaplotypeScores>`.
//!
//! - [`WindowMask`]: positions of the sub-segment's kept columns with at least
//!   `MIN_READS` reads of a base TAPS cannot explain are `Base::Unknown` in the
//!   haplotype window, which compair matches to any base. A `T` over a `C`
//!   counts only on OB and an `A` over a `G` only on OT, so methylation never
//!   masks. Without it a nearby SNV absent from every haplotype is charged to
//!   the reads carrying it, and a stray indel allele that absorbs the mismatch
//!   gets the credit.
//! - The reads are not the column's counted ones: one per fragment (the mate
//!   reaching further), MAPQ at least the column minimum, a known strand, soft
//!   clips included where the segment's reference reaches them, from the
//!   store seqair loaded under `--max-coverage`, with no further thinning.

use super::{ClusterPass, ColumnScratch, gpu};
use crate::call::pileup::indels::Scored;
use crate::{
    call::{
        pair_hmm::{self, Event, HaplotypeRead, HaplotypeScores, MAX_EVENTS, Prepared},
        phasing::FragmentId,
        process::PileupMappingParams,
    },
    metrics::PileupMetrics,
    sequence::{RastairReadExtras, Segment},
    vcf::CpgSide,
};
use seqair::bam::{RecordIdx, RecordRef, record_store::RecordStore};
use seqair_types::{Base, Offset, Pos0, SmallVec, Strand};
use std::{cmp::Reverse, ops::Range, sync::Arc};
use tracing::{trace, warn};

/// Score `events` of a cluster of neighbouring columns (`indices` into
/// `pileups`) together, against the reads `spanning` the cluster, and give
/// each column its records. One set of haplotypes for the cluster: a column's
/// allele that is only the aligner's other spelling of a neighbour's event
/// then loses to that event's haplotype instead of collecting its reads.
pub(crate) fn score_cluster(
    pileups: &mut [PileupMetrics],
    indices: &[usize],
    events: SmallVec<Event, MAX_EVENTS>,
    pass: ClusterPass,
    spanning: impl Iterator<Item = RecordIdx>,
    context: &ClusterContext<'_>,
    scratch: &mut ColumnScratch,
) {
    let prepared = prepare_cluster(context, spanning, events);
    let prepared = match (prepared, gpu::batch(&mut scratch.gpu, context.params.gpu)) {
        (Some(prepared), Some(batch)) => {
            match batch.push(indices, pass, prepared, pileups, &mut scratch.haplotypes) {
                None => return,
                unqueued => unqueued,
            }
        }
        (prepared, _) => prepared,
    };
    let scores = prepared.and_then(|prepared| prepared.align(&mut scratch.haplotypes));
    assign_scores(pileups, indices, pass, scores);
}

/// Scores the clusters [`score_cluster`] queued for the GPU; without a
/// GPU, nothing is queued. `pileups` must be the slice the clusters were
/// queued against.
pub(crate) fn flush_clusters(pileups: &mut [PileupMetrics], scratch: &mut ColumnScratch) {
    if let Some(batch) = scratch.gpu.as_mut() {
        batch.flush(pileups, &mut scratch.haplotypes);
    }
}

pub(super) fn assign_scores(
    pileups: &mut [PileupMetrics],
    indices: &[usize],
    pass: ClusterPass,
    scores: Option<HaplotypeScores>,
) {
    let scores = scores.map(Arc::new);
    for &index in indices {
        let Some(column) = pileups.get_mut(index) else { continue };
        let pos = column.pos;
        match pass {
            ClusterPass::Indels => {
                if let Some(indels) = column.indel_data.as_deref_mut() {
                    indels.scored = scores
                        .as_ref()
                        .filter(|scores| scores.indel_events_at(pos).next().is_some())
                        .map(|scores| Scored::new(Arc::clone(scores), pos));
                }
            }
            ClusterPass::Snvs => {
                if let Some(scores) =
                    scores.as_ref().filter(|s| s.base_events_at(pos).next().is_some())
                {
                    column.snv_scores = Some(Arc::clone(scores));
                }
            }
        }
    }
}

/// What every cluster of one sub-segment is scored against.
pub(crate) struct ClusterContext<'a> {
    pub(super) store: &'a RecordStore<RastairReadExtras>,
    pub(super) segment: &'a Segment,
    pub(super) params: &'a PileupMappingParams,
    pub(super) mask: WindowMask,
}

impl<'a> ClusterContext<'a> {
    /// `columns` are the sub-segment's drafted columns, which the mask is built from.
    pub(crate) fn new(
        store: &'a RecordStore<RastairReadExtras>,
        segment: &'a Segment,
        params: &'a PileupMappingParams,
        columns: &[PileupMetrics],
    ) -> Self {
        Self { store, segment, params, mask: WindowMask::new(columns) }
    }
}

/// Positions in a sub-segment where the drafted columns show a base other
/// than the reference on enough reads to be more than sequencing error. A
/// haplotype window scores such a position as unknown, matching any base: a
/// variant missing from every haplotype would otherwise be charged to the
/// reads of the haplotype that carries it, and the indel allele a stray read
/// suggests next to it gets the credit.
#[derive(Debug, Default)]
pub(crate) struct WindowMask(pub(super) Vec<Pos0>);

impl WindowMask {
    /// Reads of an unexplained base a column needs to be masked. Chosen on
    /// chr12 at 26x and 73x over 2, 3, 5 and 8: it is a count rather than a
    /// fraction because the FP it guards against is a handful of reads from
    /// a neighbouring variant, whatever the depth.
    const MIN_READS: u32 = 5;

    /// The masked positions of `columns`, ascending.
    pub(super) fn new(columns: &[PileupMetrics]) -> Self {
        // A `T` over a `C` on the top strand, or an `A` over a `G` on the
        // bottom, is what TAPS makes of a methylated base, not a variant. In
        // any context, not only a reference CpG: a neighbouring alt can
        // create the CpG (a `T>C` before this `G`), and masking its partner
        // hides the CpG from the haplotype, which then reads the converted
        // bases as reference (SNV F1 -0.008 at 26x, -0.017 at 73x).
        let unexplained = |reference: Base, alt: &crate::metrics::Alt| {
            CpgSide::converting(reference, alt.base).map_or(alt.metrics.depth, |side| {
                side.converting_and_plain(alt.metrics.strand_count).plain
            })
        };
        Self(
            columns
                .iter()
                .filter(|p| {
                    p.alts.iter().any(|alt| unexplained(p.reference_base, alt) >= Self::MIN_READS)
                })
                .map(|p| p.pos)
                .collect(),
        )
    }

    fn contains(&self, pos: Pos0) -> bool {
        self.0.binary_search(&pos).is_ok()
    }
}

/// `events`, each at its anchor, scored against the reads spanning the cluster.
/// The window is the reference under those reads, with the masked positions
/// unknown.
#[cfg(test)]
pub(super) fn haplotype_scores(
    context: &ClusterContext<'_>,
    spanning: impl Iterator<Item = RecordIdx>,
    events: SmallVec<Event, MAX_EVENTS>,
    scratch: &mut crate::call::pair_hmm::HaplotypeScratch,
) -> Option<HaplotypeScores> {
    prepare_cluster(context, spanning, events)?.align(scratch)
}

/// [`haplotype_scores`] up to the alignment.
fn prepare_cluster(
    context: &ClusterContext<'_>,
    spanning: impl Iterator<Item = RecordIdx>,
    events: SmallVec<Event, MAX_EVENTS>,
) -> Option<Prepared> {
    if events.is_empty() {
        return None;
    }
    let ClusterContext { store, segment, mask, params } = context;
    // A substitution scored as an event is no unexplained mismatch.
    let events_here: SmallVec<Pos0, MAX_EVENTS> = events
        .iter()
        .filter(|event| event.allele.as_indel().is_none())
        .map(|event| event.anchor)
        .collect();
    let reads = spanning_reads(store, spanning, segment);
    let window_start = reads.iter().map(|r| r.start).min()?;
    let window_end = reads.iter().map(|r| r.end).max()?;
    let reference = segment
        .sequence
        .get(segment.pos_to_idx(window_start).ok()?..=segment.pos_to_idx(window_end).ok()?)?;
    let genomic = |i: usize| window_start.checked_add_offset(Offset::new(i64::try_from(i).ok()?));
    let window: Vec<Base> = reference
        .iter()
        .enumerate()
        .map(|(i, &b)| {
            let masked =
                genomic(i).is_some_and(|pos| mask.contains(pos) && !events_here.contains(&pos));
            if masked { Base::Unknown } else { Base::from(b) }
        })
        .collect();
    let reads: Vec<_> = reads.iter().filter_map(|read| read.haplotype_read(window_start)).collect();
    pair_hmm::prepare(
        &window,
        window_start,
        events,
        reads,
        params.variant_calling.conversion_model(),
    )
}

/// A read spanning the column, placed on the genome: where its first and last
/// scored bases sit, and which of its stored bases are scored.
pub(super) struct PlacedRead<'store> {
    pub(super) rec: RecordRef<'store, RastairReadExtras>,
    strand: Strand,
    pub(super) start: Pos0,
    pub(super) end: Pos0,
    pub(super) bases: Range<usize>,
}

impl PlacedRead<'_> {
    /// Soft clips are scored when the segment's reference reaches them, since
    /// a read clipped at a tract is exactly the slipped alignment the pair-HMM
    /// can recover; otherwise the aligned bases alone, if those fit. A read
    /// of unknown strand is not placed: the emission has no strand to score
    /// it on.
    fn place<'store>(
        rec: RecordRef<'store, RastairReadExtras>,
        segment: &Segment,
    ) -> Option<PlacedRead<'store>> {
        let strand = rec.extra().strand?;
        let (lead, trail) = (u64::from(rec.leading_soft_clip), u64::from(rec.trailing_soft_clip));
        let segment_start = segment.range.region.start;
        let segment_end = segment_start.saturating_add(u64::try_from(segment.sequence.len()).ok()?);
        let inside = |start: u64, end: u64| start >= segment_start && end < segment_end;
        let (start, end) = (rec.pos.as_u64(), rec.end_pos.as_u64());
        let seq_len = u64::from(rec.seq_len);
        let clipped = start.checked_sub(lead).zip(end.checked_add(trail));
        let (start, end, bases) =
            if let Some((start, end)) = clipped.filter(|&(start, end)| inside(start, end)) {
                (start, end, 0..seq_len)
            } else if inside(start, end) {
                (start, end, lead..seq_len.saturating_sub(trail))
            } else {
                trace!(pos = start, "Haplotype read skipped: it leaves the segment");
                return None;
            };
        let bases = usize::try_from(bases.start).ok()?..usize::try_from(bases.end).ok()?;
        let (start, end) = (Pos0::try_from(start).ok()?, Pos0::try_from(end).ok()?);
        Some(PlacedRead { rec, strand, start, end, bases })
    }

    fn haplotype_read(&self, window_start: Pos0) -> Option<HaplotypeRead> {
        HaplotypeRead::new(
            self.rec.seq().get(self.bases.clone())?,
            self.rec.qual().get(self.bases.clone())?,
            i32::try_from((self.start - window_start).get()).ok().filter(|&offset| offset >= 0)?,
            self.strand,
            self.rec.qname_hash().map(FragmentId::of),
            self.rec.mapq,
        )
        .inspect_err(|error| warn!(%error, pos = %self.start, "Read not scored by the pair-HMM"))
        .ok()
    }
}

/// The reads covering the column with a bisulfite strand, one per fragment:
/// the mate reaching further past the column, which says more about what
/// follows an indel. The store holds only reads of the mapping quality the
/// run asks for (`RastairRecordFilter::filter_raw`).
pub(super) fn spanning_reads<'store>(
    store: &'store RecordStore<RastairReadExtras>,
    spanning: impl Iterator<Item = RecordIdx>,
    segment: &Segment,
) -> Vec<PlacedRead<'store>> {
    let mut reads: Vec<PlacedRead<'store>> = spanning
        .filter_map(|idx| store.record(idx))
        .filter_map(|rec| PlacedRead::place(rec, segment))
        .collect();
    reads.sort_by_key(|read| (read.rec.fragment_idx(), Reverse(read.end)));
    reads.dedup_by_key(|read| read.rec.fragment_idx());
    reads
}
