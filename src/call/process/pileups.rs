//! The segment's column loop, and the reads the SNV engine holds on to.
//!
//! **The SNV engine** (`--snv-engine`, on by default; the flag takes `=`, so a
//! positional can follow it) is a cascade:
//!
//! 1. The column loop scores only the indel clusters and keeps each
//!    sub-segment's engine (`PileupGuard::into_inner`) in [`HeldReads`],
//!    returned with the columns by `call::segment_columns`.
//! 2. `process_segment` scores every alt with the model's top-level forests.
//! 3. [`HeldReads::score_snv_clusters`] scores, as pair-HMM substitution
//!    events, the alts whose first score is at least `CASCADE_MIN_SCORE`
//!    (`from_seqair::cascade_events`), whatever their read count, so every
//!    alt the forests pass carries a verdict.
//! 4. The model's engine forests (`RastairFlatModel.snv_engine`) score exactly
//!    the alts the engine scored (`PileupMetrics::engine_scored`) again,
//!    reading its verdict as the six `ev_*` features.
//!
//! - Only engine-scored alts are rescored: the engine forests read
//!   `ev_scored` 0 as "the gate refused", so an alt rescored without a verdict
//!   is rejected. The first stage must be engine-off forests for the same
//!   reason.
//! - `ml train --snv-engine` collects twice, engine off and then as `call`
//!   scores (`call_indels` on), and fits the engine forests on the alts the
//!   engine scored. `call --snv-engine` builds the same clusters, indel calling being on
//!   in every output but CpG-only, so the `ev_*` features agree with the forests.
//! - A second scoring replaces the first, so `apply_ml_scores` sets
//!   `low_ml_score` rather than adding it, and `propagate_denovo_pass_flags`
//!   assigns its flag rather than only raising it.
//! - `snv_events_leave_the_indel_records_alone` (`tests/call_cli.rs`) and
//!   `the_engine_forests_see_only_what_the_cascade_rescores` (`train.rs`)
//!   guard the two invariants. `models/PROVENANCE.md` records how the bundled
//!   model was trained.

use crate::{
    call::{
        PreFilterInputs, RecordFilters, require_tags::TagRequirement,
        variant_calling::VariantCallingParams,
    },
    sequence::{ChunkRegion, Segment},
};
use color_eyre::eyre::Result;
use std::{ops::Deref, rc::Rc};

use crate::{
    call::{
        pair_hmm::{Event, MAX_EVENTS},
        pileup::from_seqair::{
            ClusterContext, ClusterPass, ColumnDraft, ColumnScratch, DenovoNeighbour,
            cascade_events, cluster_candidates, cluster_chains, flush_clusters,
            plain_reference_column, score_cluster, snv_clusters,
        },
    },
    metrics::{PileupMetrics, entropy::SlidingEntropy},
    sequence::{RastairReadExtras, ReferenceWindow, SeqairReaders},
};
use color_eyre::eyre::WrapErr as _;
use seqair::{
    bam::pileup::PileupEngine,
    reader::SegmentOptions,
};
use seqair_types::{Base, Pos0, SmallVec};
use std::{
    num::{NonZeroU32, NonZeroU64},
    ops::RangeInclusive,
    sync::Arc,
};
use tracing::{debug, instrument, warn};

#[derive(better_default::Default)]
pub struct PileupMappingParams {
    pub variant_calling: VariantCallingParams,
    pub require_tags: TagRequirement,
    pub call_indels: bool,
    /// Support an indel allele needs to be scored; see `IndelParams::min_indel_ao`.
    #[default(2)]
    pub min_indel_ao: u32,
    /// Per-segment compressed-byte budget; segments estimated above this are
    /// subdivided. `0` disables the budget.
    #[default(32 * 1024 * 1024)]
    pub segment_max_bytes: u64,
    /// Rescue the soft-clipped CpG-partner base adjacent to each alignment.
    /// See `MethylationCallingParams::rescue_soft_clip_cpg`.
    pub rescue_soft_clip_cpg: bool,
    /// Which columns are worth finishing.
    ///
    /// `None` finishes every column the reader emits, for a caller that wants
    /// the raw pileup. Production passes the record filters, so the seven
    /// columns in eight they drop never become a `PileupMetrics` at all.
    pub early_reject: Option<RecordFilters>,
    /// Score the pair-HMM on the GPU.
    pub gpu: bool,
    /// Keep each sub-segment's reads for [`HeldReads::score_snv_clusters`], the
    /// SNV engine's half of the cascade and the source of phasing votes.
    pub snv_engine: bool,
}

impl Deref for PileupMappingParams {
    type Target = VariantCallingParams;

    fn deref(&self) -> &Self::Target {
        &self.variant_calling
    }
}

pub fn get_pileups(
    readers: &mut SeqairReaders,
    region: &ChunkRegion,
    params: &PileupMappingParams,
) -> Result<(Rc<Segment>, impl Iterator<Item = PileupMetrics> + use<>)> {
    let (segment, columns, _held) = get_pileups_held(readers, region, params)?;
    Ok((segment, columns.into_iter()))
}

/// A region's sub-segment reads, kept past [`get_pileups_held`] under
/// `snv_engine` so its SNV clusters can be scored once the forest has ranked
/// their alts. Each sub-segment keeps its own engine: the reader can only lend
/// its store to one pileup at a time.
#[derive(Default)]
pub struct HeldReads {
    parts: Vec<HeldPart>,
}

struct HeldPart {
    engine: PileupEngine<RastairReadExtras>,
    /// First and last position of the sub-segment's kept columns.
    columns: RangeInclusive<u32>,
    /// Positions of its indel columns, which give an SNV cluster its context.
    indels: Vec<u32>,
}

impl HeldReads {
    /// Score the SNV clusters of the alts
    /// [`cascade_events`] admits, sub-segment by sub-segment. Columns are found
    /// by position, so `columns` may have lost some since they were built.
    pub fn score_snv_clusters(
        self,
        segment: &Segment,
        columns: &mut [PileupMetrics],
        params: &PileupMappingParams,
        forests: bool,
    ) {
        let mut scratch = ColumnScratch::default();
        for part in self.parts {
            let lo = columns.partition_point(|p| p.pos < *part.columns.start());
            let hi = columns.partition_point(|p| p.pos <= *part.columns.end());
            let context = ClusterContext::new(
                part.engine.store(),
                segment,
                params,
                columns.get(lo..hi).unwrap_or_default(),
            );
            let indels: Vec<usize> = part
                .indels
                .iter()
                .filter_map(|pos| columns.binary_search_by_key(pos, |p| p.pos).ok())
                .collect();
            let snvs: Vec<_> = (lo..hi)
                .filter_map(|i| Some((i, cascade_events(columns.get(i)?, forests))))
                .filter(|(_, events)| !events.is_empty())
                .collect();
            for cluster in snv_clusters(columns, &indels, &snvs, params.min_indel_ao) {
                let (indices, events) = (cluster.columns, cluster.events);
                let pass = ClusterPass::Snvs;
                score_on(columns, &indices, events, pass, &part.engine, &context, &mut scratch);
            }
            flush_clusters(columns, &mut scratch);
        }
    }
}

/// Score a cluster of `indices` (into `pileups`) against the reads `engine`
/// holds over them.
fn score_on(
    pileups: &mut [PileupMetrics],
    indices: &[usize],
    events: SmallVec<Event, MAX_EVENTS>,
    pass: ClusterPass,
    engine: &PileupEngine<RastairReadExtras>,
    context: &ClusterContext<'_>,
    scratch: &mut ColumnScratch,
) {
    let pos = |i: &usize| Pos0::new(pileups.get(*i)?.pos);
    let (Some(start), Some(end)) = (indices.first().and_then(pos), indices.last().and_then(pos))
    else {
        return;
    };
    let spanning = engine.records_overlapping(start, end);
    score_cluster(pileups, indices, events, pass, spanning, context, scratch);
}

/// [`get_pileups`], and under `snv_engine` the reads its SNV clusters need.
#[instrument(level = "info", skip_all)]
pub fn get_pileups_held(
    readers: &mut SeqairReaders,
    region: &ChunkRegion,
    params: &PileupMappingParams,
) -> Result<(Rc<Segment>, Vec<PileupMetrics>, HeldReads)> {
    // Build the rastair Segment (FASTA fetch only) and extract ref window.
    let segment = readers.segment(region, 2).wrap_err("Failed to fetch segment")?;
    debug!(len = segment.sequence.len(), "Processing region (seqair)");

    let ref_window = {
        let bases: Arc<[Base]> =
            segment.sequence.iter().map(|&b| Base::from(b)).collect::<Vec<_>>().into();
        let start = Pos0::try_from(region.start).wrap_err("region start out of Pos0 range")?;
        ReferenceWindow { bases, start }
    };

    // Install the reference window so compute() can do mismatch inference.
    readers.inner_mut().customize_mut().reference = Some(ref_window);
    readers.inner_mut().customize_mut().read_flags = params.read_flags.clone();
    readers.inner_mut().customize_mut().unpaired = params.unpaired;
    readers.inner_mut().customize_mut().min_mapq = params.quality.min_mapq;
    readers.inner_mut().customize_mut().tag_requirement = params.require_tags.clone();
    readers.inner_mut().customize_mut().guess_orientation = params.guess_read_orientation;
    readers.inner_mut().customize_mut().read_masking = params.read_masking.clone();
    readers.inner_mut().customize_mut().fetched_mapped = 0;

    let depth_limit = params.max_coverage.depth_limit();

    // Build the seqair Segment covering [region.start .. region.end + overfetch].
    // seqair's segment `end` is an *inclusive* 0-based position and must be
    // `<= contig_last_pos` (= contig length − 1). `region.last_position` is the
    // contig length, so clamp to `last_position - 1`: seqair rejects an end one
    // past the contig.
    let overfetch: u64 = 2;
    let contig_last_pos = region.last_position.saturating_sub(1);
    let last = region.end.saturating_add(overfetch).min(contig_last_pos);
    let start = Pos0::try_from(region.start).wrap_err("region start out of Pos0 range")?;
    let end = Pos0::try_from(last).wrap_err("region end out of Pos0 range")?;
    let len_u32 =
        u32::try_from(last.saturating_sub(region.start).saturating_add(1)).unwrap_or(u32::MAX);
    const ONE: NonZeroU32 = NonZeroU32::MIN;
    let max_len = NonZeroU32::new(len_u32.max(1)).unwrap_or(ONE);
    // Bound peak memory per worker: split this region into sub-segments whose
    // estimated compressed load stays within the byte budget (`0` disables it).
    let opts = match NonZeroU64::new(params.segment_max_bytes) {
        Some(budget) => SegmentOptions::new(max_len).with_max_bytes(budget),
        None => SegmentOptions::new(max_len).without_byte_budget(),
    };
    // Collect the (possibly several) sub-segments before piling up: `segments()`
    // borrows the reader immutably while `pileup()` needs it mutably. The
    // sub-segments are disjoint (overlap 0), but seqair fetches every read
    // overlapping each one, so no boundary reads are lost and no column is
    // emitted twice.
    let seqair_segments: Vec<_> = readers
        .inner_mut()
        .segments((region.contig.as_str(), (start..=end).into()), opts)
        .wrap_err("Failed to plan seqair segments")?
        .collect();

    let segment = Rc::new(segment);
    let mut columns = Columns::new(&segment, region, params);
    // Reused across every column of every sub-segment; see `ColumnDraft::accumulate`.
    let mut scratch = ColumnScratch::default();
    let mut held = HeldReads::default();

    for seqair_seg in &seqair_segments {
        // Fetch BAM records + FASTA into PileupEngine (compute() runs here).
        // The store is reused (cleared) per sub-segment, so peak memory is one
        // sub-segment's worth, not the whole region's.
        let mut guard = readers
            .inner_mut()
            .pileup(seqair_seg, depth_limit)
            .run()
            .wrap_err("Failed to start seqair pileup")?;
        if params.rescue_soft_clip_cpg {
            // Recover exactly the single CpG-partner base the aligner clipped.
            guard.set_soft_clip_overhang(1);
        }

        // Indel columns are scored once the sub-segment's columns are all
        // drafted, so the haplotype window knows the variants around them.
        let sub_segment_start = columns.kept.len();
        while let Some(col) = guard.pileups() {
            let pos = col.pos().as_u64();
            if !region.contains(pos) {
                continue;
            }
            columns.reach(pos)?;
            if let Some(neighbour) =
                plain_reference_column(&col, &segment, params, columns.previous)
            {
                columns.plain(neighbour);
                continue;
            }
            match ColumnDraft::accumulate(&col, segment.clone(), params, &mut scratch) {
                Ok(draft) => columns.draft(draft),
                Err(error) => {
                    warn!(error = format!("{error:#}"), pos, "Failed to get pileup, skipping");
                }
            }
        }

        let pending_indels = std::mem::take(&mut columns.pending_indels);
        let pileup_metrics = &mut columns.kept;
        let context = ClusterContext::new(
            guard.store(),
            &segment,
            params,
            pileup_metrics.get(sub_segment_start..).unwrap_or_default(),
        );
        for chain in cluster_chains(pileup_metrics, &pending_indels) {
            let columns = chain.iter().filter_map(|&i| pileup_metrics.get(i));
            let events = cluster_candidates(columns, params.min_indel_ao);
            let pass = ClusterPass::Indels;
            score_on(pileup_metrics, chain, events, pass, &guard, &context, &mut scratch);
        }
        // Before the next sub-segment's reads are loaded: a launch's prepared
        // reads held across that load cost 100 MB at 73x, and almost every
        // region is one sub-segment, so there is nothing to overlap with.
        flush_clusters(pileup_metrics, &mut scratch);
        drop(context);

        if params.snv_engine
            && let (Some(first), Some(last)) =
                (pileup_metrics.get(sub_segment_start), pileup_metrics.last())
        {
            held.parts.push(HeldPart {
                columns: first.pos..=last.pos,
                indels: pending_indels
                    .iter()
                    .filter_map(|&i| pileup_metrics.get(i).map(|p| p.pos))
                    .collect(),
                engine: guard.into_inner(),
            });
        }
    }

    let kept = columns.finish(region.end)?;
    Ok((segment, kept, held))
}

/// A region's columns, in ascending position, each kept or rejected as soon as
/// what decides it is known.
struct Columns<'p, 's> {
    segment: &'s Rc<Segment>,
    early_reject: Option<&'p RecordFilters>,
    /// Whether reference CpGs no read reaches get a column of their own.
    fills_gaps: bool,
    /// The first position of the region no column has reached yet.
    unseen: u64,
    /// `early_reject` drops seven covered columns in eight before anything
    /// downstream sees them, so this holds the survivors and not one entry per
    /// covered position.
    kept: Vec<PileupMetrics>,
    /// The kept columns with indel observations, as indices into `kept`, to
    /// be scored once their sub-segment's columns are all drafted.
    pending_indels: Vec<usize>,
    sliding_entropy: SlidingEntropy<'s>,
    /// The one-column delay that lets a column be rejected before it is
    /// finished. `set_denovo_adj` lets the *next* column rescue a reference C,
    /// so such a draft waits here until that column is known; every other
    /// verdict is settled by what precedes it. The draft carries the `before`
    /// it was judged against, because by the time its successor arrives
    /// `previous` has moved on to the draft itself.
    deferred: Option<(ColumnDraft, Option<DenovoNeighbour>)>,
    previous: Option<DenovoNeighbour>,
}

impl<'p, 's> Columns<'p, 's> {
    fn new(
        segment: &'s Rc<Segment>,
        region: &ChunkRegion,
        params: &'p PileupMappingParams,
    ) -> Self {
        // A quarter of the region is a deliberate over-reservation — a
        // `PileupMetrics` is 928 bytes and doubling into it both copies and
        // overshoots — against a chr12 keep rate of ~12 % without
        // `--cpgs-only` and far below that with it. A caller that asked for
        // every column doubles up into the rest, which costs it two copies of
        // a vector that was going to be that size anyway.
        //
        // `region.len()` has to be the *inclusive* count for this to land
        // where it is meant to.
        let kept = Vec::with_capacity(usize::try_from(region.len() / 4).unwrap_or(0));
        let early_reject = params.early_reject.as_ref();
        Self {
            segment,
            early_reject,
            fills_gaps: early_reject.is_some_and(RecordFilters::reports_uncovered_cpgs),
            unseen: region.start,
            kept,
            pending_indels: Vec::new(),
            sliding_entropy: SlidingEntropy::new(segment),
            deferred: None,
            previous: None,
        }
    }

    fn keeps(&self, inputs: PreFilterInputs) -> bool {
        self.early_reject.is_none_or(|f| f.keeps(inputs))
    }

    /// A column [`plain_reference_column`] settled without a draft: it is never
    /// kept, but it may still rescue the deferred one.
    fn plain(&mut self, neighbour: DenovoNeighbour) {
        if let Some((waiting, before)) = self.deferred.take()
            && self.keeps(waiting.pre_filter_inputs(before, Some(neighbour)))
        {
            self.keep(waiting);
        }
        self.previous = Some(neighbour);
    }

    fn draft(&mut self, draft: ColumnDraft) {
        let neighbour = draft.denovo_neighbour();
        // The deferred column comes first: output stays in ascending position
        // order.
        if let Some((waiting, before)) = self.deferred.take()
            && self.keeps(waiting.pre_filter_inputs(before, Some(neighbour)))
        {
            self.keep(waiting);
        }
        if self.keeps(draft.pre_filter_inputs(self.previous, None)) {
            if self.keep(draft).is_some_and(|kept| kept.indel_data.is_some()) {
                self.pending_indels.push(self.kept.len().saturating_sub(1));
            }
        } else if draft.awaits_successor() {
            self.deferred = Some((draft, self.previous));
        }
        self.previous = Some(neighbour);
    }

    /// Move on to the column at `pos`, first giving each reference CpG before
    /// it that no read reaches a column of its own.
    fn reach(&mut self, pos: u64) -> Result<()> {
        if self.fills_gaps {
            for gap in self.unseen..pos {
                if let Some(draft) = ColumnDraft::uncovered_cpg(gap, self.segment)? {
                    self.draft(draft);
                }
            }
        }
        self.unseen = pos.saturating_add(1);
        Ok(())
    }

    /// The kept columns, once every position before `end` has been reached.
    fn finish(mut self, end: u64) -> Result<Vec<PileupMetrics>> {
        self.reach(end)?;
        // Nothing follows the last column, so a still-deferred draft was only
        // ever waiting for a rescue that cannot come.
        if let Some((waiting, before)) = self.deferred.take()
            && self.keeps(waiting.pre_filter_inputs(before, None))
        {
            self.keep(waiting);
        }
        Ok(self.kept)
    }

    /// Finish a surviving column and give it its region entropy.
    ///
    /// The entropy is only read by the ML feature extractor, so it is computed
    /// here and not for every covered position; the sliding window is still fed
    /// ascending indices and its counts are integers, so the values are the ones
    /// the whole-region pass produced.
    ///
    /// A column kept only once its successor is known (a reference C the next
    /// column makes a de-novo CpG) is finished after its reads are gone, so it
    /// cannot be scored. Outside CpG-only output an indel column is always kept
    /// at once, so only CpG-only output loses scores this way.
    fn keep(&mut self, draft: ColumnDraft) -> Option<&mut PileupMetrics> {
        let idx = draft.idx();
        match draft.finish() {
            Ok(mut metrics) => {
                metrics.pos_metrics.extended.region_entropy = self.sliding_entropy.entropy_at(idx);
                self.kept.push(metrics);
                self.kept.last_mut()
            }
            Err(error) => {
                warn!(error = format!("{error:#}"), "Failed to finish pileup, skipping");
                None
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::ReaderParams;
    use color_eyre::eyre::Result;

    #[test]
    fn test_reading_bounds() -> Result<()> {
        // The region is 1-based and closed; the segment's core is exactly it,
        // 0-based, while its overlap reaches past both ends.
        let params = ReaderParams {
            regions: Some("chr19:6105700-6105800".parse().unwrap()),
            ..ReaderParams::test_data()
        };
        let mut readers = params.pileup_readers()?;
        let segments: Vec<_> = readers.segments(10_000, 100)?.collect();
        readers.segment(&segments[0], 0)?;

        let pileup_mapping_params = PileupMappingParams::default();
        let (segment, pileups) = get_pileups(&mut readers, &segments[0], &pileup_mapping_params)?;
        let core: Vec<_> = pileups.filter(|p| segment.is_core(u64::from(p.pos))).collect();

        assert!(!core.is_empty());
        assert_eq!(core.first().unwrap().pos, 6_105_699);
        assert_eq!(core.last().unwrap().pos, 6_105_799);

        Ok(())
    }

    /// The early rejection is exactly [`RecordFilters::pre_filter`] moved in
    /// front of the work it saves. Both halves of that claim are checked: the
    /// surviving positions are the ones the downstream `retain` would have
    /// kept, in the same order, and their region entropy is what the
    /// whole-region sliding pass produced — the window only ever sees the kept
    /// columns now.
    ///
    /// `set_denovo_adj` is the reason this is not obvious: a column with no
    /// evidence of its own is kept when its genomic neighbour carries an alt
    /// that would create the other half of a CpG, so the reference run has to
    /// go through `map_surrounding` before it filters.
    #[test]
    fn early_rejection_keeps_exactly_what_the_pre_filter_would() -> Result<()> {
        use crate::{call::RecordFilters, utils::map_surrounding};

        let params = ReaderParams {
            regions: Some("chr19:6105700-6106500".parse()?),
            ..ReaderParams::test_data()
        };
        let mut readers = params.pileup_readers()?;
        let segments: Vec<_> = readers.segments(10_000, 100)?.collect();

        let mut compared = 0usize;
        for filters in [
            RecordFilters { vcf_all: false, cpgs_only: false },
            RecordFilters { vcf_all: true, cpgs_only: false },
            RecordFilters { vcf_all: false, cpgs_only: true },
        ] {
            let unfiltered = PileupMappingParams::default();
            for chunk in &segments {
                let (_s, everything) = get_pileups(&mut readers, chunk, &unfiltered)?;
                let mut everything: Vec<_> = everything.collect();
                map_surrounding(&mut everything, super::super::set_denovo_adj, "test mapper");
                everything.retain(|p| filters.pre_filter(p));
                let expected: Vec<_> = everything
                    .iter()
                    .map(|p| (p.pos, p.pos_metrics.extended.region_entropy))
                    .collect();

                let early = PileupMappingParams {
                    early_reject: Some(filters.clone()),
                    ..Default::default()
                };
                let (_s, kept) = get_pileups(&mut readers, chunk, &early)?;
                let actual: Vec<_> =
                    kept.map(|p| (p.pos, p.pos_metrics.extended.region_entropy)).collect();

                assert_eq!(expected, actual, "early rejection diverged for {filters:?}");
                compared += expected.len();
            }
        }
        assert!(compared > 0, "the fixture produced nothing to compare");
        Ok(())
    }

    /// The cascade hands the engine only the alts the forest scored high
    /// enough: the column loop scores no SNV cluster, and once every gated alt
    /// has a stage-1 score, those at `CASCADE_MIN_SCORE` or more get an engine
    /// verdict (`ev_scored` 1) and those below none.
    #[test]
    fn the_cascade_scores_only_alts_the_forest_kept() -> Result<()> {
        use crate::{
            call::{
                pair_hmm::Allele,
                pileup::from_seqair::{CASCADE_MIN_SCORE, snv_events},
            },
            metrics::{MetricsForAlt, ml::features::shared::CommonFeatures},
        };
        use seqair_types::Probability;

        let params =
            PileupMappingParams { snv_engine: true, call_indels: true, ..Default::default() };
        let ev_scored = |column: &PileupMetrics, base: Base| {
            column.alts.iter().find(|alt| alt.base == base).map(|alt| {
                CommonFeatures::extract(&MetricsForAlt { metrics: column, alt: &alt.metrics })
                    .engine[0]
            })
        };
        let mut readers = ReaderParams::test_data().pileup_readers()?;
        let regions: Vec<_> = readers.segments(100_000, 0)?.collect();
        let (mut inside, mut outside) = (0usize, 0usize);
        for region in &regions {
            let (segment, mut columns, held) = get_pileups_held(&mut readers, region, &params)?;
            assert!(columns.iter().all(|c| c.snv_scores.is_none()), "the loop scored an SNV");
            let mut gated = Vec::new();
            for (i, column) in columns.iter_mut().enumerate() {
                for (k, event) in snv_events(column).into_iter().enumerate() {
                    let Allele::Base(base) = event.allele else { continue };
                    let kept = (i + k) % 2 == 0;
                    let score = if kept { 0.5 } else { CASCADE_MIN_SCORE / 2.0 };
                    if let Some(alt) = column.alts.iter_mut().find(|alt| alt.base == base) {
                        alt.filters.ml = Some(Probability::try_from(score)?);
                        gated.push((i, base, kept));
                    }
                }
            }
            held.score_snv_clusters(&segment, &mut columns, &params, true);
            for (i, base, kept) in gated {
                let scored = columns.get(i).and_then(|column| ev_scored(column, base));
                if kept {
                    inside += usize::from(scored == Some(1.0));
                } else {
                    assert_eq!(scored, Some(0.0), "an alt under the floor was scored");
                    outside += 1;
                }
            }
            if inside > 0 && outside > 0 {
                return Ok(());
            }
        }
        panic!("the fixture has no gated alt on both sides of the floor: {inside} / {outside}");
    }

    /// `PASS ⇒ engine-scored`: an alt the first stage scores high enough to
    /// pass carries an engine verdict whatever its reads, the one-read and
    /// strand-skewed alts the gate refuses included, so its genotype is the
    /// engine's. Only an unknown reference base is never scored.
    #[test]
    fn every_alt_the_forest_can_pass_is_engine_scored() -> Result<()> {
        use crate::call::{pair_hmm::Allele, pileup::from_seqair::snv_events};
        use seqair_types::Probability;

        let params =
            PileupMappingParams { snv_engine: true, call_indels: true, ..Default::default() };
        let reader_params = ReaderParams {
            regions: Some("chr19:6105000-6107000".parse()?),
            ..ReaderParams::test_data()
        };
        let mut readers = reader_params.pileup_readers()?;
        let regions: Vec<_> = readers.segments(100_000, 0)?.collect();
        let (mut scored, mut refused_by_gate) = (0usize, 0usize);
        for region in &regions {
            let (segment, mut columns, held) = get_pileups_held(&mut readers, region, &params)?;
            let mut passing = Vec::new();
            for (i, column) in columns.iter_mut().enumerate() {
                let gated = snv_events(column);
                let known = column.reference_base.known_index().is_some();
                let mut alts: Vec<_> = column.alts.iter_mut().collect();
                alts.sort_by_key(|alt| (std::cmp::Reverse(alt.metrics.depth), alt.base as u8));
                for alt in alts.into_iter().filter(|alt| alt.base.known_index().is_some()).take(2) {
                    alt.filters.ml = Some(Probability::try_from(0.9)?);
                    let by_gate = gated.iter().any(|e| e.allele == Allele::Base(alt.base));
                    passing.push((i, alt.base, known, by_gate));
                }
            }
            held.score_snv_clusters(&segment, &mut columns, &params, true);
            for (i, base, known, by_gate) in passing {
                let engine = columns.get(i).is_some_and(|column| column.engine_scored(base));
                assert_eq!(engine, known, "alt {base:?} at column {i}");
                scored += usize::from(engine);
                refused_by_gate += usize::from(engine && !by_gate);
            }
            if scored > 0 && refused_by_gate > 0 {
                return Ok(());
            }
        }
        panic!("no alt the gate refuses was passed: {scored} scored, {refused_by_gate} refused");
    }

    /// Byte-aware sub-segmentation must be transparent: a tiny budget (which
    /// forces the region to be split into many sub-segments) yields exactly the
    /// same pileups as an effectively-unlimited budget.
    #[test]
    fn byte_budget_subdivision_is_transparent() -> Result<()> {
        let params = ReaderParams {
            regions: Some("chr19:6105700-6105900".parse().unwrap()),
            ..ReaderParams::test_data()
        };
        let mut readers = params.pileup_readers()?;
        let segments: Vec<_> = readers.segments(10_000, 100)?.collect();

        // Everything a boundary read could move: depth, the ALT list *in order*
        // (it becomes the VCF ALT column), per-allele depth and strand split,
        // and the indel observation count. A read fetched twice, or dropped at
        // a sub-segment edge, shows up in one of these.
        #[allow(clippy::type_complexity)]
        fn fingerprint(
            pileups: impl Iterator<Item = PileupMetrics>,
        ) -> Vec<(u32, Base, usize, Vec<(Base, u32, u32, u32)>, u32)> {
            pileups
                .map(|p| {
                    let alleles = std::iter::once(&p.ref_metrics)
                        .chain(p.alts.iter().map(|alt| &alt.metrics))
                        .map(|m| (m.base, m.depth, m.strand_count.ot, m.strand_count.ob))
                        .collect();
                    let indels =
                        p.indel_data.as_ref().map_or(0, |data| data.counts.total_indel_reads());
                    (p.pos, p.reference_base, p.pos_metrics.depth as usize, alleles, indels)
                })
                .collect()
        }

        let huge = PileupMappingParams { segment_max_bytes: u64::MAX, ..Default::default() };
        let (_s1, p1) = get_pileups(&mut readers, &segments[0], &huge)?;
        let baseline = fingerprint(p1);

        // A 1-byte budget forces the region to split into many sub-segments.
        let tiny = PileupMappingParams { segment_max_bytes: 1, ..Default::default() };
        let (_s2, p2) = get_pileups(&mut readers, &segments[0], &tiny)?;
        let subdivided = fingerprint(p2);

        assert!(!baseline.is_empty());
        assert!(
            baseline.iter().any(|(_, _, _, alleles, _)| alleles.len() > 1),
            "the region must contain alts, or the ALT comparison proves nothing"
        );
        assert_eq!(baseline, subdivided, "subdivision changed the pileups");
        Ok(())
    }
}
