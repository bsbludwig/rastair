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
//!    (`from_seqair::cascade_events`, the two highest), whatever their read
//!    count, so an alt the forests pass carries a verdict unless it is a
//!    third, the reference is `N`, or no read can be scored on it.
//! 4. The model's engine forests (`RastairFlatModel.snv_engine`) score exactly
//!    the alts the engine scored (`PileupMetrics::engine_scored`) again,
//!    reading its verdict as the five `ev_*` features.
//!
//! - Only engine-scored alts are rescored: the engine forests are fitted on
//!   those alone, and an alt without a verdict has every `ev_*` at 0, a row
//!   they never saw. The first stage scores before the engine, so its forests
//!   see the same zeros in training and never split on them; both stages
//!   share one feature layout.
//! - `ml train --snv-engine` collects twice, engine off and then as `call`
//!   scores (`Engine::PairHmm`), and fits the engine forests on the alts the
//!   engine scored. `call --snv-engine` builds the same clusters, indel
//!   calling on, so the `ev_*` features agree with the forests.
//! - A second scoring replaces the first score: `extract_ml_rows` clears it
//!   on every alt the second takes, so an alt the engine forests leave
//!   unscored (a `NaN` feature) is rejected, never called on its first score.
//! - `snv_events_leave_the_indel_records_alone` (`tests/call_cli.rs`) and
//!   `the_engine_forests_see_only_what_the_cascade_rescores` (`train.rs`)
//!   guard the two invariants. `models/PROVENANCE.md` records how the bundled
//!   model was trained.

use crate::{
    call::{
        ColumnEvidence,
        pair_hmm::{Event, MAX_EVENTS},
        pileup::from_seqair::{
            CandidateFloor, ClusterContext, ClusterPass, ColumnDraft, ColumnScratch,
            DenovoNeighbour, SnvCluster, SnvColumn, cascade_events, cluster_candidates,
            cluster_chains, flush_clusters, plain_reference_column, score_cluster, snv_clusters,
        },
        require_tags::TagRequirement,
        variant_calling::{IndelParams, VariantCallingParams},
    },
    metrics::{PileupMetrics, entropy::SlidingEntropy},
    sequence::{
        ByteBudget, ChunkRegion, IndexDepth, RastairReadExtras, ReferenceWindow, Segment,
        SeqairReaders,
    },
};
use color_eyre::eyre::{Result, WrapErr as _};
use enumset::EnumSet;
use seqair::{bam::pileup::PileupEngine, reader::SegmentOptions};
use seqair_types::{Base, Pos0, SmallVec};
use std::{
    num::NonZeroU32,
    ops::{Deref, RangeInclusive},
    rc::Rc,
    sync::Arc,
};
use tracing::{debug, instrument, warn};

#[derive(better_default::Default)]
pub struct PileupMappingParams {
    pub variant_calling: VariantCallingParams,
    pub require_tags: TagRequirement,
    pub engine: Engine,
    /// Support an indel allele needs to be scored; see `IndelParams::min_indel_ao`.
    #[default(IndelParams::default().min_indel_ao)]
    pub indel_floor: CandidateFloor,
    /// Per-segment compressed-byte budget; segments estimated above this are
    /// subdivided.
    #[default(ByteBudget::new(32 * 1024 * 1024))]
    pub segment_max_bytes: ByteBudget,
    /// Rescue the soft-clipped CpG-partner base adjacent to each alignment
    /// See `MethylationCallingParams::rescue_soft_clip_cpg`.
    pub rescue_soft_clip_cpg: bool,
    pub columns: ColumnKeep,
    pub uncovered_cpgs: UncoveredCpgs,
    /// Score the pair-HMM on the GPU.
    pub gpu: bool,
}

/// Whether the column loop runs the pair-HMM.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Engine {
    /// The columns alone: training's first stage, which has no forest to rank
    /// alts yet.
    #[default]
    ColumnsOnly,
    /// The indel clusters are scored in the column loop, and each
    /// sub-segment's reads are kept for [`HeldReads::score_snv_clusters`]: the
    /// engine genotypes every called SNV and casts the phasing votes from them.
    PairHmm,
}

/// Which columns the column loop finishes into a [`PileupMetrics`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ColumnKeep {
    /// Every column the reader emits: the raw pileup.
    #[default]
    Every,
    /// Those with [`ColumnEvidence`]: the seven covered columns in eight
    /// without never become a `PileupMetrics` at all.
    Candidates,
}

/// Whether reference CpGs no read reaches get a column of their own.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum UncoveredCpgs {
    #[default]
    Skip,
    /// So that two samples' `--all` outputs list the same positions.
    Fill,
}

impl Deref for PileupMappingParams {
    type Target = VariantCallingParams;

    fn deref(&self) -> &Self::Target {
        &self.variant_calling
    }
}

#[cfg(test)]
pub fn get_pileups(
    readers: &mut SeqairReaders,
    region: &ChunkRegion,
    params: &PileupMappingParams,
) -> Result<(Rc<Segment>, impl Iterator<Item = PileupMetrics> + use<>)> {
    let (segment, columns, _held) = get_pileups_held(readers, region, params)?;
    Ok((segment, columns.into_iter()))
}

/// A region's sub-segment reads, kept past [`get_pileups_held`] under
/// [`Engine::PairHmm`] so its SNV clusters can be scored once the forest has ranked
/// their alts. Each sub-segment keeps its own engine: the reader can only lend
/// its store to one pileup at a time. So all of a region's stores are held at
/// once, and `--segment-max-bytes` does not bound its peak memory.
#[derive(Default)]
pub struct HeldReads {
    parts: Vec<HeldPart>,
    /// The region's depth, for the phasing depth guard. A read crossing a
    /// sub-segment boundary counts once per sub-segment.
    pub depth: IndexDepth,
}

struct HeldPart {
    engine: PileupEngine<RastairReadExtras>,
    /// First and last position of the sub-segment's kept columns.
    columns: RangeInclusive<Pos0>,
    /// Its indel chains, with the events the indel pass scored: an SNV
    /// cluster near one takes them as context.
    indel_chains: Vec<IndelChain>,
}

/// An indel cluster's columns, by position, and the events scored over them.
struct IndelChain {
    positions: Vec<Pos0>,
    events: SmallVec<Event, MAX_EVENTS>,
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
            let indel_chains = part
                .indel_chains
                .into_iter()
                .map(|IndelChain { positions, events }| SnvCluster {
                    columns: positions
                        .iter()
                        .filter_map(|pos| columns.binary_search_by_key(pos, |p| p.pos).ok())
                        .collect(),
                    scored: Vec::new(),
                    events,
                })
                .collect();
            let snvs: Vec<_> = (lo..hi)
                .filter_map(|column| {
                    Some(SnvColumn { column, events: cascade_events(columns.get(column)?) })
                })
                .filter(|snv| !snv.events.is_empty())
                .collect();
            for SnvCluster { columns: span, scored, events } in
                snv_clusters(columns, indel_chains, &snvs)
            {
                let pass = ClusterPass::Snvs;
                let engine = &part.engine;
                score_on(columns, &span, &scored, events, pass, engine, &context, &mut scratch);
            }
            flush_clusters(columns, &mut scratch);
        }
    }
}

/// Score a cluster over `span` (indices into `pileups`, ascending) against
/// the reads `engine` holds over them, for the columns `scored`.
#[expect(clippy::too_many_arguments, reason = "the cluster, its reads and the scratch")]
fn score_on(
    pileups: &mut [PileupMetrics],
    span: &[usize],
    scored: &[usize],
    events: SmallVec<Event, MAX_EVENTS>,
    pass: ClusterPass,
    engine: &PileupEngine<RastairReadExtras>,
    context: &ClusterContext<'_>,
    scratch: &mut ColumnScratch,
) {
    let pos = |i: &usize| Some(pileups.get(*i)?.pos);
    let (Some(start), Some(end)) = (span.first().and_then(pos), span.last().and_then(pos)) else {
        return;
    };
    let spanning = engine.records_overlapping(start, end);
    score_cluster(pileups, scored, events, pass, spanning, context, scratch);
}

/// [`get_pileups`], and under [`Engine::PairHmm`] the reads its SNV clusters
/// need.
#[instrument(level = "info", skip_all)]
pub fn get_pileups_held(
    readers: &mut SeqairReaders,
    region: &ChunkRegion,
    params: &PileupMappingParams,
) -> Result<(Rc<Segment>, Vec<PileupMetrics>, HeldReads)> {
    // Build the rastair Segment (FASTA fetch only) and extract ref window.
    let segment = readers.segment(region, 2).wrap_err("Failed to fetch segment")?;
    debug!(len = segment.sequence.len(), "Processing region");

    let ref_window = {
        let bases: Arc<[Base]> =
            segment.sequence.iter().map(|&b| Base::from(b)).collect::<Vec<_>>().into();
        let start = Pos0::try_from(region.start).wrap_err("region start out of Pos0 range")?;
        ReferenceWindow { bases, start }
    };

    // Install the reference window so compute() can do mismatch inference.
    let customize = readers.inner_mut().customize_mut();
    customize.reference = Some(ref_window);
    customize.read_flags = params.read_flags.clone();
    customize.pairing = params.pairing();
    customize.min_mapq = params.quality.min_mapq;
    customize.tag_requirement = params.require_tags.clone();
    customize.strand_source = params.strand_source();
    customize.read_masking = params.read_masking.clone();

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
    // estimated compressed load stays within the byte budget.
    let opts = match params.segment_max_bytes.get() {
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
    // A CRAM has no bins, and no run depth to compare with, so any unit does.
    let leaf_bin_len = readers.leaf_bin_len().unwrap_or(1).max(1);
    let mut last_bin = None;

    for seqair_seg in &seqair_segments {
        // Fetch BAM records + FASTA into PileupEngine (compute() runs here).
        // Under `Engine::PairHmm` each sub-segment's store is kept in `held` for
        // the SNV engine, so the region's reads are all in memory at once;
        // without it a store is reused (cleared) per sub-segment.
        let mut guard = readers
            .inner_mut()
            .pileup(seqair_seg, depth_limit)
            .run()
            .wrap_err("Failed to start seqair pileup")?;
        held.depth.mapped += guard.fetch_counts().mapped as u64;
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
            let bin = pos / leaf_bin_len;
            if last_bin.is_none_or(|last| bin > last) {
                last_bin = Some(bin);
                let bin_start = bin.saturating_mul(leaf_bin_len);
                let bin_end = bin_start.saturating_add(leaf_bin_len).min(region.end);
                held.depth.covered_bases += bin_end.saturating_sub(bin_start.max(region.start));
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
        let mut indel_chains = Vec::new();
        for chain in cluster_chains(pileup_metrics, &pending_indels) {
            let columns = chain.iter().filter_map(|&i| pileup_metrics.get(i));
            let events = cluster_candidates(columns, params.indel_floor);
            if params.engine == Engine::PairHmm {
                let positions = chain.iter().filter_map(|&i| pileup_metrics.get(i).map(|p| p.pos));
                indel_chains
                    .push(IndelChain { positions: positions.collect(), events: events.clone() });
            }
            let pass = ClusterPass::Indels;
            score_on(pileup_metrics, chain, chain, events, pass, &guard, &context, &mut scratch);
        }
        // Before the next sub-segment's reads are loaded: a launch's prepared
        // reads held across that load cost 100 MB at 73x, and almost every
        // region is one sub-segment, so there is nothing to overlap with.
        flush_clusters(pileup_metrics, &mut scratch);
        drop(context);

        if params.engine == Engine::PairHmm
            && let (Some(first), Some(last)) =
                (pileup_metrics.get(sub_segment_start), pileup_metrics.last())
        {
            held.parts.push(HeldPart {
                columns: first.pos..=last.pos,
                indel_chains,
                engine: guard.into_inner(),
            });
        }
    }

    let kept = columns.finish(region.end)?;
    Ok((segment, kept, held))
}

/// A region's columns, in ascending position, each kept or rejected as soon as
/// what decides it is known.
struct Columns<'s> {
    segment: &'s Rc<Segment>,
    keep: ColumnKeep,
    uncovered_cpgs: UncoveredCpgs,
    /// The first position of the region no column has reached yet.
    unseen: u64,
    /// [`ColumnKeep::Candidates`] drops seven covered columns in eight before
    /// anything downstream sees them, so this holds the survivors and not one
    /// entry per covered position.
    kept: Vec<PileupMetrics>,
    /// The kept columns with indel observations, as indices into `kept`, to
    /// be scored once their sub-segment's columns are all drafted.
    pending_indels: Vec<usize>,
    sliding_entropy: SlidingEntropy<'s>,
    /// The one-column delay that lets a column be rejected before it is
    /// finished. [`crate::vcf::CpgSide::beside`] lets the *next*
    /// column rescue a reference C, so such a draft waits here until that
    /// column is known; every other verdict is settled by what precedes it.
    deferred: Option<Deferred>,
    previous: Option<DenovoNeighbour>,
}

/// A draft waiting for its successor, with the `before` it was judged
/// against: by the time the successor arrives `previous` has moved on to the
/// draft itself.
struct Deferred {
    draft: ColumnDraft,
    before: Option<DenovoNeighbour>,
}

impl<'s> Columns<'s> {
    fn new(segment: &'s Rc<Segment>, region: &ChunkRegion, params: &PileupMappingParams) -> Self {
        // A quarter of the region is a deliberate over-reservation — a
        // `PileupMetrics` is over 500 bytes and doubling into it both copies and
        // overshoots — against a chr12 keep rate of ~12 %. A caller that asked for
        // every column doubles up into the rest, which costs it two copies of
        // a vector that was going to be that size anyway.
        //
        // `region.len()` has to be the *inclusive* count for this to land
        // where it is meant to.
        let kept = Vec::with_capacity(usize::try_from(region.len() / 4).unwrap_or(0));
        Self {
            segment,
            keep: params.columns,
            uncovered_cpgs: params.uncovered_cpgs,
            unseen: region.start,
            kept,
            pending_indels: Vec::new(),
            sliding_entropy: SlidingEntropy::new(segment),
            deferred: None,
            previous: None,
        }
    }

    fn keeps(&self, evidence: EnumSet<ColumnEvidence>) -> bool {
        match self.keep {
            ColumnKeep::Every => true,
            ColumnKeep::Candidates => !evidence.is_empty(),
        }
    }

    /// Keep the deferred draft, if any, when its successor `after` rescues it.
    fn settle_deferred(&mut self, after: Option<DenovoNeighbour>) {
        if let Some(Deferred { draft, before }) = self.deferred.take()
            && self.keeps(draft.evidence(before, after))
        {
            self.keep(draft);
        }
    }

    /// A column [`plain_reference_column`] settled without a draft: it is never
    /// kept, but it may still rescue the deferred one.
    fn plain(&mut self, neighbour: DenovoNeighbour) {
        self.settle_deferred(Some(neighbour));
        self.previous = Some(neighbour);
    }

    fn draft(&mut self, draft: ColumnDraft) {
        let neighbour = draft.denovo_neighbour();
        // The deferred column comes first: output stays in ascending position
        // order.
        self.settle_deferred(Some(neighbour));
        if self.keeps(draft.evidence(self.previous, None)) {
            if self.keep(draft).is_some_and(|kept| kept.indel_data.is_some()) {
                self.pending_indels.push(self.kept.len().saturating_sub(1));
            }
        } else if draft.awaits_successor() {
            self.deferred = Some(Deferred { draft, before: self.previous });
        }
        self.previous = Some(neighbour);
    }

    /// Move on to the column at `pos`, first giving each reference CpG before
    /// it that no read reaches a column of its own.
    fn reach(&mut self, pos: u64) -> Result<()> {
        if self.uncovered_cpgs == UncoveredCpgs::Fill {
            for gap in self.unseen..pos {
                let gap = Pos0::try_from(gap).wrap_err("position out of Pos0 range")?;
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
        self.settle_deferred(None);
        Ok(self.kept)
    }

    /// Finish a surviving column with its region entropy.
    ///
    /// The entropy is computed here and not for every covered position; the
    /// sliding window is still fed ascending indices and its counts are
    /// integers, so the values are the ones a whole-region pass produces.
    ///
    /// A column kept only once its successor is known (a reference C the next
    /// column makes a de-novo CpG) is finished after its reads are gone, so it
    /// cannot be scored. An indel column is always kept at once, so none loses
    /// its scores this way.
    fn keep(&mut self, draft: ColumnDraft) -> Option<&mut PileupMetrics> {
        let region_entropy = self.sliding_entropy.entropy_at(draft.idx());
        match draft.finish(region_entropy) {
            Ok(metrics) => {
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
        let core: Vec<_> = pileups.filter(|p| segment.is_core(p.pos)).collect();

        assert!(!core.is_empty());
        assert_eq!(core.first().unwrap().pos.as_u32(), 6_105_699);
        assert_eq!(core.last().unwrap().pos.as_u32(), 6_105_799);

        Ok(())
    }

    /// The early rejection keeps exactly what [`crate::call::pre_filter`]
    /// keeps of every column. Both halves of that claim are checked: the
    /// surviving positions are the same, in the same order, and their region
    /// entropy is what the whole-region sliding pass produced — the window
    /// only ever sees the kept columns.
    ///
    /// [`crate::vcf::CpgSide::beside`] is the reason this is not obvious: a
    /// column with no evidence of its own is kept when its genomic neighbour
    /// carries an alt that would create the other half of a CpG, so the
    /// column loop decides a reference C only once its successor is known.
    #[test]
    fn early_rejection_keeps_exactly_what_the_pre_filter_would() -> Result<()> {
        use crate::call::RecordFilters;

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
                crate::call::pre_filter(&mut everything);
                let expected: Vec<_> =
                    everything.iter().map(|p| (p.pos, p.pos_metrics.region_entropy)).collect();

                let early = PileupMappingParams {
                    columns: ColumnKeep::Candidates,
                    uncovered_cpgs: filters.uncovered_cpgs(),
                    ..Default::default()
                };
                let (_s, kept) = get_pileups(&mut readers, chunk, &early)?;
                let actual: Vec<_> = kept.map(|p| (p.pos, p.pos_metrics.region_entropy)).collect();

                assert_eq!(expected, actual, "early rejection diverged for {filters:?}");
                compared += expected.len();
            }
        }
        assert!(compared > 0, "the fixture produced nothing to compare");
        Ok(())
    }

    /// The cascade hands the engine only the alts the forest scored high
    /// enough: the column loop scores no SNV cluster, and once every alt has a
    /// stage-1 score, those at `CASCADE_MIN_SCORE` or more get an engine
    /// verdict and those below none.
    #[test]
    fn the_cascade_scores_only_alts_the_forest_kept() -> Result<()> {
        use crate::{
            call::pileup::from_seqair::CASCADE_MIN_SCORE,
            metrics::{MetricsForAlt, ml::features::shared::common_sections},
        };
        use seqair_types::Probability;

        let params = PileupMappingParams { engine: Engine::PairHmm, ..Default::default() };
        let engine = |column: &PileupMetrics, base: Base| {
            let alt = column.alts.iter().find(|alt| alt.base == base)?;
            let (_, section) =
                common_sections(&MetricsForAlt { metrics: column, alt: &alt.metrics });
            Some((column.engine_scored(base), section.engine))
        };
        let reader_params = ReaderParams {
            regions: Some("chr19:6105000-6107000".parse()?),
            ..ReaderParams::test_data()
        };
        let mut readers = reader_params.pileup_readers()?;
        let regions: Vec<_> = readers.segments(100_000, 0)?.collect();
        let (mut inside, mut outside) = (0usize, 0usize);
        for region in &regions {
            let (segment, mut columns, held) = get_pileups_held(&mut readers, region, &params)?;
            assert!(columns.iter().all(|c| c.snv_scores.is_none()), "the loop scored an SNV");
            let mut first = Vec::new();
            for (i, column) in columns.iter_mut().enumerate() {
                let known = column.alts.iter_mut().filter(|alt| alt.base.known_index().is_some());
                for (k, alt) in known.enumerate() {
                    let kept = (i + k) % 2 == 0;
                    let score = if kept { 0.5 } else { CASCADE_MIN_SCORE / 2.0 };
                    alt.verdict = crate::metrics::AltVerdict::Scored(Probability::try_from(score)?);
                    first.push((i, alt.base, kept));
                }
            }
            held.score_snv_clusters(&segment, &mut columns, &params);
            for (i, base, kept) in first {
                let scored = columns.get(i).and_then(|column| engine(column, base));
                if kept {
                    inside += usize::from(scored.is_some_and(|(scored, _)| scored));
                } else {
                    assert_eq!(
                        scored,
                        Some((false, [0.0; 5])),
                        "an alt under the floor was scored"
                    );
                    outside += 1;
                }
            }
            if inside > 0 && outside > 0 {
                return Ok(());
            }
        }
        panic!("the fixture has no alt on both sides of the floor: {inside} / {outside}");
    }

    /// `PASS ⇒ engine-scored`: an alt the first stage scores high enough to
    /// pass carries an engine verdict whatever its reads, one-read alts
    /// included, so its genotype is the engine's. Only an unknown reference
    /// base, and a column's third such alt, are never scored.
    #[test]
    fn every_alt_the_forest_can_pass_is_engine_scored() -> Result<()> {
        use seqair_types::Probability;

        let params = PileupMappingParams { engine: Engine::PairHmm, ..Default::default() };
        let reader_params = ReaderParams {
            regions: Some("chr19:6105000-6107000".parse()?),
            ..ReaderParams::test_data()
        };
        let mut readers = reader_params.pileup_readers()?;
        let regions: Vec<_> = readers.segments(100_000, 0)?.collect();
        let (mut scored, mut one_read) = (0usize, 0usize);
        for region in &regions {
            let (segment, mut columns, held) = get_pileups_held(&mut readers, region, &params)?;
            let mut passing = Vec::new();
            for (i, column) in columns.iter_mut().enumerate() {
                let known = column.reference_base.known_index().is_some();
                let mut alts: Vec<_> = column.alts.iter_mut().collect();
                alts.sort_by_key(|alt| (std::cmp::Reverse(alt.metrics.depth), alt.base as u8));
                for alt in alts.into_iter().filter(|alt| alt.base.known_index().is_some()).take(2) {
                    alt.verdict = crate::metrics::AltVerdict::Scored(Probability::try_from(0.9)?);
                    passing.push((i, alt.base, known, alt.metrics.depth == 1));
                }
            }
            held.score_snv_clusters(&segment, &mut columns, &params);
            for (i, base, known, single) in passing {
                let engine = columns.get(i).is_some_and(|column| column.engine_scored(base));
                assert_eq!(engine, known, "alt {base:?} at column {i}");
                scored += usize::from(engine);
                one_read += usize::from(engine && single);
            }
            if scored > 0 && one_read > 0 {
                return Ok(());
            }
        }
        panic!("no one-read alt was passed: {scored} scored, {one_read} with one read");
    }

    /// A stretch without reads, like an assembly gap, leaves the segment's
    /// depth alone: `test.bam`'s chr19 reads end at 6,143,150, inside the leaf
    /// bin ending at 6,144,000, so both regions hold the same reads over the
    /// same bases.
    #[test]
    fn an_uncovered_stretch_does_not_dilute_the_segment_depth() -> Result<()> {
        let held = |region: &str| -> Result<(u64, u64)> {
            let params =
                ReaderParams { regions: Some(region.parse()?), ..ReaderParams::test_data() };
            let mut readers = params.pileup_readers()?;
            let segments: Vec<_> = readers.segments(100_000, 0)?.collect();
            let [segment] = segments.as_slice() else {
                panic!("{region} is one segment, got {}", segments.len());
            };
            let (_, _, held) =
                get_pileups_held(&mut readers, segment, &PileupMappingParams::default())?;
            Ok((held.depth.mapped, held.depth.covered_bases))
        };
        let short = held("chr19:6140001-6145000")?;
        let (mapped, covered) = short;
        assert!(mapped > 0 && covered == 4_000, "{short:?}");
        assert_eq!(held("chr19:6140001-6150000")?, short);
        Ok(())
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
                    (
                        p.pos.as_u32(),
                        p.reference_base,
                        p.pos_metrics.depth as usize,
                        alleles,
                        indels,
                    )
                })
                .collect()
        }

        let huge = PileupMappingParams {
            segment_max_bytes: ByteBudget::new(u64::MAX),
            ..Default::default()
        };
        let (_s1, p1) = get_pileups(&mut readers, &segments[0], &huge)?;
        let baseline = fingerprint(p1);

        // A 1-byte budget forces the region to split into many sub-segments.
        let tiny =
            PileupMappingParams { segment_max_bytes: ByteBudget::new(1), ..Default::default() };
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
