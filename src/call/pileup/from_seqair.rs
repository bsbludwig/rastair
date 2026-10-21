//! Columns from seqair's pileup engine, accumulated into
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
//! - Tract runs are measured one base past the anchor, where a left-aligned
//!   indel starts; only `ref_features::indel_tract_runs_at` knows this.
//! - Overlap dedup keeps one read per fragment and can keep the mate that does
//!   not span the indel, so the dropped mate stands in:
//!   `build_indel_observation(view).or_else(mate)`. The fallback is keyed on
//!   the observation, because `build_indel_observation` rejects indels near a
//!   read end or on reads with too many non-TAPS mismatches; seqair's
//!   `pair_indel` ("own indel wins") cannot express that, so it is not used.
//!   On chr12 this is worth +2.3 insertion / +1.5 deletion F1.
//! - `soft_clip_count` and the noisy-reference count describe a fragment:
//!   whatever either mate shows counts once, and a fragment is noisy
//!   reference only when it contributed no observation and a mate that is a
//!   reference alignment slips (terminal repeat or soft clip).
//!
//! Testing: `tests/data/test.bam` cannot produce an indel call (no two of its
//! indel reads agree on an allele), so indel tests build their own reads; a
//! fixture segment must not start at 0, where position bugs and the right
//! answer coincide (use `segment_at`); and a unit test must set
//! `params.call_indels` to reach the indel branch at all.

use super::indels::{IndelAllele, IndelAlleleCounts, IndelCounts, Scored};
use crate::{
    call::{
        PreFilterInputs,
        pair_hmm::{
            self, Allele, Cluster, Event, HaplotypeRead, HaplotypeScores, HaplotypeScratch,
            MAX_EVENTS, Prepared,
        },
        phasing::{self, ColumnPhase, FragmentId, PhaseObservation},
        process::PileupMappingParams,
    },
    metrics::{
        Alt, AltFilters, Filters, FormsDenovo, PairedCounts, PerBaseAccumulators, PileupMetrics,
        ReadKey, RecordTags, alt_forms_denovo, order_alts,
    },
    sequence::{RastairReadExtras, Segment},
    utils::SequenceContext,
    vcf::InCpG,
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use seqair::bam::{
    RecordIdx, RecordRef,
    pileup::{AlignmentView, Indel, PileupColumn, PileupOp},
    record_store::RecordStore,
};
use seqair_types::{Base, QPos, SmallVec, Strand, SumOfSquares};
use std::{cmp::Reverse, ops::Range, rc::Rc, sync::Arc};
use tracing::{instrument, trace, warn};

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
    /// Gathered here because the buffer the observations come from is reset at
    /// the next column.
    phase: Option<Box<ColumnPhase>>,
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

        let dedup_overlaps = !params.keep_overlapping_reads;
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
                None => observed(&view, params, reference_base, &context),
            };
            let Some(Observed { base, baseq, qpos }) = resolved else {
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
                    scratch,
                )
            {
                continue;
            }
            total_depth += 1;

            let qual = seqair_types::BaseQuality::from_byte(baseq);
            if params.phase && phasing::accepts_read(view.mapq, qual) {
                match view.qname_hash().map(FragmentId::of) {
                    Some(fragment) => scratch.phase_observations.push(PhaseObservation {
                        fragment,
                        base,
                        strand,
                        qual,
                    }),
                    None => scratch.saw_nameless_read = true,
                }
            }

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
                let mate = dedup_overlaps
                    .then(|| {
                        counted_mate(
                            column,
                            &view,
                            params,
                            reference_base,
                            &context,
                            &mut scratch.mates,
                        )
                    })
                    .flatten();
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

        // Without `--phase` nothing was collected, so nothing is a candidate.
        let phase = phasing::candidate_observations(&scratch.phase_observations, reference_base)
            .map(|observations| Box::new(ColumnPhase::Observed(observations)));

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
            phase,
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
        let cpg =
            InCpG::new(reference_base, base_at(idx.checked_sub(1)), base_at(idx.checked_add(1)));
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
            phase: None,
        }))
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
            phase,
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
            phase,
            snv_scores: None,
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
        scratch: &mut ColumnScratch,
    ) -> Result<PileupMetrics> {
        let draft = ColumnDraft::accumulate(column, Rc::clone(&segment), params, scratch)?;
        let mut metrics = draft.finish()?;
        let pos = column.pos();
        let spanning = column.records_overlapping(pos, pos);
        let context = ClusterContext::new(column.store(), &segment, params, &[]);
        let column = std::slice::from_mut(&mut metrics);
        let events = cluster_candidates(column.iter(), params.min_indel_ao);
        score_cluster(column, &[0], events, ClusterPass::Indels, spanning, &context, scratch);
        flush_clusters(column, scratch);
        Ok(metrics)
    }
}

/// Columns closer than this are one cluster: their alleles are often the
/// aligner's different spellings of one tract event, or two events one
/// chromosome carries, and either needs them scored together. 6, 12 and 20 bp
/// are alike on chr12 at both depths, 20 a shade better at 73x.
const CLUSTER_GAP: u32 = 20;

/// Which of a column's results a cluster's scores become.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ClusterPass {
    /// The indel records: a cluster of indel events only.
    Indels,
    /// The SNV engine's features: a cluster with substitution events, and any
    /// indel events in it only as context. It never writes an indel record, so
    /// an SNV event cannot change an indel call.
    Snvs,
}

/// `columns` (indices into `pileups`, ascending) cut into chains of
/// neighbours within [`CLUSTER_GAP`].
pub(crate) fn cluster_chains<'a>(
    pileups: &[PileupMetrics],
    columns: &'a [usize],
) -> Vec<&'a [usize]> {
    let pos = |i: usize| pileups.get(i).map(|p| p.pos);
    columns
        .chunk_by(|&a, &b| pos(a).zip(pos(b)).is_some_and(|(a, b)| b.abs_diff(a) <= CLUSTER_GAP))
        .collect()
}

/// An SNV cluster to score: its columns (indices into the pileups,
/// ascending) and its events in record order, the indel events first.
#[derive(Debug, Default, PartialEq)]
pub(crate) struct SnvCluster {
    pub(crate) columns: Vec<usize>,
    pub(crate) events: SmallVec<Event, MAX_EVENTS>,
}

/// The SNV clusters of a sub-segment. An SNV column within [`CLUSTER_GAP`] of
/// a chain of indel columns joins the chain's indel events while they admit
/// its events, so it is scored against the haplotypes of its indel
/// neighbours; the other SNV columns chain among themselves, cut wherever a
/// column's events no longer fit. `snvs` are SNV columns with their events,
/// ascending; a column's events stay together.
pub(crate) fn snv_clusters(
    pileups: &[PileupMetrics],
    indel_columns: &[usize],
    snvs: &[(usize, SmallVec<Event, 2>)],
    min_indel_ao: u32,
) -> Vec<SnvCluster> {
    let pos = |i: usize| pileups.get(i).map(|p| p.pos);
    let near =
        |a: usize, b: usize| pos(a).zip(pos(b)).is_some_and(|(a, b)| b.abs_diff(a) <= CLUSTER_GAP);
    let mut joined = vec![false; snvs.len()];
    let mut clusters = Vec::new();
    for chain in cluster_chains(pileups, indel_columns) {
        let (Some(&first), Some(&last)) = (chain.first(), chain.last()) else { continue };
        let indels = cluster_candidates(chain.iter().filter_map(|&i| pileups.get(i)), min_indel_ao);
        let mut cluster = SnvCluster { columns: chain.to_vec(), events: indels };
        let mut any = false;
        for ((column, snv), joined) in snvs.iter().zip(&mut joined) {
            let fits = |events: &SmallVec<Event, MAX_EVENTS>| {
                let mut trial = events.clone();
                snv.iter().all(|event| {
                    let admitted = Cluster::admits(&trial, event);
                    trial.push(event.clone());
                    admitted
                })
            };
            if !*joined
                && (near(first, *column) || near(*column, last) || (first..=last).contains(column))
                && fits(&cluster.events)
            {
                cluster.events.extend(snv.iter().cloned());
                if !cluster.columns.contains(column) {
                    cluster.columns.push(*column);
                }
                *joined = true;
                any = true;
            }
        }
        if any {
            cluster.columns.sort_unstable();
            clusters.push(cluster);
        }
    }

    let rest: Vec<&(usize, SmallVec<Event, 2>)> =
        snvs.iter().zip(&joined).filter(|(_, joined)| !**joined).map(|(snv, _)| snv).collect();
    for chain in rest.chunk_by(|a, b| near(a.0, b.0)) {
        let mut current = SnvCluster::default();
        for (column, snv) in chain {
            if current.events.len() + snv.len() > MAX_EVENTS {
                clusters.push(std::mem::take(&mut current));
            }
            current.columns.push(*column);
            current.events.extend(snv.iter().cloned());
        }
        if !current.events.is_empty() {
            clusters.push(current);
        }
    }
    clusters
}

/// Allele fraction the heterozygote test uses, under one half so that an
/// allele the informative strand under-observes still passes.
const HET_GATE_FRACTION: f64 = 0.4;
/// A column whose informative strand leaves a heterozygote this improbable
/// admits no SNV event.
const HET_GATE_P: f64 = 0.01;
/// Fewer informative reads than this exclude nothing.
const HET_GATE_MIN_READS: u32 = 3;

/// The column's substitution events, most reads first, at most two (as for
/// any anchor). An alt is admitted on Ben Wulf's gate: at least two reads (and
/// 5 % of the column, as an indel allele needs), and on the strand TAPS leaves
/// informative for it, a count a heterozygote could produce. Without the gate
/// every methylated `CpG` and every one-read error would be a cluster.
pub(crate) fn snv_events(column: &PileupMetrics) -> SmallVec<Event, 2> {
    if column.reference_base.known_index().is_none() {
        return SmallVec::new();
    }
    let floor = 2u32.max((column.pos_metrics.depth * 5).div_ceil(100));
    let mut admitted: SmallVec<&Alt, 4> = column
        .alts
        .iter()
        .filter(|alt| alt.base.known_index().is_some() && alt.base != column.reference_base)
        .filter(|alt| alt.metrics.depth >= floor && admits_het(column, alt))
        .collect();
    admitted.sort_by_key(|alt| (Reverse(alt.metrics.depth), alt.base as u8));
    admitted
        .into_iter()
        .take(2)
        .map(|alt| Event { anchor: column.pos, allele: Allele::Base(alt.base) })
        .collect()
}

/// The lowest first-stage forest score at which the cascade hands an alt the
/// gate admits to the engine. The bundled forests put most true SNVs at
/// .995-.999, so there is no upper edge; .005 keeps 70-110 % of the gain of
/// scoring every gated alt for a quarter of its clusters on 26x slices, .01
/// loses some.
pub(crate) const CASCADE_MIN_SCORE: f64 = 0.005;

/// [`snv_events`] of the alts the forest scored at [`CASCADE_MIN_SCORE`] or more.
pub(crate) fn cascade_events(column: &PileupMetrics) -> SmallVec<Event, 2> {
    let mut events = snv_events(column);
    events.retain(|event| {
        matches!(event.allele, Allele::Base(base) if column
            .alt_filters(base)
            .and_then(|filters| filters.ml)
            .is_some_and(|ml| *ml >= CASCADE_MIN_SCORE))
    });
    events
}

/// Whether the strand TAPS conversion leaves informative for `alt` does not
/// rule a heterozygote out: a `T` over a `C` may be a converted methylated
/// `C` on OT and an `A` over a `G` on OB, and a true `C` or `G` allele is
/// itself converted there.
fn admits_het(column: &PileupMetrics, alt: &Alt) -> bool {
    let informative = match (column.reference_base, alt.base) {
        (Base::C, Base::T) | (_, Base::C) => Some(Strand::OB),
        (Base::G, Base::A) | (_, Base::G) => Some(Strand::OT),
        _ => None,
    };
    let on = |metrics: &crate::metrics::AlleleMetrics| match informative {
        Some(Strand::OT) => metrics.strand_count.ot,
        Some(_) => metrics.strand_count.ob,
        None => metrics.depth,
    };
    let reads: u32 = column.ref_alts_metrics().map(on).sum();
    reads < HET_GATE_MIN_READS
        || binomial_cdf(on(&alt.metrics), reads, HET_GATE_FRACTION) >= HET_GATE_P
}

/// `P(X <= k)` for `X ~ Binomial(n, p)`, summed upwards from the zero term.
fn binomial_cdf(k: u32, n: u32, p: f64) -> f64 {
    let mut term = (1.0 - p).powi(i32::try_from(n).unwrap_or(i32::MAX));
    let mut total = term;
    for i in 1..=k.min(n) {
        term *= f64::from(n - i + 1) / f64::from(i) * p / (1.0 - p);
        total += term;
    }
    total.min(1.0)
}

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

fn assign_scores(
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
    store: &'a RecordStore<RastairReadExtras>,
    segment: &'a Segment,
    params: &'a PileupMappingParams,
    mask: WindowMask,
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
pub(crate) struct WindowMask(Vec<u32>);

impl WindowMask {
    /// Reads of an unexplained base a column needs to be masked. Chosen on
    /// chr12 at 26x and 73x over 2, 3, 5 and 8: it is a count rather than a
    /// fraction because the FP it guards against is a handful of reads from
    /// a neighbouring variant, whatever the depth.
    const MIN_READS: u32 = 5;

    /// The masked positions of `columns`, ascending.
    fn new(columns: &[PileupMetrics]) -> Self {
        // A `T` over a `C` on the top strand, or an `A` over a `G` on the
        // bottom, is what TAPS makes of a methylated base, not a variant.
        let unexplained = |reference: Base, alt: &crate::metrics::Alt| match (reference, alt.base) {
            (Base::C, Base::T) => alt.metrics.strand_count.ob,
            (Base::G, Base::A) => alt.metrics.strand_count.ot,
            _ => alt.metrics.depth,
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

    fn contains(&self, pos: u32) -> bool {
        self.0.binary_search(&pos).is_ok()
    }
}

/// The cluster's events: the alleles the most reads carry, of those with
/// enough to be called, as many as [`Cluster::admits`], in record order. Ties
/// go by the allele itself, so neither the choice nor the order depends on
/// which read the engine yielded first.
pub(crate) fn cluster_candidates<'a>(
    columns: impl Iterator<Item = &'a PileupMetrics>,
    min_ao: u32,
) -> SmallVec<Event, MAX_EVENTS> {
    let mut alleles: SmallVec<(u32, &IndelAlleleCounts), 8> = columns
        .filter_map(|column| {
            // Two reads of a stutter allele are noise at 73x and the only
            // evidence at 14x, so the floor also grows with the column: 5 %
            // of its depth drops a quarter of the candidate events and 17 %
            // of the pair-HMM cells at 73x and changes no column of the chr12
            // scores at either depth; 2 % drops nothing.
            let floor = min_ao.max((column.pos_metrics.depth * 5).div_ceil(100));
            Some((column.pos, floor, &column.indel_data.as_deref()?.counts))
        })
        .flat_map(|(pos, floor, counts)| {
            counts.alleles.iter().filter(move |a| a.total >= floor).map(move |a| (pos, a))
        })
        .collect();
    alleles.sort_by_key(|(pos, a)| (Reverse(a.total), *pos, allele_order(&a.allele)));
    let mut events: SmallVec<Event, MAX_EVENTS> = SmallVec::new();
    for (anchor, counted) in alleles {
        let event = Event { anchor, allele: counted.allele.clone().into() };
        if Cluster::admits(&events, &event) {
            events.push(event);
        }
    }
    events
}

fn allele_order(allele: &IndelAllele) -> (u8, SmallVec<u8, 4>) {
    match allele {
        IndelAllele::Insertion(bases) => (0, bases.iter().map(|b| *b as u8).collect()),
        IndelAllele::Deletion(bases) => (1, bases.iter().map(|b| *b as u8).collect()),
    }
}

/// `events`, each at its anchor, scored against the reads spanning the cluster.
/// The window is the reference under those reads, with the masked positions
/// unknown.
#[cfg(test)]
fn haplotype_scores(
    context: &ClusterContext<'_>,
    spanning: impl Iterator<Item = RecordIdx>,
    events: SmallVec<Event, MAX_EVENTS>,
    scratch: &mut HaplotypeScratch,
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
    let ClusterContext { store, segment, params, mask } = context;
    // A substitution scored as an event is no unexplained mismatch.
    let events_here: SmallVec<u32, MAX_EVENTS> = events
        .iter()
        .filter(|event| event.allele.as_indel().is_none())
        .map(|event| event.anchor)
        .collect();
    let reads = spanning_reads(store, spanning, segment, params);
    let window_start = reads.iter().map(|r| r.start).min()?;
    let window_end = reads.iter().map(|r| r.end).max()?;
    let segment_start = segment.range.region.start;
    let offset = |genomic: u64| usize::try_from(genomic.checked_sub(segment_start)?).ok();
    let reference = segment.sequence.get(offset(window_start)?..=offset(window_end)?)?;
    let genomic =
        |i: usize| u32::try_from(window_start.saturating_add(i as u64)).unwrap_or(u32::MAX);
    let window: Vec<Base> = reference
        .iter()
        .enumerate()
        .map(|(i, &b)| {
            let pos = genomic(i);
            if mask.contains(pos) && !events_here.contains(&pos) {
                Base::Unknown
            } else {
                Base::from(b)
            }
        })
        .collect();
    let reads: Vec<_> = reads.iter().filter_map(|read| read.haplotype_read(window_start)).collect();
    pair_hmm::prepare(&window, window_start, events, reads)
}

/// A read spanning the column, placed on the genome: where its first and last
/// scored bases sit, and which of its stored bases are scored.
struct PlacedRead<'store> {
    rec: RecordRef<'store, RastairReadExtras>,
    start: u64,
    end: u64,
    bases: Range<usize>,
}

impl PlacedRead<'_> {
    /// Soft clips are scored when the segment's reference reaches them, since
    /// a read clipped at a tract is exactly the slipped alignment the pair-HMM
    /// can recover; otherwise the aligned bases alone, if those fit.
    fn place<'store>(
        rec: RecordRef<'store, RastairReadExtras>,
        segment: &Segment,
    ) -> Option<PlacedRead<'store>> {
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
        Some(PlacedRead { rec, start, end, bases })
    }

    fn haplotype_read(&self, window_start: u64) -> Option<HaplotypeRead> {
        HaplotypeRead::new(
            self.rec.seq().get(self.bases.clone())?,
            self.rec.qual().get(self.bases.clone())?,
            i32::try_from(self.start.checked_sub(window_start)?).ok()?,
            self.rec.extra().strand,
            self.rec.qname_hash().map(FragmentId::of),
            self.rec.mapq,
        )
        .inspect_err(|error| warn!(%error, pos = self.start, "Read not scored by the pair-HMM"))
        .ok()
    }
}

/// The reads covering the column with a bisulfite strand and the mapping
/// quality the column asks for, one per fragment: the mate reaching further
/// past the column, which says more about what follows an indel.
fn spanning_reads<'store>(
    store: &'store RecordStore<RastairReadExtras>,
    spanning: impl Iterator<Item = RecordIdx>,
    segment: &Segment,
    params: &PileupMappingParams,
) -> Vec<PlacedRead<'store>> {
    let mut reads: Vec<PlacedRead<'store>> = spanning
        .filter_map(|idx| store.record(idx))
        .filter(|rec| rec.mapq >= params.quality.min_mapq && rec.extra().strand != Strand::Unknown)
        .filter_map(|rec| PlacedRead::place(rec, segment))
        .collect();
    reads.sort_by_key(|read| (read.rec.fragment_idx(), Reverse(read.end)));
    reads.dedup_by_key(|read| read.rec.fragment_idx());
    reads
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
    scratch: &mut ColumnScratch,
) -> bool {
    if !view.in_mate_overlap() {
        return false;
    }
    let this_idx = view.alignment().record_idx();
    let Some(mate_idx) = view.alignment().mate_idx() else { return false };

    if this_idx > mate_idx {
        return scratch.mate_drops.binary_search(&this_idx).is_ok();
    }

    // A mate that is absent from this column, or that fails a filter here,
    // never formed a pair — this read stands on its own.
    let Some(slot) = column.position_of(mate_idx) else { return false };
    let Some(mate) = scratch.mates.resolve(column, slot, params, reference_base, context) else {
        return false;
    };
    let Some(Observed { base: mate_base, .. }) = mate.observed else {
        return false;
    };

    // The name-based collector resolved the pair when it reached the *later*
    // read: it dropped that read if the bases agreed or it was read 2, and
    // otherwise dropped the earlier one.
    if base != mate_base && !mate.second_in_template {
        return true;
    }
    if let Err(slot) = scratch.mate_drops.binary_search(&mate_idx) {
        scratch.mate_drops.insert(slot, mate_idx);
    }
    false
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
    reference_base: Base,
    context: &SequenceContext,
) -> Option<Observed> {
    let baseq = view.qual()?.get()?;
    let base = view.base()?;
    let qpos = view.qpos()?;
    if baseq < params.quality.min_baseq {
        return None;
    }
    if !passes_read_masking(view, reference_base, context) {
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
    /// The current column's phase observations, before the candidate test
    /// decides whether they are kept.
    phase_observations: Vec<PhaseObservation>,
    /// Whether a phase observation was dropped because its read has no qname.
    pub(crate) saw_nameless_read: bool,
    haplotypes: HaplotypeScratch,
    /// `None` until the first cluster, and after it when the pair-HMM is
    /// scored on the CPU.
    gpu: Option<gpu::Batch>,
}

impl ColumnScratch {
    fn begin(&mut self, depth: usize) {
        self.mate_drops.clear();
        self.phase_observations.clear();
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
/// quality filtering, read-end masking and, for a rescued clip, a CpG-partner
/// test each time.
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
        reference_base: Base,
        context: &SequenceContext,
    ) -> Option<MateObservation> {
        if let Some(known) = self.slots.get(slot).filter(|s| s.generation == self.generation) {
            return Some(*known);
        }
        let mate = column.alignment_at(slot)?;
        let entry = MateObservation {
            generation: self.generation,
            second_in_template: mate.flags.is_second_in_template(),
            observed: observed(&mate, params, reference_base, context),
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
    reference_base: Base,
    context: &SequenceContext,
    mates: &mut MateObservations,
) -> Option<AlignmentView<'a, 'eng, RastairReadExtras>> {
    // Same reject `mate_of` applies before its search: outside the pair's
    // overlap the mate cannot be in this column.
    if !view.in_mate_overlap() {
        return None;
    }
    let slot = column.position_of(view.mate_idx()?)?;
    let entry = mates.resolve(column, slot, params, reference_base, context)?;
    entry.observed.is_some().then(|| column.alignment_at(slot)).flatten()
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

/// The pair-HMM on the GPU: a region's clusters queued into launches, each
/// run by the worker's own GPU thread while the worker prepares the next.
mod gpu {
    use super::{ClusterPass, HaplotypeScratch, PileupMetrics, Prepared, assign_scores};
    use compair::gpu::{GpuAligner, GpuContext, GpuError, GpuPairs, GpuScore, KernelOptions, Wait};
    use std::{
        cell::OnceCell,
        mem,
        ops::Range,
        sync::{
            Arc, OnceLock,
            atomic::{AtomicU64, Ordering},
            mpsc,
        },
        thread,
    };
    use tracing::{debug, info, warn};

    /// Pairs a launch holds at most. A launch costs ~1.5 ms over its pairs,
    /// so smaller launches lose that; larger ones hold more prepared reads
    /// in memory, twice over while one runs and the next fills. 4,096 and
    /// 8,192 time the same on chr12 at 73x; 4,096 holds 120 MB less.
    const MAX_LAUNCH_PAIRS: usize = 4_096;

    fn context() -> Option<&'static Arc<GpuContext>> {
        static CONTEXT: OnceLock<Option<Arc<GpuContext>>> = OnceLock::new();
        CONTEXT
            .get_or_init(|| match GpuContext::new() {
                Ok(context) => {
                    info!(adapter = %context.adapter_info().name, "The pair-HMM is scored on the GPU");
                    Some(context)
                }
                Err(error) => {
                    warn!(%error, "No GPU for the pair-HMM; it is scored on the CPU");
                    None
                }
            })
            .as_ref()
    }

    type Launched = (GpuPairs, Result<Vec<GpuScore>, GpuError>);

    /// A worker's aligner on a thread of its own: `submit` a plan, `collect`
    /// it back with its scores. The aligner's buffers grow to the largest
    /// launch and stay, so it outlives the region's [`Batch`].
    struct Lane {
        submit: mpsc::Sender<GpuPairs>,
        collect: mpsc::Receiver<Launched>,
    }

    thread_local! {
        static LANE: OnceCell<Option<Lane>> = const { OnceCell::new() };
    }

    impl Lane {
        fn spawn(context: &Arc<GpuContext>) -> std::io::Result<Self> {
            let (submit, launches) = mpsc::channel::<GpuPairs>();
            let (results, collect) = mpsc::channel();
            // Sleep between polls rather than spin a core while a launch runs.
            let options = KernelOptions { wait: Wait::Block, ..KernelOptions::default() };
            let context = Arc::clone(context);
            let aligner = move || GpuAligner::with_options(Arc::clone(&context), options);
            let name = thread::current().name().unwrap_or("worker").to_owned();
            thread::Builder::new().name(format!("{name}-gpu")).spawn(move || {
                let mut current = aligner();
                for pairs in launches {
                    let scores = current.align(&pairs);
                    // A failed launch can leave the aligner's buffers mapped;
                    // a fresh one costs a regrowth and keeps the lane alive.
                    if scores.is_err() {
                        current = aligner();
                    }
                    if results.send((pairs, scores)).is_err() {
                        break;
                    }
                }
            })?;
            Ok(Self { submit, collect })
        }

        fn with<R>(f: impl FnOnce(Option<&Lane>) -> R) -> R {
            LANE.with(|lane| {
                f(lane
                    .get_or_init(|| {
                        Lane::spawn(context()?)
                            .inspect_err(|error| {
                                warn!(%error, "No GPU thread for the pair-HMM; it is scored on the CPU");
                            })
                            .ok()
                    })
                    .as_ref())
            })
        }
    }

    /// The region's batch, created on first use; `None` when the pair-HMM
    /// is scored on the CPU.
    pub(super) fn batch(slot: &mut Option<Batch>, gpu: bool) -> Option<&mut Batch> {
        if gpu && slot.is_none() && Lane::with(|lane| lane.is_some()) {
            *slot = Some(Batch::default());
        }
        slot.as_mut()
    }

    /// A queued cluster: its columns, its prepared state, and where its
    /// scores land in the launch.
    struct Pending {
        columns: Vec<usize>,
        pass: ClusterPass,
        prepared: Prepared,
        range: Range<usize>,
    }

    #[derive(Default)]
    pub(in super::super) struct Batch {
        filling: GpuPairs,
        pending: Vec<Pending>,
        /// The clusters of the launch on the lane, if one is.
        in_flight: Option<Vec<Pending>>,
    }

    static LAUNCHES: AtomicU64 = AtomicU64::new(0);
    static PAIRS: AtomicU64 = AtomicU64::new(0);
    static RESCORED: AtomicU64 = AtomicU64::new(0);

    impl Batch {
        /// Queues a cluster, launching the batch once it is full; hands the
        /// cluster back when it cannot be queued, to score on the CPU.
        pub(super) fn push(
            &mut self,
            columns: &[usize],
            pass: ClusterPass,
            prepared: Prepared,
            pileups: &mut [PileupMetrics],
            scratch: &mut HaplotypeScratch,
        ) -> Option<Prepared> {
            match prepared.push(&mut self.filling) {
                Ok(range) => {
                    self.pending.push(Pending { columns: columns.to_vec(), pass, prepared, range });
                    if self.filling.len() >= MAX_LAUNCH_PAIRS {
                        self.launch(pileups, scratch);
                    }
                    None
                }
                Err(error) => {
                    warn!(%error, "A cluster did not fit the GPU launch; it is scored on the CPU");
                    Some(prepared)
                }
            }
        }

        /// Scores everything queued, waiting for the GPU.
        pub(super) fn flush(
            &mut self,
            pileups: &mut [PileupMetrics],
            scratch: &mut HaplotypeScratch,
        ) {
            if !self.pending.is_empty() {
                self.launch(pileups, scratch);
            }
            if self.in_flight.is_some() {
                self.filling = self.collect(pileups, scratch);
            }
        }

        /// Submits the filling plan, after collecting the launch before it.
        fn launch(&mut self, pileups: &mut [PileupMetrics], scratch: &mut HaplotypeScratch) {
            let spare = if self.in_flight.is_some() {
                self.collect(pileups, scratch)
            } else {
                GpuPairs::new()
            };
            let plan = mem::replace(&mut self.filling, spare);
            let pending = mem::take(&mut self.pending);
            match Lane::with(|lane| lane.map(|lane| lane.submit.send(plan))) {
                Some(Ok(())) => self.in_flight = Some(pending),
                _ => {
                    warn!("The pair-HMM's GPU thread is gone; its clusters are scored on the CPU");
                    settle(pending, None, pileups, scratch);
                }
            }
        }

        /// Waits for the launch on the lane, assigns its scores, and returns
        /// its plan, cleared, for the next launch.
        fn collect(
            &mut self,
            pileups: &mut [PileupMetrics],
            scratch: &mut HaplotypeScratch,
        ) -> GpuPairs {
            let pending = self.in_flight.take().unwrap_or_default();
            let Some((mut plan, scores)) = Lane::with(|lane| lane?.collect.recv().ok()) else {
                warn!("The pair-HMM's GPU thread is gone; its clusters are scored on the CPU");
                settle(pending, None, pileups, scratch);
                return GpuPairs::new();
            };
            let scores = scores
                .inspect_err(|error| {
                    warn!(%error, "The pair-HMM's GPU launch failed; its clusters are scored on the CPU");
                })
                .ok();
            if let Some(scores) = &scores {
                let rescored = scores.iter().filter(|score| score.trusted().is_none()).count();
                let total =
                    |counter: &AtomicU64, n: u64| counter.fetch_add(n, Ordering::Relaxed) + n;
                debug!(
                    launches = total(&LAUNCHES, 1),
                    pairs = total(&PAIRS, scores.len() as u64),
                    rescored = total(&RESCORED, rescored as u64),
                    launch_pairs = scores.len(),
                    "Pair-HMM GPU launch (running totals)"
                );
            }
            settle(pending, scores.as_deref(), pileups, scratch);
            plan.clear();
            plan
        }
    }

    /// A batch dropped with a launch in flight (a region that failed) takes
    /// the launch's result off the lane, or the next region would collect it.
    impl Drop for Batch {
        fn drop(&mut self) {
            if self.in_flight.is_some()
                && !Lane::with(|lane| lane.is_some_and(|lane| lane.collect.recv().is_ok()))
            {
                debug!("A pair-HMM launch was lost with its GPU thread");
            }
        }
    }

    /// Gives each cluster its scores from the launch, or scores it on the
    /// CPU where the launch has none.
    fn settle(
        pending: Vec<Pending>,
        scores: Option<&[GpuScore]>,
        pileups: &mut [PileupMetrics],
        scratch: &mut HaplotypeScratch,
    ) {
        for Pending { columns, pass, prepared, range } in pending {
            let result = match scores.and_then(|scores| scores.get(range)) {
                Some(scores) => prepared.finish_gpu(scores),
                None => prepared.align(scratch),
            };
            assign_scores(pileups, &columns, pass, result);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::pileup::indels::IndelData;
    use crate::call::process::PileupMappingParams;
    use crate::call::variant_calling::{ReadMaskParams, ReadMaskSetting};
    use crate::sequence::RastairRecordFilter;
    use crate::sequence::{ChunkRegion, Region, Segment};
    use crate::utils::default;
    use seqair::bam::DepthCap;
    use seqair::bam::cigar::{CigarOp, CigarOpType};
    use seqair::bam::pileup::PileupEngine;
    use seqair::bam::record_store::{
        CustomizeRecordStore, FilterRawFields, RecordStore, SlimRecord,
    };
    use seqair::bam::{BamHeader, MateFields};
    use seqair::reader::{ResolveTid, Tid};
    use seqair_types::{BamFlags, Pos0, Strand};
    use std::num::NonZeroU32;

    /// The only reference the test reads are on; mates are only linked when
    /// both records name the same one.
    fn tid0() -> Tid {
        0u32.resolve_tid(&BamHeader::synthetic(1, 1_000)).unwrap()
    }

    /// A mate on the test reference at `pos`; a negative `pos` means no mate position.
    fn mate_at(pos: i32) -> MateFields {
        MateFields {
            tid: Some(tid0()),
            pos: u32::try_from(pos).ok().and_then(Pos0::new),
            template_len: None,
        }
    }

    /// The parts of `RastairRecordFilter` these tests need: the MAPQ filter,
    /// strand from the flags, and the read-end mask window.
    #[derive(Default, Clone)]
    struct TestExtras {
        masking: ReadMaskParams,
        min_mapq: u8,
    }

    impl TestExtras {
        fn new(params: &PileupMappingParams) -> Self {
            Self { masking: params.read_masking.clone(), min_mapq: params.quality.min_mapq }
        }
    }

    impl CustomizeRecordStore for TestExtras {
        type Extra = RastairReadExtras;

        fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
            fields.mapq >= self.min_mapq
        }

        fn compute(
            &mut self,
            rec: &SlimRecord,
            _store: &RecordStore<RastairReadExtras>,
        ) -> RastairReadExtras {
            let strand = Strand::from(rec.flags);
            RastairReadExtras {
                strand,
                mask: self
                    .masking
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
    fn store_of(
        reads: &[TestRead],
        params: &PileupMappingParams,
    ) -> RecordStore<RastairReadExtras> {
        let mut extras = TestExtras::new(params);
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
                    Some(tid0()),
                    mate_at(mate_pos),
                    &mut extras,
                )
                .unwrap();
        }
        // Linking (and sorting) is `prepare_for_pileup`'s job at the point the
        // engine is built, so reads may be listed in any order.
        store
    }

    /// The overlapping-pair rule by read name: group the column's surviving
    /// alignments by qname and, on the second one, drop whichever
    /// `resolve_pair` chose. An independent oracle — it reaches the answer by
    /// matching names, which is exactly what the mate links do not do.
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
                if baseq < params.quality.min_baseq {
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
        let store = store_of(reads, &params);
        let last = seq.len().saturating_sub(1) as u32;
        let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, last));
        if params.rescue_soft_clip_cpg {
            engine.set_soft_clip_overhang(1);
        }

        let mut scratch = ColumnScratch::default();
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
                pm.pos_metrics.depth, expected_depth,
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

    /// The closed position span seqair's pileup engine and readers take.
    fn span(first: u32, last: u32) -> core::range::RangeInclusive<Pos0> {
        (Pos0::new(first).unwrap()..=Pos0::new(last).unwrap()).into()
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
    /// survive that, through the mate.
    ///
    /// Fixture: mates agree on the anchor base, so the rule drops the right
    /// mate — which is the one carrying the deletion.
    #[test]
    fn a_dropped_mate_does_not_take_the_fragments_indel_with_it() {
        // Segment offsets:  0123456789
        let seq = b"AAAACGTTGG";
        let seg = segment(seq);
        let params = PileupMappingParams { call_indels: true, ..default() };

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
            let store = store_of(&reads, &params);
            let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 9));
            let mut scratch = ColumnScratch::default();
            let mut out = Vec::new();
            while let Some(col) = engine.pileups() {
                if col.pos().as_u64() != 3 {
                    continue;
                }
                let pm =
                    PileupMetrics::from_seqair(&col, seg.clone(), params, &mut scratch).unwrap();
                out = pm
                    .indel_data
                    .map(|d| d.counts.alleles.iter().map(|a| (a.allele.clone(), a.total)).collect())
                    .unwrap_or_default();
            }
            out
        };

        let deletion = IndelAllele::Deletion([Base::C, Base::G].into_iter().collect());
        assert_eq!(
            alleles_at_anchor(&params),
            vec![(deletion.clone(), 1)],
            "the deduplicated column must still report the dropped mate's deletion"
        );

        // With dedup off both mates are their own fragment, and only one of
        // them spans the deletion — so the count is unchanged, not doubled.
        let mut kept = PileupMappingParams { call_indels: true, ..default() };
        kept.variant_calling.keep_overlapping_reads = true;
        assert_eq!(
            alleles_at_anchor(&kept),
            vec![(deletion, 1)],
            "keeping both mates must not turn one fragment's indel into two"
        );
    }

    /// A deletion's REF bases must come from the reference *at the deletion*,
    /// which is `genomic_pos + 1 - segment.start` into the segment's sequence.
    ///
    /// Reading it at `view.qpos()`, a read-local offset, instead of the genomic
    /// `pos` saturates the subtraction to 0 and returns the segment's opening
    /// bases — real-looking sequence from the wrong locus. A segment starting
    /// at 0 cannot see this, because there the clamp and the correct answer
    /// coincide for a read at offset 0.
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
        let params = PileupMappingParams { call_indels: true, ..default() };

        let mut extras = TestExtras::new(&params);
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
                Some(tid0()),
                MateFields::default(),
                &mut extras,
            )
            .unwrap();

        let mut engine = PileupEngine::new(
            store.prepare_for_pileup().input,
            span(START as u32, START as u32 + 9),
        );
        let mut allele = None;
        while let Some(col) = engine.pileups() {
            if col.pos().as_u64() != START + 3 {
                continue;
            }
            let metrics = PileupMetrics::from_seqair(
                &col,
                seg.clone(),
                &params,
                &mut ColumnScratch::default(),
            )
            .unwrap();
            let data = metrics.indel_data.expect("the anchor column carries the deletion");
            let [counts] = data.counts.alleles.as_ref() else {
                panic!("one deletion allele: {:?}", data.counts.alleles)
            };
            allele = Some(counts.allele.clone());
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
        let mut extras = TestExtras::new(&params);

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
                    Some(tid0()),
                    MateFields::default(),
                    &mut extras,
                )
                .unwrap();
            store
        };

        let mut metrics_at = |overhang: u32| -> Option<PileupMetrics> {
            let mut engine =
                PileupEngine::new(build_store().prepare_for_pileup().input, span(0, 5));
            engine.set_soft_clip_overhang(overhang);
            let mut out = None;
            while let Some(col) = engine.pileups() {
                if col.pos() == Pos0::new(2).unwrap() {
                    out = Some(
                        PileupMetrics::from_seqair(
                            &col,
                            seg.clone(),
                            &params,
                            &mut ColumnScratch::default(),
                        )
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
        let mut extras = TestExtras::new(&params);

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
                    Some(tid0()),
                    MateFields::default(),
                    &mut extras,
                )
                .unwrap();
            store
        };

        let mut metrics_at = |overhang: u32| -> Option<PileupMetrics> {
            let mut engine =
                PileupEngine::new(build_store().prepare_for_pileup().input, span(0, 5));
            engine.set_soft_clip_overhang(overhang);
            let mut out = None;
            while let Some(col) = engine.pileups() {
                if col.pos() == Pos0::new(3).unwrap() {
                    out = Some(
                        PileupMetrics::from_seqair(
                            &col,
                            seg.clone(),
                            &params,
                            &mut ColumnScratch::default(),
                        )
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
        let mut extras = TestExtras::new(&params);

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
                Some(tid0()),
                MateFields::default(),
                &mut extras,
            )
            .unwrap();

        let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 5));
        engine.set_soft_clip_overhang(1);
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(3).unwrap() {
                let pm = PileupMetrics::from_seqair(
                    &col,
                    seg.clone(),
                    &params,
                    &mut ColumnScratch::default(),
                )
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
        let mut extras = TestExtras::new(&params);

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
                Some(tid0()),
                MateFields::default(),
                &mut extras,
            )
            .unwrap();

        let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 5));
        engine.set_soft_clip_overhang(1);
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(2).unwrap() {
                let pm = PileupMetrics::from_seqair(
                    &col,
                    seg.clone(),
                    &params,
                    &mut ColumnScratch::default(),
                )
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
        let mut extras = TestExtras::new(&params);

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
                Some(tid0()),
                MateFields::default(),
                &mut extras,
            )
            .unwrap();

        let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 5));
        engine.set_soft_clip_overhang(1);
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(2).unwrap() {
                let pm = PileupMetrics::from_seqair(
                    &col,
                    seg.clone(),
                    &params,
                    &mut ColumnScratch::default(),
                )
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
        let mut extras = TestExtras::new(&params);

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
                Some(tid0()),
                mate_at(3),
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
                Some(tid0()),
                mate_at(2),
                &mut extras,
            )
            .unwrap();

        // `Readers::pileup` links mates after fetching; this test drives the
        // engine directly, so it links by hand.
        let stats = store.link_mates();
        assert_eq!(stats.pairs, 1, "the fixture's mates must link");
        let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 5));
        engine.set_soft_clip_overhang(1);

        let mut checked = false;
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(2).unwrap() {
                // The engine presents both the aligned mate and the rescued clip.
                assert_eq!(col.depth(), 2, "both mates present at the CpG-C before dedup");
                let pm = PileupMetrics::from_seqair(
                    &col,
                    seg.clone(),
                    &params,
                    &mut ColumnScratch::default(),
                )
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
        let mut extras = TestExtras::new(&params);

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
                Some(tid0()),
                mate_at(3),
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
                Some(tid0()),
                mate_at(2),
                &mut extras,
            )
            .unwrap();

        // `Readers::pileup` links mates after fetching; this test drives the
        // engine directly, so it links by hand.
        let stats = store.link_mates();
        assert_eq!(stats.pairs, 1, "the fixture's mates must link");
        let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 5));
        engine.set_soft_clip_overhang(1);

        let mut checked = false;
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(2).unwrap() {
                let pm = PileupMetrics::from_seqair(
                    &col,
                    seg.clone(),
                    &params,
                    &mut ColumnScratch::default(),
                )
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

    // ── differential: mate links vs. pairing by read name ───────────────────

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
        let store = store_of(&reads, &params);
        let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 19));
        engine.set_soft_clip_overhang(1);
        let mut scratch = ColumnScratch::default();
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

    /// `--max-coverage` counts reads that pass the record filter. Two MAPQ 0
    /// reads arrive first at the column; with a cap of 2 they must not take
    /// the two slots, or the column comes out empty.
    #[test]
    fn coverage_cap_counts_filtered_reads_only() {
        let reads =
            [(b"junk1", 0u8), (b"junk2", 0), (b"good1", 60), (b"good2", 60), (b"good3", 60)];
        let params = dedup_params();
        let filter = RastairRecordFilter {
            read_flags: params.read_flags.clone(),
            min_mapq: params.quality.min_mapq,
            ..Default::default()
        };
        let mut capped = DepthCap::with_cap(filter, NonZeroU32::new(2));
        let mut store = RecordStore::<RastairReadExtras>::new();
        for (name, mapq) in reads {
            store
                .push_fields(
                    Pos0::new(2).unwrap(),
                    Pos0::new(9).unwrap(),
                    BamFlags::from(OT_FIRST),
                    mapq,
                    8,
                    0,
                    name,
                    &[CigarOp::new(CigarOpType::Match, 8)],
                    &[Base::A; 8],
                    &[40u8; 8],
                    &[],
                    Some(tid0()),
                    mate_at(500),
                    &mut capped,
                )
                .unwrap();
        }
        let seg = segment(REF);
        let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 9));
        let mut scratch = ColumnScratch::default();
        let mut depth = None;
        while let Some(col) = engine.pileups() {
            if col.pos().as_u64() != 5 {
                continue;
            }
            let pm = PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut scratch).unwrap();
            depth = Some(pm.pos_metrics.depth);
        }
        assert_eq!(depth, Some(2), "the cap applies to reads that pass the MAPQ filter");
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
        let store = store_of(&reads, &params);
        let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 19));
        let mut scratch = ColumnScratch::default();
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

    // ── phase observations ──────────────────────────────────────────────────

    /// A read whose bases are exactly the reference over `[pos, pos + len)`, so
    /// that `with_base_at` introduces the only mismatch in the column.
    fn ref_read(qname: &[u8], pos: u32, len: usize, flags: u16) -> TestRead {
        let start = pos as usize;
        let mut read = TestRead::matching(qname, pos, len, Base::A, flags);
        read.bases = REF
            .get(start..start + len)
            .expect("read lies inside the reference")
            .iter()
            .map(|&b| Base::from(b))
            .collect();
        read
    }

    fn phasing_params() -> PileupMappingParams {
        PileupMappingParams { phase: true, ..dedup_params() }
    }

    /// Build the column at `pos`, handing back the scratch too so a test can see
    /// the counters that live for a whole region.
    fn phase_column(
        reads: &[TestRead],
        params: &PileupMappingParams,
        pos: u32,
    ) -> (PileupMetrics, ColumnScratch) {
        let seg = segment(REF);
        let store = store_of(reads, params);
        let last = u32::try_from(REF.len() - 1).expect("a short reference");
        let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, last));
        let mut scratch = ColumnScratch::default();
        let mut metrics = None;
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(pos).unwrap() {
                metrics = Some(
                    PileupMetrics::from_seqair(&col, seg.clone(), params, &mut scratch).unwrap(),
                );
            }
        }
        (metrics.expect("column at the requested position"), scratch)
    }

    fn observations(pm: &PileupMetrics) -> Option<&[crate::call::phasing::PhaseObservation]> {
        match pm.phase.as_deref() {
            Some(ColumnPhase::Observed(observations)) => Some(observations),
            _ => None,
        }
    }

    /// Reference position 8 is an A outside any CpG.
    const REF_A: u32 = 8;

    #[test]
    fn no_column_keeps_observations_without_the_flag() {
        let reads = vec![
            ref_read(b"g0", 4, 8, OT_FIRST).with_base_at(4, Base::G),
            ref_read(b"g1", 4, 8, OB_FIRST).with_base_at(4, Base::G),
            ref_read(b"r0", 4, 8, OT_FIRST),
        ];
        let (pm, _) = phase_column(&reads, &dedup_params(), REF_A);
        assert!(pm.alt(Base::G).is_some(), "the column really does carry alts");
        assert!(pm.phase.is_none(), "observations must cost nothing without --phase");
    }

    #[test]
    fn a_het_column_keeps_one_observation_per_read() {
        let reads = vec![
            ref_read(b"g0", 4, 8, OT_FIRST).with_base_at(4, Base::G),
            ref_read(b"g1", 4, 8, OB_FIRST).with_base_at(4, Base::G),
            ref_read(b"r0", 4, 8, OT_FIRST),
        ];
        let (pm, _) = phase_column(&reads, &phasing_params(), REF_A);
        let obs = observations(&pm).expect("het column is a phase candidate");
        assert_eq!(obs.len(), 3, "every kept read contributes, reference reads included");
        assert_eq!(obs.iter().filter(|o| o.base == Base::G).count(), 2);

        let mut fragments: Vec<_> = obs.iter().map(|o| o.fragment).collect();
        fragments.sort_unstable();
        fragments.dedup();
        assert_eq!(fragments.len(), 3, "three templates, three fragment ids");
    }

    /// Both mates of a template carry the same fragment id, so a het one mate
    /// covers can be linked to a het only the other mate reaches.
    #[test]
    fn both_mates_of_a_pair_share_a_fragment_id() {
        let reads = vec![
            ref_read(b"p0", 0, 8, OT_FIRST).with_base_at(4, Base::G),
            ref_read(b"p1", 0, 8, OB_FIRST).with_base_at(4, Base::G),
            ref_read(b"p0", 8, 8, OT_SECOND).with_base_at(4, Base::G),
            ref_read(b"p1", 8, 8, OB_SECOND).with_base_at(4, Base::G),
        ];
        let params = phasing_params();

        let fragments_at = |pos: u32| {
            let (pm, _) = phase_column(&reads, &params, pos);
            let mut ids: Vec<_> = observations(&pm)
                .unwrap_or_else(|| panic!("column {pos} is a phase candidate"))
                .iter()
                .map(|o| o.fragment)
                .collect();
            ids.sort_unstable();
            ids
        };

        let left = fragments_at(4);
        assert_eq!(left.len(), 2);
        assert_eq!(left, fragments_at(12), "the same two fragments span both columns");
    }

    #[test]
    fn phase_quality_thresholds_are_stricter_than_the_calling_ones() {
        let reads = vec![
            ref_read(b"g0", 4, 8, OT_FIRST).with_base_at(4, Base::G),
            ref_read(b"g1", 4, 8, OT_FIRST).with_base_at(4, Base::G),
            ref_read(b"lowbaseq", 4, 8, OB_FIRST).with_base_at(4, Base::G).with_qual(10),
            ref_read(b"lowmapq", 4, 8, OB_FIRST).with_base_at(4, Base::G).with_mapq(5),
        ];
        let params = phasing_params();
        assert!(
            params.quality.min_mapq <= 5 && params.quality.min_baseq <= 10,
            "the calling filters must accept both reads, or this proves nothing"
        );

        let (pm, _) = phase_column(&reads, &params, REF_A);
        assert_eq!(pm.pos_metrics.depth, 4);
        let obs = observations(&pm).expect("candidate");
        assert_eq!(obs.len(), 2, "the low-quality reads contribute depth but not phase");
    }

    /// A read without a name has no template identity; grouping such reads
    /// would link every het in the region, so they are dropped and flagged.
    #[test]
    fn reads_without_a_name_contribute_no_observation() {
        let reads = vec![
            ref_read(b"g0", 4, 8, OT_FIRST).with_base_at(4, Base::G),
            ref_read(b"g1", 4, 8, OT_FIRST).with_base_at(4, Base::G),
            ref_read(b"", 4, 8, OB_FIRST).with_base_at(4, Base::G),
        ];
        let (pm, scratch) = phase_column(&reads, &phasing_params(), REF_A);
        assert_eq!(pm.pos_metrics.depth, 3, "the nameless read still counts towards depth");
        assert_eq!(observations(&pm).map(<[_]>::len), Some(2));
        assert!(scratch.saw_nameless_read);
    }

    /// Overlapping mates are collapsed by the dedup first, so one fragment
    /// cannot observe a column twice.
    #[test]
    fn overlapping_mates_observe_a_column_once() {
        let reads = vec![
            ref_read(b"pair", 4, 8, OT_FIRST).with_base_at(4, Base::G),
            ref_read(b"pair", 6, 8, OT_SECOND).with_base_at(2, Base::G),
            ref_read(b"other", 4, 8, OB_FIRST).with_base_at(4, Base::G),
        ];
        let (pm, _) = phase_column(&reads, &phasing_params(), REF_A);
        let obs = observations(&pm).expect("candidate");
        assert_eq!(obs.len(), 2, "the pair speaks once");
        assert_ne!(obs.first().map(|o| o.fragment), obs.get(1).map(|o| o.fragment));
    }

    /// Non-repetitive reference for the pair-HMM tests, placed off zero.
    const HAPLOTYPE_START: u32 = 1000;
    const HAPLOTYPE_SEQ: &[u8] = b"ACGTTAGCATCGGATCCATGACTGATCGGCTAGCTTACGGATCCTAGGCTAAGCTTGACT";
    /// The column the pair-HMM tests' insertions sit at.
    const HAPLOTYPE_ANCHOR: u32 = HAPLOTYPE_START + 29;

    /// `<len>M` at `pos` copying [`HAPLOTYPE_SEQ`], its own fragment.
    fn haplotype_read(qname: &[u8], pos: u32, len: usize, flags: u16) -> TestRead {
        let from = (pos - HAPLOTYPE_START) as usize;
        let mut read = TestRead::matching(qname, pos, len, Base::A, flags);
        read.bases = HAPLOTYPE_SEQ[from..from + len].iter().map(|&b| Base::from(b)).collect();
        read
    }

    /// A read over offsets 5..55 of [`HAPLOTYPE_SEQ`] with `inserted` after the
    /// anchor, 25 aligned bases either side.
    fn inserting_read(qname: &[u8], flags: u16, inserted: &[u8]) -> TestRead {
        let bases: Vec<Base> = HAPLOTYPE_SEQ[5..30]
            .iter()
            .chain(inserted)
            .chain(&HAPLOTYPE_SEQ[30..55])
            .map(|&b| Base::from(b))
            .collect();
        TestRead {
            qname: qname.to_vec(),
            pos: HAPLOTYPE_START + 5,
            flags,
            quals: vec![40; bases.len()],
            bases,
            mapq: 60,
            cigar: vec![
                CigarOp::new(CigarOpType::Match, 25),
                CigarOp::new(CigarOpType::Insertion, u32::try_from(inserted.len()).unwrap()),
                CigarOp::new(CigarOpType::Match, 25),
            ],
            mate_pos: -1,
        }
    }

    /// The column at `pos` over [`HAPLOTYPE_SEQ`].
    fn haplotype_column(
        reads: &[TestRead],
        params: &PileupMappingParams,
        pos: u32,
    ) -> PileupMetrics {
        let seg = segment_at(u64::from(HAPLOTYPE_START), HAPLOTYPE_SEQ);
        let store = store_of(reads, params);
        let last = HAPLOTYPE_START + u32::try_from(HAPLOTYPE_SEQ.len()).unwrap() - 1;
        let mut engine =
            PileupEngine::new(store.prepare_for_pileup().input, span(HAPLOTYPE_START, last));
        let mut scratch = ColumnScratch::default();
        let mut metrics = None;
        while let Some(col) = engine.pileups() {
            if col.pos().as_u64() == u64::from(pos) {
                metrics = Some(
                    PileupMetrics::from_seqair(&col, seg.clone(), params, &mut scratch).unwrap(),
                );
            }
        }
        metrics.expect("the column")
    }

    fn fragment_of(qname: &[u8]) -> Option<FragmentId> {
        FragmentId::new(seqair::bam::record_store::qname_hash(qname))
    }

    /// Three reads insert `GA`, three do not; all six are scored, each
    /// under its own fragment.
    #[test]
    fn an_indel_column_is_scored_with_fragment_ids() {
        let mut reads: Vec<TestRead> =
            (0..3).map(|i| inserting_read(format!("ins{i}").as_bytes(), 0, b"GA")).collect();
        reads
            .extend((0..3).map(|i| {
                haplotype_read(format!("ref{i}").as_bytes(), HAPLOTYPE_START + 5, 50, 16)
            }));
        let scoring = PileupMappingParams { call_indels: true, ..default() };

        let pm = haplotype_column(&reads, &scoring, HAPLOTYPE_ANCHOR);
        let data = pm.indel_data.as_deref().expect("an insertion column");
        let scores = &data.scored.as_ref().expect("scores").haplotypes;
        let inserted = IndelAllele::Insertion([Base::G, Base::A].into_iter().collect());
        assert_eq!(scores.cluster().events().len(), 1);
        assert_eq!(scores.cluster().events()[0].allele, inserted.into());
        let mut ids: Vec<_> = scores.reads.iter().map(|r| r.fragment).collect();
        ids.sort_unstable();
        let mut expected: Vec<_> = reads.iter().map(|r| fragment_of(&r.qname)).collect();
        expected.sort_unstable();
        assert_eq!(ids, expected);
        assert!(ids.iter().all(Option::is_some));
        let (low, high) = scores.verdict.best.haplotypes();
        assert_eq!((low.index(), high.index()), (0, 1));
    }

    /// An allele below `--min-indel-ao` is not a candidate, so a column with
    /// no other has nothing to score.
    #[test]
    fn an_allele_below_the_minimum_support_is_not_scored() {
        let mut reads = vec![inserting_read(b"ins", 0, b"GA")];
        reads
            .extend((0..3).map(|i| {
                haplotype_read(format!("ref{i}").as_bytes(), HAPLOTYPE_START + 5, 50, 16)
            }));
        let params = PileupMappingParams { call_indels: true, ..default() };
        let pm = haplotype_column(&reads, &params, HAPLOTYPE_ANCHOR);
        assert!(pm.indel_data.expect("an insertion column").scored.is_none());
    }

    fn deleted(bases: &[u8]) -> IndelAllele {
        IndelAllele::Deletion(bases.iter().map(|&b| Base::from(b)).collect())
    }

    /// A column at `pos` whose reads carry `alleles`, each on `total` reads.
    fn indel_column(pos: u32, alleles: &[(IndelAllele, u32)]) -> PileupMetrics {
        let mut column = phasing::fixtures::column(pos, Base::A, &[], None, &[]);
        let alleles = alleles
            .iter()
            .map(|(allele, total)| IndelAlleleCounts { allele: allele.clone(), total: *total })
            .collect();
        let counts = IndelCounts { alleles };
        column.indel_data = Some(Box::new(IndelData { counts, scored: None }));
        column
    }

    fn events_of(columns: &[PileupMetrics]) -> Vec<(u32, IndelAllele)> {
        cluster_candidates(columns.iter(), 2)
            .into_iter()
            .filter_map(|e| Some((e.anchor, e.allele.as_indel()?.clone())))
            .collect()
    }

    /// Most reads first, then the earlier anchor, then insertions before
    /// deletions; an allele below `--min-indel-ao` is no event, and the
    /// fourth is left out.
    #[test]
    fn a_cluster_takes_its_best_supported_events_in_order() {
        let insertion = IndelAllele::Insertion([Base::A].into_iter().collect());
        let columns = [
            indel_column(100, &[(deleted(b"T"), 5), (insertion.clone(), 5), (deleted(b"TT"), 1)]),
            indel_column(105, &[(deleted(b"G"), 5)]),
            indel_column(110, &[(deleted(b"TT"), 8)]),
        ];
        assert_eq!(
            events_of(&columns),
            [(110, deleted(b"TT")), (100, insertion), (100, deleted(b"T"))]
        );
    }

    /// A third allele at one anchor could not share a pair with the other
    /// two, so the next anchor's event takes its place.
    #[test]
    fn a_cluster_takes_at_most_two_events_per_anchor() {
        let columns = [
            indel_column(100, &[(deleted(b"T"), 9), (deleted(b"TT"), 8), (deleted(b"TTT"), 7)]),
            indel_column(120, &[(deleted(b"G"), 3)]),
        ];
        assert_eq!(
            events_of(&columns),
            [(100, deleted(b"T")), (100, deleted(b"TT")), (120, deleted(b"G"))]
        );
    }

    /// A column with `alt` on `ot`/`ob` reads and the reference on
    /// `ref_ot`/`ref_ob`, depth their sum.
    fn snv_column(
        pos: u32,
        reference: Base,
        alt: Base,
        [ot, ob]: [u32; 2],
        [ref_ot, ref_ob]: [u32; 2],
    ) -> PileupMetrics {
        let mut column = mismatch_column(pos, reference, alt, ot, ob);
        column.ref_metrics.strand_count.ot = ref_ot;
        column.ref_metrics.strand_count.ob = ref_ob;
        column.ref_metrics.depth = ref_ot + ref_ob;
        column.pos_metrics.depth = ot + ob + ref_ot + ref_ob;
        column
    }

    fn admitted(column: &PileupMetrics) -> Vec<Base> {
        snv_events(column)
            .iter()
            .filter_map(|e| match e.allele {
                Allele::Base(b) => Some(b),
                Allele::Indel(_) => None,
            })
            .collect()
    }

    /// Two reads and a het the informative strand allows: a `T` over a `C`
    /// counts on OB only, where TAPS converts nothing, so six OT `T`s are a
    /// methylated `CpG` and no event, while three OB `T`s among twelve are.
    #[test]
    fn an_snv_event_needs_two_reads_and_a_possible_het() {
        let (c, t, a) = (Base::C, Base::T, Base::A);
        assert_eq!(admitted(&snv_column(10, c, t, [6, 0], [0, 12])), [] as [Base; 0]);
        assert_eq!(admitted(&snv_column(10, c, t, [0, 3], [6, 9])), [t]);
        assert_eq!(admitted(&snv_column(10, c, a, [1, 0], [6, 6])), [] as [Base; 0], "one read");
        assert_eq!(admitted(&snv_column(10, c, a, [1, 1], [6, 6])), [a]);
        // 3 of 70 is below 5 % of the column.
        assert_eq!(admitted(&snv_column(10, c, a, [2, 1], [33, 34])), [] as [Base; 0]);
        let mut unknown = snv_column(10, c, a, [3, 3], [3, 3]);
        unknown.reference_base = Base::Unknown;
        assert!(snv_events(&unknown).is_empty());
    }

    fn snv_at(pos: u32) -> SmallVec<Event, 2> {
        [Event { anchor: pos, allele: Allele::Base(Base::G) }].into_iter().collect()
    }

    /// Indel columns chain while each is within `CLUSTER_GAP` of the last.
    #[test]
    fn indel_columns_chain_within_the_cluster_gap() {
        let columns = [
            indel_column(100, &[(deleted(b"T"), 5)]),
            indel_column(115, &[(deleted(b"G"), 4)]),
            indel_column(135, &[(deleted(b"G"), 4)]),
            indel_column(200, &[(deleted(b"TT"), 3)]),
        ];
        let chains: Vec<&[usize]> = vec![&[0, 1, 2], &[3]];
        assert_eq!(cluster_chains(&columns, &[0, 1, 2, 3]), chains);
    }

    /// An SNV near an indel chain is scored with the chain's indel events as
    /// context; the rest chain among themselves, cut where a fourth event
    /// would not fit, and one that would not fit beside the indels goes with
    /// the rest.
    #[test]
    fn snvs_join_a_nearby_indel_chain_or_chain_on_their_own() {
        let columns: Vec<PileupMetrics> = [90, 100, 110, 125, 135, 140, 145, 150, 400]
            .into_iter()
            .map(|pos| {
                if pos == 100 {
                    indel_column(100, &[(deleted(b"T"), 5), (deleted(b"TT"), 4)])
                } else {
                    mismatch_column(pos, Base::A, Base::G, 2, 2)
                }
            })
            .collect();
        let snvs: Vec<_> = [0usize, 2, 3, 4, 5, 6, 7, 8]
            .into_iter()
            .map(|i| (i, snv_at(columns[i].pos)))
            .collect();
        let clusters = snv_clusters(&columns, &[1], &snvs, 2);
        let summary: Vec<(Vec<u32>, usize)> = clusters
            .iter()
            .map(|c| (c.columns.iter().map(|&i| columns[i].pos).collect(), c.events.len()))
            .collect();
        assert_eq!(
            summary,
            [
                (vec![90, 100], 3),
                (vec![110, 125, 135], 3),
                (vec![140, 145, 150], 3),
                (vec![400], 1),
            ]
        );
        assert!(clusters[0].events[0].allele.as_indel().is_some(), "indel events rank first");
        assert!(clusters[0].events[1].allele.as_indel().is_some());
    }

    /// A column with `ot` and `ob` reads of `alt` over `reference`.
    fn mismatch_column(pos: u32, reference: Base, alt: Base, ot: u32, ob: u32) -> PileupMetrics {
        let mut column = phasing::fixtures::column(
            pos,
            reference,
            &[(alt, crate::metrics::AltCall::Uncalled)],
            None,
            &[],
        );
        for alt in &mut column.alts {
            alt.metrics.strand_count.ot = ot;
            alt.metrics.strand_count.ob = ob;
            alt.metrics.depth = ot + ob;
        }
        column
    }

    /// Five reads of a base TAPS cannot explain mask a column: a `T` over a
    /// `C` only on the bottom strand, an `A` over a `G` only on the top, any
    /// other base on either; the conversions TAPS makes do not.
    #[test]
    fn the_window_mask_takes_mismatches_taps_does_not_explain() {
        use Base::{A, C, G, T};
        let columns = [
            mismatch_column(10, C, T, 0, 5),
            mismatch_column(11, C, T, 5, 0),
            mismatch_column(12, G, A, 5, 0),
            mismatch_column(13, G, A, 0, 5),
            mismatch_column(14, A, C, 3, 2),
            mismatch_column(15, A, C, 2, 2),
        ];
        assert_eq!(WindowMask::new(&columns).0, [10, 12, 14]);
    }

    /// The scores of a 1 bp deletion at [`HAPLOTYPE_ANCHOR`] under `mask`.
    fn deletion_scores(reads: &[TestRead], mask: WindowMask) -> Option<HaplotypeScores> {
        let seg = segment_at(u64::from(HAPLOTYPE_START), HAPLOTYPE_SEQ);
        let params = PileupMappingParams { call_indels: true, ..default() };
        let store = store_of(reads, &params);
        let last = HAPLOTYPE_START + u32::try_from(HAPLOTYPE_SEQ.len()).unwrap() - 1;
        let mut engine =
            PileupEngine::new(store.prepare_for_pileup().input, span(HAPLOTYPE_START, last));
        let events: SmallVec<Event, MAX_EVENTS> =
            [Event { anchor: HAPLOTYPE_ANCHOR, allele: deleted(b"T").into() }]
                .into_iter()
                .collect();
        while let Some(col) = engine.pileups() {
            let here = col.pos();
            if here.as_u64() == u64::from(HAPLOTYPE_ANCHOR) {
                let spanning = col.records_overlapping(here, here);
                let scratch = &mut HaplotypeScratch::default();
                let context =
                    ClusterContext { store: col.store(), segment: &seg, params: &params, mask };
                return haplotype_scores(&context, spanning, events, scratch);
            }
        }
        None
    }

    /// Reads ending on an SNV after the anchor whose base is the one after it:
    /// the reference charges them a mismatch, the deletion of the SNV's base
    /// explains them. Masked, the SNV's position matches any base, and the
    /// deletion has nothing left to explain.
    #[test]
    fn a_masked_position_matches_any_base() {
        let snv = HAPLOTYPE_ANCHOR + 1;
        let offset = (snv - HAPLOTYPE_START) as usize;
        assert_eq!(HAPLOTYPE_SEQ[offset + 1], b'A', "the base after the SNV");
        let mut reads: Vec<TestRead> = (0..6)
            .map(|i| {
                let start = HAPLOTYPE_START + 5 + i;
                let len = (snv - start + 1) as usize;
                haplotype_read(format!("snv{i}").as_bytes(), start, len, 0)
                    .with_base_at(len - 1, Base::A)
            })
            .collect();
        reads
            .extend((0..4).map(|i| {
                haplotype_read(format!("ref{i}").as_bytes(), HAPLOTYPE_START + 5, 50, 16)
            }));
        let event = pair_hmm::EventIdx::new(0).unwrap();

        let unmasked = deletion_scores(&reads, WindowMask::default()).expect("scores");
        assert_eq!(unmasked.verdict.copies_of(event), 1, "{:?}", unmasked.verdict);
        let masked = deletion_scores(&reads, WindowMask(vec![snv])).expect("scores");
        assert_eq!(masked.verdict.copies_of(event), 0, "{:?}", masked.verdict);
    }

    /// A read at `pos` of [`HAPLOTYPE_SEQ`]'s bases, clipped `lead` bases
    /// before and `trail` after `aligned` aligned ones.
    fn clipped_read(qname: &[u8], pos: u32, lead: u32, aligned: u32, trail: u32) -> TestRead {
        // The clipped bases are the reference's where it has them; which bases
        // a clip off the segment carries does not matter to its placement.
        let from =
            usize::try_from(pos.saturating_sub(HAPLOTYPE_START).saturating_sub(lead)).unwrap();
        let len = usize::try_from(lead + aligned + trail).unwrap();
        let mut cigar = vec![CigarOp::new(CigarOpType::Match, aligned)];
        if lead > 0 {
            cigar.insert(0, CigarOp::new(CigarOpType::SoftClip, lead));
        }
        if trail > 0 {
            cigar.push(CigarOp::new(CigarOpType::SoftClip, trail));
        }
        let bases: Vec<Base> =
            HAPLOTYPE_SEQ.iter().cycle().skip(from).take(len).map(|&b| Base::from(b)).collect();
        TestRead {
            cigar,
            quals: vec![40; bases.len()],
            bases,
            ..TestRead::matching(qname, pos, len, Base::A, 0)
        }
    }

    /// Each read the pair-HMM would score at `pos` of [`HAPLOTYPE_SEQ`] placed
    /// at `start`: its placed span, the stored bases it scores, and its name.
    fn placed_reads_at(
        reads: &[TestRead],
        start: u32,
        pos: u32,
    ) -> Vec<(u64, u64, Range<usize>, Vec<u8>)> {
        let seg = segment_at(u64::from(start), HAPLOTYPE_SEQ);
        let params = PileupMappingParams::default();
        let store = store_of(reads, &params);
        let last = start + u32::try_from(HAPLOTYPE_SEQ.len()).unwrap() - 1;
        let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(start, last));
        let mut placed = Vec::new();
        while let Some(col) = engine.pileups() {
            if col.pos().as_u64() == u64::from(pos) {
                let here = col.pos();
                placed =
                    spanning_reads(col.store(), col.records_overlapping(here, here), &seg, &params)
                        .iter()
                        .map(|r| (r.start, r.end, r.bases.clone(), r.rec.qname().to_vec()))
                        .collect();
            }
        }
        placed
    }

    /// A soft clip the segment's reference reaches is scored with the read,
    /// on either end; the last reference base is inside, `end` being inclusive.
    #[test]
    fn a_clip_that_fits_the_segment_is_placed_with_the_read() {
        let lead = clipped_read(b"lead", HAPLOTYPE_START + 10, 3, 20, 0);
        assert_eq!(
            placed_reads_at(&[lead], HAPLOTYPE_START, HAPLOTYPE_START + 15),
            [(1007, 1029, 0..23, b"lead".to_vec())]
        );
        // Aligned to 1037..=1056, the trailing clip ends on the last base, 1059.
        let trail = clipped_read(b"trail", HAPLOTYPE_START + 37, 0, 20, 3);
        assert_eq!(
            placed_reads_at(&[trail], HAPLOTYPE_START, HAPLOTYPE_START + 45),
            [(1037, 1059, 0..23, b"trail".to_vec())]
        );
    }

    /// A clip that would run off the segment is dropped, and the read is
    /// scored on its aligned bases alone.
    #[test]
    fn a_clip_outside_the_segment_leaves_the_aligned_bases() {
        let lead = clipped_read(b"lead", HAPLOTYPE_START + 1, 3, 20, 0);
        assert_eq!(
            placed_reads_at(&[lead], HAPLOTYPE_START, HAPLOTYPE_START + 5),
            [(1001, 1020, 3..23, b"lead".to_vec())]
        );
        let trail = clipped_read(b"trail", HAPLOTYPE_START + 38, 0, 20, 3);
        assert_eq!(
            placed_reads_at(&[trail], HAPLOTYPE_START, HAPLOTYPE_START + 45),
            [(1038, 1057, 0..20, b"trail".to_vec())]
        );
    }

    /// At the contig start a clip can reach before position 0; the read is
    /// then placed on its aligned bases, not shifted onto 0.
    #[test]
    fn a_clip_before_the_contig_start_leaves_the_aligned_bases() {
        let lead = clipped_read(b"lead", 1, 3, 20, 0);
        assert_eq!(placed_reads_at(&[lead], 0, 5), [(1, 20, 3..23, b"lead".to_vec())]);
    }

    /// One read per fragment, whichever mate reaches further, and a read
    /// without a mate is its own fragment.
    #[test]
    fn a_fragment_is_scored_once_from_the_mate_reaching_further() {
        let mate = |qname: &[u8], pos: u32, len: u32, flags: u16| TestRead {
            flags,
            ..clipped_read(qname, pos, 0, len, 0)
        };
        let reads = [
            mate(b"right", HAPLOTYPE_START, 20, 99),
            mate(b"right", HAPLOTYPE_START + 10, 30, 147),
            mate(b"left", HAPLOTYPE_START + 2, 40, 99),
            mate(b"left", HAPLOTYPE_START + 12, 20, 147),
            mate(b"solo", HAPLOTYPE_START + 5, 20, 0),
        ];
        let mut placed = placed_reads_at(&reads, HAPLOTYPE_START, HAPLOTYPE_START + 15);
        placed.sort_by(|a, b| a.3.cmp(&b.3));
        assert_eq!(
            placed,
            [
                (1002, 1041, 0..40, b"left".to_vec()),
                (1010, 1039, 0..30, b"right".to_vec()),
                (1005, 1024, 0..20, b"solo".to_vec()),
            ]
        );
    }
}
