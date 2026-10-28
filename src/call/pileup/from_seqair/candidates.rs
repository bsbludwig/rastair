//! Which alleles a cluster scores, and which columns form one.
//!
//! Indel columns are not scored while they are built. After a sub-segment's
//! column loop, kept indel columns within `CLUSTER_GAP` of the previous one
//! form a cluster (chained, so it can be longer).
//!
//! - [`cluster_candidates`]: a cluster's events are its columns' alleles
//!   [`CandidateFloor`] admits (`--min-indel-ao`), ranked by read count, then
//!   anchor, then `allele_order` (insertions first), as many as
//!   `Cluster::admits`: at most `MAX_EVENTS`, and two per anchor, since a
//!   third could not share a pair with the other two and would take a
//!   neighbour's slot without ever being called. That ranking is record
//!   order: `Cluster::events`, `EventIdx`, the calls and the `ALT` column keep
//!   it. An allele not ranked in is never called.
//! - [`cascade_events`] admits an SNV alt to the SNV engine (see
//!   `call::process::pileups`): of the alts the first-stage forest scores at
//!   [`CASCADE_MIN_SCORE`] or more, the two it scores highest. An SNV within
//!   `CLUSTER_GAP` of an indel chain joins the chain's events as context ([`snv_clusters`]); such a cluster
//!   ([`ClusterPass::Snvs`]) fills `snv_scores` and never writes an indel
//!   record.
//! - Each SNV column is scored with its nearest neighbours ([`snv_clusters`]),
//!   so its record does not depend on which side they lie. A record's `GL`
//!   holds the cluster's other events at their calls: an event scored beside
//!   a linked neighbour has `0/0` far less likely than one scored alone
//!   (PL 255 against about 90 at 32x), and a chain cut left to right would
//!   give that to the events on one side of the cut only. Columns with one
//!   neighbourhood share a cluster, so only chains of more events than a
//!   cluster holds cost more.

use crate::call::pileup::indels::{IndelAllele, IndelAlleleCounts};
use crate::{
    call::pair_hmm::{Allele, Cluster, Event, MAX_EVENTS},
    metrics::{Alt, PileupMetrics},
};
use seqair_types::{Pos0, SmallVec};
use std::cmp::Reverse;

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
        .chunk_by(|&a, &b| {
            pos(a).zip(pos(b)).is_some_and(|(a, b)| b.as_u32().abs_diff(a.as_u32()) <= CLUSTER_GAP)
        })
        .collect()
}

/// An SNV cluster to score: its columns (indices into the pileups,
/// ascending), whose spanning reads it is scored on, the SNV columns that
/// take its scores, and its events in record order, the indel events first.
#[derive(Debug, Default, PartialEq)]
pub(crate) struct SnvCluster {
    pub(crate) columns: Vec<usize>,
    pub(crate) scored: Vec<usize>,
    pub(crate) events: SmallVec<Event, MAX_EVENTS>,
}

/// A column the cascade hands to the engine, by its index, with its events.
#[derive(Debug)]
pub(crate) struct SnvColumn {
    pub(crate) column: usize,
    pub(crate) events: SmallVec<Event, 2>,
}

/// The SNV clusters of a sub-segment, one per SNV column, those of equal
/// events merged. An SNV column within [`CLUSTER_GAP`] of one of the
/// `indel_chains` (each with the events the indel pass scored) whose events
/// it admits takes those as context, the nearest such chain; the other SNV
/// columns chain among themselves. Then each column takes its nearest
/// neighbours among the SNV columns of its chain or indel context while
/// they fit, both of two at one distance or neither, so a column's cluster
/// does not depend on which side its neighbours lie. `snvs` are SNV columns
/// with their events, ascending; a column's events stay together.
pub(crate) fn snv_clusters(
    pileups: &[PileupMetrics],
    indel_chains: Vec<SnvCluster>,
    snvs: &[SnvColumn],
) -> Vec<SnvCluster> {
    let pos = |i: usize| pileups.get(i).map(|p| p.pos.as_u32());
    let distance = |a: usize, b: usize| pos(a).zip(pos(b)).map(|(a, b)| a.abs_diff(b));
    let near = |a: usize, b: usize| distance(a, b).is_some_and(|d| d <= CLUSTER_GAP);
    let fits = |events: &mut SmallVec<Event, MAX_EVENTS>, snv: &SmallVec<Event, 2>| {
        snv.iter().all(|event| {
            let admitted = Cluster::admits(events, event);
            events.push(event.clone());
            admitted
        })
    };
    let context: Vec<Option<&SnvCluster>> = snvs
        .iter()
        .map(|snv| {
            indel_chains
                .iter()
                .filter(|chain| {
                    let (Some(&first), Some(&last)) = (chain.columns.first(), chain.columns.last())
                    else {
                        return false;
                    };
                    (near(first, snv.column)
                        || near(snv.column, last)
                        || (first..=last).contains(&snv.column))
                        && fits(&mut chain.events.clone(), &snv.events)
                })
                .min_by_key(|chain| {
                    chain.columns.iter().filter_map(|&c| distance(c, snv.column)).min()
                })
        })
        .collect();
    // Columns of one group are each other's neighbours: those of an indel
    // context by its chain, the others by chaining within `CLUSTER_GAP`.
    let mut groups: Vec<(Option<usize>, usize)> = Vec::with_capacity(snvs.len());
    let mut last_free: Option<usize> = None;
    let mut free_chain = 0;
    for (snv, context) in snvs.iter().zip(&context) {
        let group = match context {
            Some(chain) => (chain.columns.first().copied(), 0),
            None => {
                if !last_free.is_some_and(|last| near(last, snv.column)) {
                    free_chain += 1;
                }
                last_free = Some(snv.column);
                (None, free_chain)
            }
        };
        groups.push(group);
    }

    let mut clusters: Vec<SnvCluster> = Vec::new();
    for ((owner, context), group) in snvs.iter().zip(&context).zip(&groups) {
        let mut room = context.map(|chain| chain.events.clone()).unwrap_or_default();
        room.extend(owner.events.iter().cloned());
        let mut neighbours: Vec<&SnvColumn> = snvs
            .iter()
            .zip(&groups)
            .filter(|&(snv, g)| g == group && snv.column != owner.column)
            .map(|(snv, _)| snv)
            .collect();
        let from_owner = |snv: &&SnvColumn| distance(snv.column, owner.column);
        neighbours.sort_by_key(from_owner);
        let mut members = vec![owner];
        for tied in neighbours.chunk_by(|a, b| from_owner(a) == from_owner(b)) {
            let mut trial = room.clone();
            if tied.iter().all(|snv| fits(&mut trial, &snv.events)) {
                room = trial;
                members.extend(tied);
            }
        }
        members.sort_by_key(|snv| snv.column);
        let mut events = context.map(|chain| chain.events.clone()).unwrap_or_default();
        events.extend(members.iter().flat_map(|snv| snv.events.iter().cloned()));
        let mut columns = context.map(|chain| chain.columns.clone()).unwrap_or_default();
        columns.extend(members.iter().map(|snv| snv.column));
        columns.sort_unstable();
        columns.dedup();
        match clusters.iter_mut().find(|c| c.events == events && c.columns == columns) {
            Some(same) => same.scored.push(owner.column),
            None => clusters.push(SnvCluster { columns, scored: vec![owner.column], events }),
        }
    }
    clusters
}

/// The lowest first-stage forest score at which the cascade hands an alt to
/// the engine. The bundled forests put most true SNVs at .995-.999, so there
/// is no upper edge; .005 keeps 70-110 % of the gain of scoring every gated
/// alt for a quarter of its clusters on 26x slices, .01 loses some.
pub(crate) const CASCADE_MIN_SCORE: f64 = 0.005;

/// The events the cascade scores: the alts the first stage scored at
/// [`CASCADE_MIN_SCORE`] or more, whatever their reads, the two it scored
/// highest (a diploid holds no more); none at an unknown reference base.
pub(crate) fn cascade_events(column: &PileupMetrics) -> SmallVec<Event, 2> {
    if column.reference_base.known_index().is_none() {
        return SmallVec::new();
    }
    let score = |alt: &Alt| alt.verdict.ml().map_or(0.0, |ml| *ml);
    let mut admitted: SmallVec<&Alt, 4> = column
        .alts
        .iter()
        .filter(|alt| alt.base.known_index().is_some())
        .filter(|alt| score(alt) >= CASCADE_MIN_SCORE)
        .collect();
    admitted.sort_by(|a, b| score(b).total_cmp(&score(a)).then(a.base.cmp(&b.base)));
    admitted
        .into_iter()
        .take(2)
        .map(|alt| Event { anchor: column.pos, allele: Allele::Base(alt.base) })
        .collect()
}

/// Reads a candidate indel allele needs to become a pair-HMM event:
/// `min_reads`, and 5 % of the column's depth. Two reads of a stutter allele are noise at 73x
/// and the only evidence at 14x, so the floor grows with the column: for
/// indels, 5 % drops a quarter of the candidate events and 17 % of the
/// pair-HMM cells at 73x and changes no column of the chr12 scores at either
/// depth; 2 % drops nothing. The alleles count fragments, not reads, unless
/// `--keep-overlapping-reads`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(transparent)]
pub struct CandidateFloor {
    pub min_reads: u32,
}

impl std::fmt::Display for CandidateFloor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.min_reads.fmt(f)
    }
}

impl std::str::FromStr for CandidateFloor {
    type Err = std::num::ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self { min_reads: s.parse()? })
    }
}

impl CandidateFloor {
    pub(crate) fn at(self, depth: u32) -> u32 {
        self.min_reads.max((depth * 5).div_ceil(100))
    }
}

/// The cluster's events: the alleles the most reads carry, of those
/// [`CandidateFloor`] admits, as many as [`Cluster::admits`], in record order. Ties
/// go by the allele itself, so neither the choice nor the order depends on
/// which read the engine yielded first.
pub(crate) fn cluster_candidates<'a>(
    columns: impl Iterator<Item = &'a PileupMetrics>,
    floor: CandidateFloor,
) -> SmallVec<Event, MAX_EVENTS> {
    let mut alleles: SmallVec<(Pos0, &IndelAlleleCounts), 8> = columns
        .filter_map(|column| {
            let floor = floor.at(column.pos_metrics.depth);
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
