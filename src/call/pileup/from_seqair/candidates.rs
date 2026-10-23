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
//!   `call::process::pileups`): every alt the first-stage forest scores at
//!   [`CASCADE_MIN_SCORE`] or more; without forests, [`snv_events`]:
//!   `CandidateFloor::SNV`, and a count a het could produce on the strand
//!   TAPS leaves informative. An SNV within `CLUSTER_GAP` of an indel chain
//!   joins the chain's events as context ([`snv_clusters`]); such a cluster
//!   ([`ClusterPass::Snvs`]) fills `snv_scores` and never writes an indel
//!   record.

use crate::call::pileup::indels::{IndelAllele, IndelAlleleCounts};
use crate::{
    call::pair_hmm::{Allele, Cluster, Event, MAX_EVENTS},
    metrics::{Alt, PileupMetrics},
};
use seqair_types::{Base, SmallVec, Strand};
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
    indel_floor: CandidateFloor,
) -> Vec<SnvCluster> {
    let pos = |i: usize| pileups.get(i).map(|p| p.pos);
    let near =
        |a: usize, b: usize| pos(a).zip(pos(b)).is_some_and(|(a, b)| b.abs_diff(a) <= CLUSTER_GAP);
    let mut joined = vec![false; snvs.len()];
    let mut clusters = Vec::new();
    for chain in cluster_chains(pileups, indel_columns) {
        let (Some(&first), Some(&last)) = (chain.first(), chain.last()) else { continue };
        let indels = cluster_candidates(chain.iter().filter_map(|&i| pileups.get(i)), indel_floor);
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
/// any anchor). An alt is admitted on Ben Wulf's gate: [`CandidateFloor::SNV`],
/// and on the strand TAPS leaves
/// informative for it, a count a heterozygote could produce. Without the gate
/// every methylated `CpG` and every one-read error would be a cluster.
pub(crate) fn snv_events(column: &PileupMetrics) -> SmallVec<Event, 2> {
    substitution_events(column, |alt| gate_admits(column, alt))
}

fn gate_admits(column: &PileupMetrics, alt: &Alt) -> bool {
    CandidateFloor::SNV.admits(alt.metrics.depth, column.pos_metrics.depth)
        && admits_het(column, alt)
}

/// The column's substitution events `admit` takes, most reads first, at most
/// two; none at an unknown reference base.
fn substitution_events(column: &PileupMetrics, admit: impl Fn(&Alt) -> bool) -> SmallVec<Event, 2> {
    if column.reference_base.known_index().is_none() {
        return SmallVec::new();
    }
    let mut admitted: SmallVec<&Alt, 4> = column
        .alts
        .iter()
        .filter(|alt| alt.base.known_index().is_some() && alt.base != column.reference_base)
        .filter(|alt| admit(alt))
        .collect();
    admitted.sort_by_key(|alt| (Reverse(alt.metrics.depth), alt.base as u8));
    admitted
        .into_iter()
        .take(2)
        .map(|alt| Event { anchor: column.pos, allele: Allele::Base(alt.base) })
        .collect()
}

/// The lowest first-stage forest score at which the cascade hands an alt to
/// the engine. The bundled forests put most true SNVs at .995-.999, so there
/// is no upper edge; .005 keeps 70-110 % of the gain of scoring every gated
/// alt for a quarter of its clusters on 26x slices, .01 loses some.
pub(crate) const CASCADE_MIN_SCORE: f64 = 0.005;

/// The events the cascade scores: with `forests`, every alt the first stage
/// scored at [`CASCADE_MIN_SCORE`] or more, whatever its reads, so every alt
/// the forests can pass carries the engine's genotype; without, the alts
/// [`snv_events`] admits.
pub(crate) fn cascade_events(column: &PileupMetrics, forests: bool) -> SmallVec<Event, 2> {
    if !forests {
        return snv_events(column);
    }
    substitution_events(column, |alt| alt.filters.ml.is_some_and(|ml| *ml >= CASCADE_MIN_SCORE))
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

/// Reads a candidate allele needs to become a pair-HMM event: `min_reads`, and
/// 5 % of the column's depth. Two reads of a stutter allele are noise at 73x
/// and the only evidence at 14x, so the floor grows with the column: for
/// indels, 5 % drops a quarter of the candidate events and 17 % of the
/// pair-HMM cells at 73x and changes no column of the chr12 scores at either
/// depth; 2 % drops nothing. Indel alleles count fragments, SNV alts reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CandidateFloor {
    pub min_reads: u32,
}

impl CandidateFloor {
    /// The SNV gate's floor without forests (`snv_events`).
    pub(crate) const SNV: Self = Self { min_reads: 2 };

    pub(crate) fn at(self, depth: u32) -> u32 {
        self.min_reads.max((depth * 5).div_ceil(100))
    }

    pub(crate) fn admits(self, reads: u32, depth: u32) -> bool {
        reads >= self.at(depth)
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
    let mut alleles: SmallVec<(u32, &IndelAlleleCounts), 8> = columns
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
