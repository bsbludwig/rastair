//! Columns from seqair's pileup engine, accumulated into
//! [`PileupMetrics`](crate::metrics::PileupMetrics), and the pair-HMM
//! clusters scored over them.
//!
//! - [`column`]: one column's reads, through [`ColumnDraft`], into metrics;
//!   driving the engine, overlap dedup, mates, indel observations.
//! - [`candidates`]: which indel alleles and SNV alts become pair-HMM events,
//!   and which columns chain into one cluster.
//! - [`cluster`]: a cluster's reads and haplotype window
//!   ([`ClusterContext`], `WindowMask`), scored and handed to its columns.
//! - [`gpu`]: the pair-HMM on the GPU (`--gpu`).
//!
//! Testing: `tests/data/test.bam` cannot produce an indel call (no two of its
//! indel reads agree on an allele), so indel tests build their own reads or a
//! `Scenario`; a fixture segment must not start at 0, where position bugs and
//! the right answer coincide (use `segment_at`); and a unit test must set
//! `params.call_indels` to reach the indel branch at all. The cluster pieces
//! have unit tests here and in `pair_hmm`; the CLI scenario
//! `the_pair_hmm_calls_two_events_of_one_cluster` fails with one event per
//! haplotype.

mod candidates;
mod cluster;
mod column;
mod gpu;
#[cfg(test)]
mod tests;

#[cfg(test)]
pub(crate) use candidates::{CASCADE_MIN_SCORE, snv_events};
pub(crate) use candidates::{
    CandidateFloor, ClusterPass, cascade_events, cluster_candidates, cluster_chains, snv_clusters,
};
pub(crate) use cluster::{ClusterContext, flush_clusters, score_cluster};
pub(crate) use column::{ColumnDraft, ColumnScratch, DenovoNeighbour, plain_reference_column};
