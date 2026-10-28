//! Tests for VCF output from rastair call, one module per scenario.
//!
//! `pileups!` builds synthetic columns (`[base1 base2 ...] Strand`) and
//! `vcf_assert!` checks the output (`(Ref Alt...) PASS/FAIL Field=value`);
//! both and helpers like `set_pass`/`set_fail` are in `super::utils`.
//! `reprocess()` runs the calls again after a fixture is modified: alt calls,
//! the CpGs of called alts, the `0/0` genotype and methylation.
//! `Pileup::metrics` builds each column with the production `ColumnDraft`
//! (its doc says what a fixture leaves out). `test_call` runs `call`'s cascade
//! on them: the forests, the SNV engine (`score_snv_lines`, on every column's
//! two most-read alts), and the engine forests on the alts the first stage
//! scored at `CASCADE_MIN_SCORE`.

mod basic;
mod ben_edge_cases;
mod cpgs;
mod denovo;
mod edge_cases;
mod genotyping;
mod indels;
mod phasing;
mod real_hg96_examples;
