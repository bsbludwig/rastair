//! Tests for VCF output from rastair call, one module per scenario.
//!
//! `pileups!` builds synthetic columns (`[base1 base2 ...] Strand`) and
//! `vcf_assert!` checks the output (`(Ref Alt...) PASS/FAIL Field=value`);
//! both and helpers like `set_pass`/`set_fail` are in `super::utils`.
//! `reprocess()` recomputes strand info, genotypes, alt calls and methylation
//! after a fixture is modified. `Pileup::metrics` turns each column's reads
//! into seqair records and builds its metrics with the production
//! `ColumnDraft`: one fragment per read, no base-quality or read-end filter.

mod basic;
mod ben_edge_cases;
mod cpgs;
mod denovo;
mod edge_cases;
mod genotyping;
mod indels;
mod phasing;
mod real_hg96_examples;
