//! Tests for VCF output from rastair call, one module per scenario.
//!
//! `pileups!` builds synthetic columns (`[base1 base2 ...] Strand`) and
//! `vcf_assert!` checks the output (`(Ref Alt...) PASS/FAIL Field=value`);
//! both and helpers like `set_pass`/`set_fail` are in `super::utils`.
//! `reprocess()` recomputes strand info, genotypes, alt calls and methylation
//! after a fixture is modified. The columns go through the test-only
//! `Pileup` → `PileupMetrics::new` path, not seqair.

mod basic;
mod ben_edge_cases;
mod cpgs;
mod denovo;
mod edge_cases;
mod genotyping;
mod indels;
mod phasing;
mod real_hg96_examples;
