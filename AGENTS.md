# What is Rastair?

Rastair is a Rust CLI that detects genetic variants and methylated positions at the same time from short-read TAPS sequencing data.

TAPS converts methylated C to T, while unmethylated C is converted to U and then read as C.
So methylation shows as T reads at a reference C on the OT strand, and as A reads at a reference G on the OB strand.
De-novo CpG positions arise when X→C or X→G SNPs occur.
C→T and G→A SNPs are confounded with methylation, which is what makes variant calling hard.
Rastair combines thresholds and machine learning to decide which variants are real.

Pileups are read and VCF/BCF is written with seqair.
rust-htslib is still used for BAM rewriting, `call-reads` and reading VCF/BCF.

## Structure

`call` runs per segment in `src/call.rs`:

1. `process::get_pileups` builds each column from seqair (`src/call/pileup/from_seqair.rs`), drops the columns the record filters would reject, and finishes the rest into `PileupMetrics`.
2. `process_pre_built_metrics` sets de-novo adjacency and strand info and applies the pre-filter.
3. `process_collected_pileups`: ML scores → threshold filters → de-novo pass propagation → alt calls → genotype (`src/call/variant_calling/genotype.rs`) → methylation (`src/metrics/methylation.rs`).

Area-specific knowledge lives in module doc comments; read the one for the code you touch:

| area                                           | read                                                            |
| ---------------------------------------------- | --------------------------------------------------------------- |
| columns, the seqair engine, indel observations | `src/call/pileup/from_seqair/column.rs` (index: `from_seqair.rs`) |
| candidates, clusters, window mask, GPU lane    | `src/call/pileup/from_seqair/{candidates,cluster,gpu}.rs`       |
| pair-HMM: haplotypes, verdict, dosage, CpGs    | `src/call/pair_hmm.rs`, `src/call/pair_hmm/scoring.rs`          |
| indel records                                  | `src/call/variant_calling/indel_calling.rs`                     |
| SNV engine cascade                             | `src/call/process/pileups.rs`                                   |
| genotype rules, GQ                             | `src/call/variant_calling/genotype.rs`                          |
| phasing                                        | `src/call/phasing.rs`, its `problem` and `solver` modules       |
| measuring accuracy, comparing binaries         | `docs/measuring-accuracy.md`                                    |
| readers, strand assignment                     | `src/sequence/readers.rs`                                       |
| BAM rewriting, XM vs MM/ML                     | `src/bam.rs`                                                    |
| external tool tests                            | `tests/bam_external_tools.rs`, `Dockerfile.ci`                  |
| ML feature layout                              | `src/metrics/ml/features.rs`                                    |
| training a model                               | `src/train.rs`                                                  |
| VCF version and field cardinality              | `src/vcf/schema.rs`                                             |
| FILTER codes                                   | `src/vcf/filter.rs`                                             |
| QC report                                      | `src/mbias.rs`                                                  |
| VCF tests and `pileups!`                       | `src/call/tests/vcf_tests.rs`, `src/call/pileup.rs` (`Pileup::metrics`) |
| CI and release secrets                         | `.github/workflows/README.md`                                   |

User documentation is in `docs/src`; `docs/src/cli.md` is generated from the clap doc comments
with `cargo run -- internal cli-docs docs/src/cli.md`.

# Rust coding guidelines

Write expert-level Rust with modern features and idioms. Specific, well-named types are the main
way to ensure correctness and to introduce abstraction.

- Prioritize correctness and clarity. Speed is secondary unless stated otherwise.
- Comments explain _why_, where the reason is tricky or non-obvious; no organizational or summarizing comments.
  They describe the code as it is, not earlier implementations or fixed bugs.
  In doc comments, no parameter or return sections; only surprising constraints.
- Prefer existing files unless the code is a new logical component. Never use `mod.rs`: a module is `src/some_module.rs`.
- Avoid creative additions unless explicitly requested.

## Error handling and logging

- Model the full error space; encode correctness in types and prefer compile-time guarantees.
- Use `color_eyre` for errors and `tracing` for logging.
- No `unwrap()` or indexing: propagate with `?`, use `get()`; `ensure!`/`bail!` where types cannot help.
- Never discard errors with `let _ =`: propagate, log with `warn!(?error, "...")` when ignoring on purpose, or handle with `match` / `if let Err(...)`.

## Testing

- Write unit tests for the critical and complex parts when you add them or find bugs in them, and integration tests for critical workflows (like `tests/call_cli.rs`).
- Run `cargo test` and `cargo clippy --all-targets`.
- `cargo xtask insta` runs the tests and updates snapshots; verify the updated content. The `call` VCF snapshots carry the config in their header, so a new CLI option moves them all.

# Interactivity guidelines

When asked to implement something, ask for clarification if needed. If you think there is a better way, suggest it and explain why, but don't implement it without approval.

# Release version bump checklist

Update all user-facing version strings together; the release workflow refuses a tag that is not `v` + `[package].version`, so a forgotten bump fails fast.

- `[package].version` in `Cargo.toml` and the `rastair` entry in `Cargo.lock`
- the version in `docs/src/cli.md`
- tag references in `README.md` (e.g. `version-X.Y.Z`)
- `##rastairVersion=...` in the VCF snapshots under `tests/snapshots/`

# Keep this updated

When you learn something about developing, finding code or debugging, write it down: in the doc comment of the module it concerns, or here if it cuts across the codebase.
