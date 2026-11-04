# What is Rastair?

Rastair is a Rust CLI that detects genetic variants and methylated positions at the same time from short-read TAPS sequencing data.

TAPS converts methylated C to T, while unmethylated C is converted to U and then read as C.
So methylation shows as T reads at a reference C on the OT strand, and as A reads at a reference G on the OB strand.
De-novo CpG positions arise when X→C or X→G SNPs occur.
C→T and G→A SNPs are confounded with methylation, which is what makes variant calling hard.
Rastair combines machine learning with a TAPS-aware pair-HMM to decide which variants are real.

Pileups are read and VCF/BCF is written with seqair.
rust-htslib does the BAM rewriting, `per-read` and the reading of VCF/BCF.

## Structure

`call` runs per segment (`process_segment` in `src/call.rs`):

1. `process::get_pileups_held` builds each column from seqair (`src/call/pileup/from_seqair.rs`), drops the columns with no alt, indel or CpG (`ColumnEvidence`, `Columns::keeps` in `src/call/process/pileups.rs`), finishes the rest into `PileupMetrics`, and calls the indels.
2. The forests score every alt; the SNV engine scores the ones they kept, and the engine forests score those again (`src/call/process/pileups.rs`).
3. `process_collected_pileups`: `call_columns` (alt calls with the engine's genotype → CpGs of called alts → `0/0` genotype, `src/call/variant_calling/genotype.rs` → methylation, `src/metrics/methylation.rs`), then phasing, then the overlap is trimmed.

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
| readers, strand assignment                     | `src/sequence/readers.rs`, `docs/src/calling/strand-guesser.md` |
| BAM rewriting, XM vs MM/ML                     | `src/bam.rs`                                                    |
| external tool tests                            | `tests/bam_external_tools.rs`, `Dockerfile.ci`                  |
| ML feature layout                              | `src/metrics/ml/features.rs`                                    |
| training a model                               | `src/train.rs`                                                  |
| VCF version and field cardinality              | `src/vcf/schema.rs`                                             |
| FILTER codes                                   | `src/vcf/filter.rs`                                             |
| QC report                                      | `src/mbias.rs`                                                  |
| VCF tests and `pileups!`                       | `src/call/tests/vcf_tests.rs`, `src/call/tests/pileup.rs` (`Pileup::metrics`) |
| CI and release secrets                         | `.github/workflows/README.md`                                   |

User documentation is in `docs/src`; `docs/src/cli.md` is generated from the clap doc comments
with `cargo run -- internal cli-docs docs/src/cli.md`, and `docs/src/formats/vcf-fields.md` from
`src/vcf/schema.rs` with `cargo run -- internal vcf-docs docs/src/formats/vcf-fields.md`
(`cargo xtask docs` runs both, then mdbook).

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
- Run `cargo nextest run`, `cargo test --doc` (nextest skips doctests) and `cargo clippy --all-targets`.
- Prefer property-based tests with `hegel` (`hegeltest`) where applicable: `#[hegel::test] fn f(tc: TestCase)` with `tc.draw(gs::..)`, `#[hegel::composite]` for custom generators, plain `assert!`/`assert_eq!`, `tc.assume(..)` to reject a case. Case counts live in `hegel.toml`; an expensive test pins its own with `pinned::cases(N)` (`src/pinned.rs` says why and how). Prefer constructing inputs in the shape you need over rejecting them — a rejection loop makes hegel's smallest generated input large and its shrinking useless. Draw a length as an integer when it matters, since hegel keeps a vector's own size short.
- `cargo xtask insta` runs the tests and updates snapshots; verify the updated content. The `call` VCF snapshots carry the config in their header, so a new CLI option moves them all.

# Interactivity guidelines

When asked to implement something, ask for clarification if needed. If you think there is a better way, suggest it and explain why, but don't implement it without approval.

# Release version bump checklist

Update all user-facing version strings together; the release workflow refuses a tag that is not `v` + `[package].version`, so a forgotten bump fails fast.

- `[package].version` in `Cargo.toml` and the `rastair` entry in `Cargo.lock`
- `version` in `CITATION`
- the version in `docs/src/cli.md` (regenerate it)
- tag references in `README.md` (e.g. `version-X.Y.Z`)
- `##rastairVersion=...` in the VCF snapshots under `tests/snapshots/` (`cargo xtask insta`)
- the `Unreleased` heading of `docs/src/CHANGELOG.md` becomes `Version X.Y.Z (date)`

# Keep this updated

When you learn something about developing, finding code or debugging, write it down: in the doc comment of the module it concerns, or here if it cuts across the codebase.
