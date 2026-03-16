# What is Rastair?

Rastair is a CLI application written in Rust that allows
the simultaneous detection of genetic variants and methylated positions
from short-read sequencing data created using the TAPS method.

## Methylation and variant calling

TAPS converts methylated C to T, while unmethylated C is converted to U and then read as C.
Thus, methylation is evidenced by having a C reference position show T reads on the OT strand, G refs show A reads on OB strand.
In addition, de-novo CpG postions are possible when X->C or X->G SNPs occur.
Variant calling is complicated by the fact that C->T and G->A SNPs are confounded with methylation.

Rastair uses a combination of thresholding and machine learning to determine true variants.
Rastair's main feature is `call` which processes pileup data
in multiple steps and produces VCF records with variant and methylation calls.
Rastair uses htslib via rust-htslib for reading/writing BAM and VCF/BCF files.

## Structure

Rastair is structured as a CLI application using `clap` for argument parsing.
The main functionality is implemented in the `rastair` crate,
with submodules for different components like pileup processing, variant calling, and methylation analysis.
The `xtask` crate is used for auxiliary tasks like testing and benchmarking.

The core processing pipeline is implemented in `src/call.rs` in the `process_region` function, which processes pileups through several stages:
calculate pileup metrics → set de-novo adjacency flags → add ML metrics → apply threshold filters → propagate de-novo pass flags → set alt calls → estimate genotype → call methylation.
Methylation calling logic is in `src/metrics/methylation.rs` with separate functions for reference C/G positions (`ref_c`, `ref_g`) and de-novo CpG creation (`ref_t_to_c`, `ref_a_to_g`, etc.).
Genotype estimation happens before methylation calling in `src/call/variant_calling/genotype.rs`.
The results are stored in `PileupMetrics.pos_metrics.extended.genotype` and `.methylated`.

# Rust coding guidelines

Rust code is to be written in expert-level Rust. Use the most modern features and idioms.
Specific adn well-named types are the main way to ensure correctness and introduce abstraction.

## General style

* Prioritize code correctness and clarity. Speed and efficiency are secondary priorities unless otherwise specified.
* Do not write organizational comments or ones that summarize the code.
  * Comments should only be written in order to explain "why" the code is written in some way in the case there is a reason that is tricky / non-obvious.
  * In doc comments, do not write parameters and return type sections. Only add susprising constraints.
* Prefer implementing functionality in existing files unless it is a new logical component. Avoid creating many small files.
* Never use files with `mod.rs` paths - modules are always in `src/some_module.rs` instead of `src/some_module/mod.rs`.
* Avoid creative additions unless explicitly requested

## Error handling and logging

* Model the full error space—no shortcuts or simplified error handling. Use the type system to encode correctness constraints. Prefer compile-time guarantees over runtime checks where possible.
* Use `color_eyre` for error handling and reporting 
* Avoid using functions that panic like `unwrap()`, instead use mechanisms like `?` to propagate errors.
* Don't use indexing operations, prefer methods like `get()` that return `Option` types.
* If you can't ensure correctness via the type system, use `ensure!` or `bail!` macros from `color_eyre` to handle unexpected states.
* Never silently discard errors with `let _ =` on fallible operations. Always handle errors appropriately:
  - Propagate errors with `?` when the calling function should handle them
  - Call `warn!(?error, "<what went wrong>")` or similar when you need to ignore errors but want visibility
  - Use explicit error handling with `match` or `if let Err(...)` when you need custom logic
* Use `tracing` for logging

## Testing

* Write comprehensive unit tests for the most critical and complex parts of the codebase when you either add them or encounter bugs in them
* Write integration tests for critical workflows and components, e.g. like the ones in `tests/call_cli.rs`
* Run the tests with `cargo test`.
* Use `cargo xtask insta` to run tests and update any snapshot tests. You need to verify the updated content is correct!

### VCF Tests

VCF tests are in `src/call/tests/vcf_tests/` with separate modules for different scenarios (cpgs.rs, denovo.rs, basic.rs).
Tests use the `pileups!` macro to create synthetic read data with format `[base1 base2 ...] Strand`, and `vcf_assert!` macro to check expected VCF output with format `(Ref Alt...) PASS/FAIL Field=value`.
Test utilities in `src/call/tests/utils.rs` provide the `pileups!` macro for creating test data, and helper functions like `set_pass`/`set_fail` for modifying alt calls with ML scores.
The `reprocess()` function recalculates methylation_strand_info, genotypes, alt calls, and methylation values after modifications.

# Interactivity guidelines

When you are asked to implement something, always ask for clarifications if needed.
If you are unsure about the requirements, ask for more details.
If you think there is a better way to implement something, suggest it and explain your reasoning, but don't implement it immediately without approval.

# Key data flow details

## Pileup construction and `Base::Unknown`

`Pileup` objects are constructed in `src/call/pileup/from_hts.rs` via `Pileup::from_hts()`.
The `reference_base` comes from the FASTA sequence via `Base::from(u8)`, which maps any non-ACGT byte (e.g. `N`) to `Base::Unknown`.
There is **no filtering** to skip pileups with zero reads or Unknown reference bases before `PileupMetrics::new()` is called.

Important implications:
* `pileup.reference_base` can be `Base::Unknown` at N-positions in the reference — code must handle this gracefully (return default metrics), not treat it as an error.
* A pileup can have zero reads after filtering (all reads removed by quality/flag/overlap filters) — the zero-depth allele path is a real code path, not dead code.
* `Base::known_index()` maps A/C/G/T → `Some(0..3)` and Unknown → `None`. Use it to safely index into per-base arrays without needing an Unknown slot.

## Single-pass accumulator pattern

When computing grouped statistics from a collection of items (e.g. per-base metrics from reads), prefer a single-pass accumulator over collect-then-compute:
1. Create an accumulator struct with `Default` that holds incremental state (e.g. `RmsAccumulator`, running counts).
2. Feed items in one loop via an `add(&mut self, item)` method.
3. Finalize with `finish(self) -> Result<T>` that **takes `self` by value** to prevent accidental double-use.
4. When grouping by key (e.g. per-base), use `[Accumulator; N]` indexed by a method like `Base::known_index()` rather than named fields — this eliminates match arms for invalid variants and works naturally with const arrays like `Base::KNOWN`.
5. To extract a single group's accumulator, use a `take(&mut self, key) -> Option<Accumulator>` method via `mem::take` — `None` signals "not applicable" (e.g. Unknown base) rather than an error.

# Keep this updated

**Important:** Whenever you learned something new about how to develop features, find code, or how to debug issues, you **must** add it to this document.
This is the single source of truth for how to work on this codebase, and it must be kept up-to-date with any new insights or changes.
If you find yourself asking "How do I do X?" and you figure it out, add a section here so that the next person doesn't have to ask the same question.
