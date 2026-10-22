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
Rastair reads its pileups and writes VCF/BCF with seqair, our own crate (`../seqair`); rust-htslib
remains for BAM rewriting, `call-reads` and reading VCF/BCF.

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

- Prioritize code correctness and clarity. Speed and efficiency are secondary priorities unless otherwise specified.
- Do not write organizational comments or ones that summarize the code.
  - Comments should only be written in order to explain "why" the code is written in some way in the case there is a reason that is tricky / non-obvious.
  - In doc comments, do not write parameters and return type sections. Only add susprising constraints.
- Prefer implementing functionality in existing files unless it is a new logical component. Avoid creating many small files.
- Never use files with `mod.rs` paths - modules are always in `src/some_module.rs` instead of `src/some_module/mod.rs`.
- Avoid creative additions unless explicitly requested

## Error handling and logging

- Model the full error space—no shortcuts or simplified error handling. Use the type system to encode correctness constraints. Prefer compile-time guarantees over runtime checks where possible.
- Use `color_eyre` for error handling and reporting
- Avoid using functions that panic like `unwrap()`, instead use mechanisms like `?` to propagate errors.
- Don't use indexing operations, prefer methods like `get()` that return `Option` types.
- If you can't ensure correctness via the type system, use `ensure!` or `bail!` macros from `color_eyre` to handle unexpected states.
- Never silently discard errors with `let _ =` on fallible operations. Always handle errors appropriately:
  - Propagate errors with `?` when the calling function should handle them
  - Call `warn!(?error, "<what went wrong>")` or similar when you need to ignore errors but want visibility
  - Use explicit error handling with `match` or `if let Err(...)` when you need custom logic
- Use `tracing` for logging

## Testing

- Write comprehensive unit tests for the most critical and complex parts of the codebase when you either add them or encounter bugs in them
- Write integration tests for critical workflows and components, e.g. like the ones in `tests/call_cli.rs`
- Run the tests with `cargo test`.
- Use `cargo xtask insta` to run tests and update any snapshot tests. You need to verify the updated content is correct!

### VCF Tests

VCF tests are in `src/call/tests/vcf_tests/` with separate modules for different scenarios (cpgs.rs, denovo.rs, basic.rs).
Tests use the `pileups!` macro to create synthetic read data with format `[base1 base2 ...] Strand`, and `vcf_assert!` macro to check expected VCF output with format `(Ref Alt...) PASS/FAIL Field=value`.
Test utilities in `src/call/tests/utils.rs` provide the `pileups!` macro for creating test data, and helper functions like `set_pass`/`set_fail` for modifying alt calls with ML scores.
The `reprocess()` function recalculates methylation_strand_info, genotypes, alt calls, and methylation values after modifications.

## BAM rewriting

The BAM rewrite pipeline is in `src/bam.rs` with tag generation in `src/bam/base_modification.rs`.
There are two modes: `legacy` (XR/XG/XM tags, SEQ unchanged) and `standard` (MM/ML tags, SEQ rewritten T→C / A→G).

Critical invariant: **XM and MM/ML tags encode per-read methylation**, not position-level calls.
A CpG with beta=0.3 (`methylated: false` in `RastairCall`) still has individually methylated reads
that must show `Z` in XM and appear in MM/ML. The `methylated` field in `RastairCall::Cpg` only
controls whether the position is recognized as a CpG site, not whether individual reads are methylated.
Per-read methylation is determined solely by the observed base: T at OT C = methylated, A at OB G = methylated.

### MM/ML vs XM paired-read asymmetry

MM/ML tags only encode modifications at C bases in the stored SEQ. For paired reads overlapping a CpG,
only one mate has C at that position (the other has G on the complementary strand). XM tags annotate
both mates. This means tools reading MM/ML (like modkit) see roughly half the reads that XM-based
counting does. Methylation **fractions** agree, but exact counts differ ~2:1.

When comparing legacy (XM) and standard (MM/ML) output, always compare fractions, not absolute counts,
and require minimum coverage to avoid noise at low-coverage positions.

### External tool tests

Tests are in `tests/bam_external_tools.rs` behind `--features external-tool-tests`.
CI runs them in a dedicated `external-tools` job — separate from `test` so third-party
CLI drift does not redden the main test signal, but on a plain runner rather than in
Docker, sharing the `test` job's cargo cache. On Linux:

```bash
export PATH="$(.github/scripts/install-external-tools.sh):$PATH"
cargo test --features external-tool-tests
```

Each test self-skips when its tool is missing, so this is harmless elsewhere.

On **macOS** neither modkit nor the Bismark tarball has a build, so use `Dockerfile.ci`:

```bash
docker build -f Dockerfile.ci -t rastair-ext-tests .
docker run --rm -v "$(pwd):/rastair" rastair-ext-tests
```

That image is a fallback, not a mirror of CI, and nothing builds it automatically — which
is how the `tabix` bug below survived. It has no R, so `tests/mbias_report.rs` self-skips
there while CI renders the report; and it takes bismark/modkit from bioconda rather than
the versions `install-external-tools.sh` pins. When the two disagree, the `external-tools`
job is the source of truth. Header comment in the file has the details.

Two tool-version traps, both encoded in `install-external-tools.sh`:

- Use the **Perl** Bismark (`v0.25.x`), not the `bismark-rust-v3.x` rewrite — the rewrite aborts with
  "not yet implemented in this build: paired-end extraction ... PE arrives in Phase C", and the test BAM is paired.
- `modkit summary` dropped `--no-sampling`; `--sampling-frac 1` is the equivalent.

Cross-validation tests use `RASTAIR_TEST_MIN_COVERAGE` env var (default 5) to set the minimum read
coverage at a position before comparing fractions between tools. Lower values check more positions
but are noisier due to paired-read asymmetry.

# Interactivity guidelines

When you are asked to implement something, always ask for clarifications if needed.
If you are unsure about the requirements, ask for more details.
If you think there is a better way to implement something, suggest it and explain your reasoning, but don't implement it immediately without approval.

# Key data flow details

## Pileup construction and `Base::Unknown`

Columns are built in `src/call/pileup/from_seqair.rs` (`ColumnDraft::accumulate`, then `finish`).
`Pileup`, `SimpleRead` and `PileupMetrics::new` are test fixtures only (`#[cfg(test)]`), for `pileups!`.
The `reference_base` comes from the FASTA sequence via `Base::from(u8)`, which maps any non-ACGT byte (e.g. `N`) to `Base::Unknown`.
There is **no filtering** to skip columns with zero reads or Unknown reference bases before they are finished.

Important implications:

- `pileup.reference_base` can be `Base::Unknown` at N-positions in the reference — code must handle this gracefully (return default metrics), not treat it as an error.
- A pileup can have zero reads after filtering (all reads removed by quality/flag/overlap filters) — the zero-depth allele path is a real code path, not dead code.
- `Base::known_index()` maps A/C/G/T → `Some(0..3)` and Unknown → `None`. Use it to safely index into per-base arrays without needing an Unknown slot.

## Single-pass accumulator pattern

When computing grouped statistics from a collection of items (e.g. per-base metrics from reads), prefer a single-pass accumulator over collect-then-compute:

1. Create an accumulator struct with `Default` that holds incremental state (e.g. `RmsAccumulator`, running counts).
2. Feed items in one loop via an `add(&mut self, item)` method.
3. Finalize with `finish(self) -> Result<T>` that **takes `self` by value** to prevent accidental double-use.
4. When grouping by key (e.g. per-base), use `[Accumulator; N]` indexed by a method like `Base::known_index()` rather than named fields — this eliminates match arms for invalid variants and works naturally with const arrays like `Base::KNOWN`.
5. To extract a single group's accumulator, use a `take(&mut self, key) -> Option<Accumulator>` method via `mem::take` — `None` signals "not applicable" (e.g. Unknown base) rather than an error.

## Read orientation modes

The main pileup-based `call` path assigns OT/OB per record when seqair loads it (`RastairRecordFilter::compute` in `src/sequence/readers.rs`).
The default `VariantCallingParams.read_orientation=flags` path still uses `strand_from_flags()`.
The opt-in `--guess-read-orientation` mode does **not** require reference CpG annotation. Instead it scans each read over `aligned_pairs_full()` and only looks at read positions where the observed base mismatches the reference:

- at each mismatch, inspect both 2 bp windows that include that read base: current+next and previous+current
- count `TG` motifs and `CA` motifs in the read sequence as the BAM stores it (reference-oriented)
- `TG > CA` means OT, `CA > TG` means OB
- ties / no evidence: split deterministically from a hash of qname + start + flags so repeated runs stay reproducible

Implementation detail: the orientation is computed once per record as seqair loads it and kept in `RastairReadExtras`, because scoring mismatches per pileup column would rescan the full read for every covered base.

Current scope: this new evidence-based OT/OB assignment only affects the main pileup-based `call` path. `call-reads` and BAM rewriting still use their existing orientation logic.

For BAM-backed regression tests that compare strand-assignment modes, `tests/call_cli.rs` can write plain BED output with `call --cpgs-only --bed <path>` and compare per-CpG `(start, strand)` records via the BED columns `beta_est`, `unmod`, and `mod`. This is a convenient way to inspect differences before choosing hard thresholds.

## Driving seqair's pileup engine

`PileupEngine::new` takes a `PileupInput`, and the only way to make one is
`store.prepare_for_pileup()` — which returns `Prepared { input, stats }`, sorts
the store by position and links its mates. So a test that builds a store by hand
no longer has to remember either step, and can list its reads in any order.

That type exists because both preconditions used to fail silently and badly.
Measured on a three-read fixture: pushing them out of position order produced
columns for **one** of the three — the engine never reached the other two and no
column reported a gap. And on an unlinked store every alignment reports
`mate_idx() == None` *and* `in_mate_overlap() == false`, which is exactly what a
read with no mate looks like, so overlap dedup quietly does nothing and a test
asserting it passes while proving nothing.

Two related API notes:

- **`PileupColumn::mate_of(&view)`** gives the view's linked mate when it is also
  in this column — the whole of what `counted_mate` needs. It replaced
  `pair_indel`, a query that folded in the precedence "the view's own indel wins,
  the mate is never consulted"; that cannot serve a caller whose own filters may
  reject the view's indel. Keep pairwise rules on this side of the boundary.
- **`PileupEngine::reclaim_allocation`** (was `take_store`) returns an *empty*
  store keeping its slab capacity, for the next region. It is not a way to read
  the pileup's input back.

## Indel observations (`from_seqair.rs`)

Each rule below guards a *class* of mistake:

- **A deletion's REF comes from the deletion site.** `view.qpos()` is a read offset, not the
  column's genomic `pos`; seqair's `QPos` newtype keeps the two apart.
- **Tract runs are measured one base past the anchor**, where a left-aligned indel starts. The
  convention lives only in `ref_features::indel_tract_runs_at`.
- **Overlap dedup must not drop the fragment's only indel.** Dedup keeps one read per fragment by
  base agreement and template order, and can keep the mate that does *not* span the indel, so the
  dropped mate stands in: `build_indel_observation(view).or_else(mate)`. The fallback is keyed on
  the observation, not on the indel, because `build_indel_observation` rejects an indel near a read
  end or on a read with too many non-TAPS mismatches. seqair's `PileupColumn::pair_indel` ("own
  indel wins, the mate is never consulted") cannot express that, so rastair does not use it. On
  chr12 this rule is worth +2.3 insertion / +1.5 deletion F1.
- **One terminal-repeat predicate**, `has_terminal_repeat`, counted in whole units (4 for a
  homopolymer, 3 for a dinucleotide repeat), computed once per read when seqair loads it.
- **Fragment-level counts.** `soft_clip_count` and the noisy-reference count describe a fragment:
  whatever either mate shows counts once. A fragment is noisy reference only when it contributed
  no observation and a mate that is a reference alignment slips (terminal repeat or soft clip);
  `IndelCounts::clean_depth` subtracts a fragment counted on both sides twice otherwise.

Testing blind spots, still easy to reproduce:

- **`tests/data/test.bam` cannot produce an indel call.** It has 86 indel-carrying reads, but no
  two agree on an allele, so nothing passes `--min-indel-ao` at any threshold. An indel test has to
  build its own BAM; `indel_ref_alleles_come_from_the_deletion_site` is the CLI test that holds a
  multi-base REF.
- **A fixture segment starting at 0 hides every position bug**, because `pos - segment_start` and
  the clamped wrong answer coincide there. Use `segment_at(start, seq)` with a non-zero start; the
  CLI fixture puts its deletion at 24,137 for the same reason.
- **A `from_seqair` unit test must set `params.call_indels`**, or the indel branch is unreachable
  from it.

## ML feature layout (`src/metrics/ml/features/`)

Each model's feature vector is defined by a `#[repr(C)]` struct of `f32` / `[f32; N]`
fields built with the `define_features!` macro in `src/metrics/ml/features.rs`.
**The struct field order IS the feature vector order**, so there are no hand-counted
`buf[start..end]` index ranges anymore.

- The macro generates, per struct: `FEATURES` (from `size_of`), `names()`/`extend_names()`,
  and `as_row(&self) -> &[f32]` (zero-copy via `bytemuck::cast_slice`; the struct is `Pod`
  because all fields are `f32` and there is no padding).
- Field kinds in the macro: `scalar name;` (one feature, named after the ident),
  `array name: N = ["..", ..];` (N features with explicit per-slot names), and
  `flatten name: Type;` (embeds a nested feature struct and delegates its names).
- `CommonFeatures` (in `shared.rs`) is the shared extractor; its layout is split into
  `CommonSectionA` (33) + `CommonSectionB` (18) because the alt-based models interleave
  model-specific scalars (e.g. `alt_score`) _between_ the two halves. Build them via
  `CommonSectionA::from_common(&common)`.
- Model structs: `CpgFeatures` (55), `DenovoCpgFeatures` (56), `OthersFeatures` (54),
  `InsertionFeatures` (34), `DeletionFeatures` (38). Each has an `extract()` returning the
  struct; `FeatureCalculator::calculate_*` wraps `as_row()` into an `Array2`.

**Feature order is frozen by every trained model.** Reordering a field silently corrupts
predictions. Two tests guard this in `features.rs`: `feature_counts_are_stable` (pins the
counts) and `feature_name_layout` (an insta snapshot of every `name→index` mapping —
this replaces the old "never change the order" comments; update it via `cargo xtask insta`
only after verifying a layout change is intentional).

Feature names flow to training output via `FeatureCalculator::feature_names() -> FeatureNames`.
`train.rs` uses them for the `--feature-analytics` importance CSVs (`index\tfeature\timportance`)
and the `--export-features` TSV headers, so both exports agree by construction.

## Training a model (`rastair ml train`)

```bash
rastair ml train -r hg38.fa.gz taps.bam HG001_benchmark.vcf.gz \
  -R HG001_benchmark.bed -l "chr1 chr6 chr11" --seed <n> -o models/mymodel.rff.mpk.lz4 -@ 8
```

- Positionals are **BAM first, then the truth VCF**. The CLI defaults are
  the recipe of the bundled model (`models/rastair_with_indels.rff.mpk.lz4`, `include_bytes!`'d
  by `src/call/ml.rs`), so a bare command line retrains it.
- The truth VCF needs a **`.csi`** (a `.tbi` alone is refused), and `-R` must be a plain-text
  BED, not bgzipped. Without `-R`, candidates outside the truth set's confident regions are
  labelled negative.
- Train on GIAB HG001 and evaluate against Platinum Genomes, holding out the evaluation
  chromosome. **Give it enough genome**: all five forests must train or no model is written.
- Collection is bounded: rows are `f32` in one flat buffer, each model keeps per class the
  examples with the smallest uniform keys (a keyed reservoir, capped at the draw plus a
  holdout), and segments are folded with `reduce`, not collected. `KeySource` seeds one stream
  per segment *and* model from `--seed`, so the sample does not depend on scheduling, and one
  model's pool does not move another's draw.
- The draw takes at most four fifths of a class, so small models keep a holdout; the holdout
  keeps the **population's** class ratio, not the reservoir's near-balanced one, or Platt
  scaling would shift what an ML threshold means.
- SNV and indel models are drawn separately (`--n-*` vs `--indel-n-*`): SNV candidates are
  almost all negative, indel candidates arrive filtered and mostly positive. The log line
  `Collected training examples` reports pool, kept, requested and drawn per model.
- `--export-features` writes a 1-based `pos`, like a VCF, and no allele column.

## VCF header version and cardinality

Output is **`##fileformat=VCFv4.5`**, and seqair writes that unconditionally —
`VcfHeader::FILE_FORMAT`, with no setter — so rastair does not ask for it. That
is the version defining the fields we emit: `M5mC`, `DPM5mC` and `ADM5mC` are
VCF 4.5 reserved FORMAT keys, aliases for the ChEBI-numbered `M27551C` family,
and 4.3 defines none of them.

Declaring 4.5 is not cosmetic. From VCF 4.4 the **first allele's `GT` phase bit
is read** rather than ignored, so an encoder that leaves it unset makes htslib
render a phased `0|1` as `/0|1`. seqair sets it now (spec rule
`vcf_record.gt_first_phase`); if phased output ever comes back looking like
that, this is why.

**Keep those three at `Number=.` anyway.** VCF 4.5 pairs them with `Number=M`
("one value for each possible base modification for the corresponding ChEBI
ID"), and seqair still has a `Number::BaseModification` variant that emits it —
but noodles rejects `M`, and **declaring 4.5 does not help.** Measured on one of
our own files with only the header rewritten:

| `##fileformat` | `Number` on M5mC | aardvark |
| --- | --- | --- |
| VCFv4.3 | `M` | rejects the whole file |
| VCFv4.5 | `M` | **still rejects** |
| VCFv4.5 | `.` | reads it |

The rejection is `invalid FORMAT: ID=M5mC: invalid number`, and it kills the
*file*, not the line — so every noodles-based tool, PacBio's `aardvark`
included, sees nothing. `.` states the same cardinality in a way every reader
accepts. Revisit when noodles implements the 4.5 cardinalities; the two tests in
`src/vcf/schema.rs` pin both halves of this.

## VCF FILTER is a set (`RastairFilter` / `Filters`)

`Filters` (`src/metrics/pileup_metrics.rs`) is an `EnumSet<RastairFilter>` — a
`u16` bitset, pinned by `#[enumset(repr = "u16")]` on the enum — plus the
`other_pos_in_denovo_passes` override, which is *not* a FILTER code and so does
not live in the set. `add`/`merge` are `insert`/`|=`; there is no dedup to do by
hand and no `Deref` to a list any more.

**The enum's declaration order is load-bearing twice.** `RastairFilter as usize`
indexes `Schema::filter`'s `[FilterId; COUNT]` table (so the discriminants must
stay `0..COUNT` — this is why `enumset` fits and `enumflags2`, which wants
power-of-two discriminants, does not), and a set iterates in discriminant order,
which is the order the FILTER column prints. Reordering variants is an output
change.

**Nothing snapshots a non-PASS FILTER column.** Rejected records are only
emitted under `--all` (`emit_rejected_record`, gated by `RecordFilters`), and
every committed VCF snapshot is 100 % `PASS`. To see FILTER output at all:

```bash
cargo build && ./target/debug/rastair call --fasta-file=tests/data/test.fasta.gz \
  tests/data/test.bam --all | grep -v '^#' | awk -F'\t' '$7!="PASS"{print $7}' \
  | sort | uniq -c | sort -rn
```

That blind spot hid a real defect until 2026-09-08: FILTER used to be built by
appending three lists, so a code could land twice — *every* non-PASS record on
`tests/data/test.bam` carried `low_ml_score;low_ml_score`. Use the command above
when touching filter emission.

Still open, found while fixing that: `emit_rejected_record` adds `low_ml_score`
when `alt.filters.ml < ml_threshold`, and `None < Some(_)` in Rust — so a record
whose ML was *skipped* (`pre_ml`) is also labelled `low_ml_score`. Fixing it
changes `--all` output beyond a reordering, so it was left alone.

## Release version bump checklist

When bumping Rastair's release version, update all user-facing version strings together:

- Root crate version in `Cargo.toml` (`[package].version`)
- Root package entry in `Cargo.lock` (`name = "rastair"`)
- CLI docs version in `docs/src/cli.md`
- README example tag references in `README.md` (e.g. `version-X.Y.Z`)
- Snapshot VCF header lines in `tests/snapshots/` containing `##rastairVersion=...`

The release workflow refuses to build a tag whose name does not equal
`v` + `[package].version` from `Cargo.toml`, so a forgotten bump fails fast.

## CI (GitHub Actions)

CI lives in `.github/`. See `.github/README.md` for the secrets/variables a release needs.

## CLI docs generation

The command-line reference at `docs/src/cli.md` is generated from clap doc comments.
Use the hidden command:

- `cargo run -- internal cli-docs docs/src/cli.md`

Toolchain note: `rust-toolchain.toml` pins the compiler, so `cargo run` picks it up automatically.

## QC report M-bias orientation

In `scripts/QC_report.Rmd`, OT/OB assignment for the M-bias table must use the same pair-orientation logic as Rust:

- OT if `bitwAnd(flag, 96) == 96` (F1R2) or `bitwAnd(flag, 144) == 144` (R2F1)
- OB otherwise

Using `80/160` (first+reverse / second+mate_reverse) swaps OT and OB labels and flips the wrong mate.

To plot/read cutoffs in read 5'->3' coordinates, flip positions for reverse-aligned mates only:

- OT + `Second`
- OB + `First`

## QC report (`rastair mbias`) architecture and testing

The `mbias` subcommand (`src/mbias.rs`) has **no analysis logic of its own** — it only shells out to `scripts/mbias.R`, which renders `scripts/QC_report.Rmd`. All M-bias cutoff math, plotting, and the per-contig `{chrom}_cutoffs.txt` outputs live in the R code. Fixes to cutoff/plot behaviour go in the `.Rmd`, not Rust.

- **Per-contig × per-group cutoffs.** The `plot_mbias` chunk computes cutoffs per `(chr, read_pair, strand)` group (up to 4 groups/contig). A group needs ≥`MIN_MBIAS_OBS` (3) covered read positions. Sparse groups (tiny alt/decoy contigs) used to `stop()` and abort the **entire** report. Now the behaviour depends on whether the run was scoped: for a **genome-wide run** (no `--region`) the chunk skips the whole sparse contig (no plot, no cutoffs file) and `warning()`s; when a **`--region`/chromosome was explicitly requested** a sparse contig is still a hard `stop()` (the user asked for exactly that data). `calculate_cutoff()` is also total (returns `left=0,right=0` instead of `stop()`).
- **`--plot-fp` is opt-in.** In `mbias.R`, `plot_fp` must be `isTRUE(args$plot_fp)`, not `!is.na(...)` — argparser `flag=TRUE` args default to `FALSE` (not `NA`), so `!is.na()` is always TRUE and forces the false-positives plot on, which aborts any `--bed`-only run lacking a vcf/bam.
- **Render path needs vcf/bam only for some chunks.** A `--bed`-only render skips V-bias/GC/FP chunks (gated by `params$plot_vbias`/`plot_gc`/`plot_fp`); `mbias.rs` auto-adds `--no-vbias`/`--no-gc` when no `--reference` is given.

### Testing the report

The render is exercised by `tests/mbias_report.rs`, gated behind the `external-tool-tests` feature and **self-skipping** when `Rscript` (+ `rmarkdown`/`argparser`/`data.table`/`ggplot2`), `tabix`, or `bgzip` are missing. It writes a synthetic per-read BED (header mirrors `PerRead::HEADER` in `src/bed/per_read/format.rs`) with one healthy and one sparse contig, then asserts the healthy contig gets a `*_cutoffs.txt` and the sparse one does not.

- **macOS caveat:** the `.Rmd` loads the no-region input via `zcat <bgz>`, which on macOS only handles `.Z`, not `.gz`. The test therefore only renders cleanly on Linux/Docker (where `external-tool-tests` are meant to run). To run it locally on macOS, put a `zcat` shim that execs `gzip -dc` early on `PATH`.
- To verify an `.Rmd` change quickly without the `mbias.R`/argparser wrapper, render directly: `Rscript -e "rmarkdown::render('scripts/QC_report.Rmd', params=list(input_bgz=..., output_dir=..., region=NA, plot_vbias=FALSE, plot_gc=FALSE))"` (pass `region=NA`, not NULL).

# Keep this updated

**Important:** Whenever you learned something new about how to develop features, find code, or how to debug issues, you **must** add it to this document.
This is the single source of truth for how to work on this codebase, and it must be kept up-to-date with any new insights or changes.
If you find yourself asking "How do I do X?" and you figure it out, add a section here so that the next person doesn't have to ask the same question.
