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
With `--phase`, `phasing::phase_segment` runs *before* the overlap trim, given the segment's core
predicate so it only writes calls for the sites the segment emits.

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
- Run the tests with `cargo test` **and** `cargo test --features experimental-seqair`. Phasing is seqair-only, so its tests are `#[cfg(feature = "experimental-seqair")]` and a plain `cargo test` never runs them; the htslib build instead tests that `--phase` is refused.
- Run `cargo clippy --all-targets` under both feature sets too.
- Use `cargo xtask insta` to run tests and update any snapshot tests. You need to verify the updated content is correct!

### VCF Tests

VCF tests are in `src/call/tests/vcf_tests/` with separate modules for different scenarios (cpgs.rs, denovo.rs, basic.rs, phasing.rs).
Tests use the `pileups!` macro to create synthetic read data with format `[base1 base2 ...] Strand`, and `vcf_assert!` macro to check expected VCF output with format `(Ref Alt...) PASS/FAIL Field=value`.
Test utilities in `src/call/tests/utils.rs` provide the `pileups!` macro for creating test data, and helper functions like `set_pass`/`set_fail` for modifying alt calls with ML scores.
The `reprocess()` function recalculates methylation_strand_info, genotypes, alt calls, and methylation values after modifications.

An htslib-shaped `Pileup` has no template identity, so `pileups!` fixtures cannot carry phase
observations; phasing VCF tests set the phase call by hand. The solver is tested end to end with
scenarios (see "Testing phasing").

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
`bcftools_reads_phased_genotypes` additionally needs `experimental-seqair`.

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

`Pileup` objects are constructed in `src/call/pileup/from_hts.rs` via `Pileup::from_hts()`.
The `reference_base` comes from the FASTA sequence via `Base::from(u8)`, which maps any non-ACGT byte (e.g. `N`) to `Base::Unknown`.
There is **no filtering** to skip pileups with zero reads or Unknown reference bases before `PileupMetrics::new()` is called.

Important implications:

- `pileup.reference_base` can be `Base::Unknown` at N-positions in the reference — code must handle this gracefully (return default metrics), not treat it as an error.
- A pileup can have zero reads after filtering (all reads removed by quality/flag/overlap filters) — the zero-depth allele path is a real code path, not dead code.
- `Base::known_index()` maps A/C/G/T → `Some(0..3)` and Unknown → `None`. Use it to safely index into per-base arrays without needing an Unknown slot.

`PileupMetrics` is pinned in size (`pileup_metrics_stays_small`): a region holds one per
covered base and the pipeline walks that vec many times. Per-feature data a run may not use
lives behind one `Option<Box<_>>` (`phase: Option<Box<ColumnPhase>>`,
`indel_data: Option<Box<IndelData>>`), so an unphased run pays one pointer.

## Single-pass accumulator pattern

When computing grouped statistics from a collection of items (e.g. per-base metrics from reads), prefer a single-pass accumulator over collect-then-compute:

1. Create an accumulator struct with `Default` that holds incremental state (e.g. `RmsAccumulator`, running counts).
2. Feed items in one loop via an `add(&mut self, item)` method.
3. Finalize with `finish(self) -> Result<T>` that **takes `self` by value** to prevent accidental double-use.
4. When grouping by key (e.g. per-base), use `[Accumulator; N]` indexed by a method like `Base::known_index()` rather than named fields — this eliminates match arms for invalid variants and works naturally with const arrays like `Base::KNOWN`.
5. To extract a single group's accumulator, use a `take(&mut self, key) -> Option<Accumulator>` method via `mem::take` — `None` signals "not applicable" (e.g. Unknown base) rather than an error.

## Read orientation modes

The main pileup-based `call` path assigns OT/OB in `src/call/pileup/from_hts.rs` before `PileupMetrics::new()` sees a `SimpleRead`.
The default `VariantCallingParams.read_orientation=flags` path still uses `strand_from_flags()`.
The opt-in `--guess-read-orientation` mode does **not** require reference CpG annotation. Instead it scans each read over `aligned_pairs_full()` and only looks at read positions where the observed base mismatches the reference:

- at each mismatch, inspect both 2 bp windows that include that read base: current+next and previous+current
- count `TG` motifs and `CA` motifs in the htslib/reference-oriented read sequence
- `TG > CA` means OT, `CA > TG` means OB
- ties / no evidence: split deterministically from a hash of qname + start + flags so repeated runs stay reproducible

Implementation detail: `src/call/process/pileups.rs` keeps a per-segment `ReadOrientationCache`, because `alignment_to_read()` is called once per pileup column and mismatch scoring would otherwise rescan the full read for every covered base.

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

## Indel parity between the two pileup backends

`from_hts.rs` and `from_seqair.rs` build the same `PileupMetrics` from two
different readers. **SNVs are byte-identical across them** — F1, recall,
precision and both aardvark genotype-error counts — which is the control that
makes any remaining difference attributable to indel-specific code rather than
to read handling. Indels were *not*: on chr12 the seqair backend scored 4-5 F1
points below htslib, from five separate divergences, all now fixed.

**How the divergence is measured.** Two numbers, and the second is the more
useful one:

```bash
CARGO_TARGET_DIR=target-hts cargo build --release          # htslib
cargo build --release --features experimental-seqair       # seqair
# ...call chr12 with --experimental-indels=ml --vcf-all-fields on each, then:
bcftools view -f PASS -i 'GT="alt"' calls.bcf -Ou \
  | bcftools norm -f hg38.fa.gz -m -any -Oz -o q.vcf.gz
aardvark compare --reference hg38.fa.gz --truth-vcf truth.vcf.gz \
  --truth-sample NA12878 --query-vcf q.vcf.gz --query-sample sample \
  --regions PG_ConfidentRegions_hg38.bed.gz --output-dir out/
# and, far more sensitive than F1, the site-level disagreement:
bcftools isec -p isec <(bcftools view -v indels hts.vcf.gz) <(bcftools view -v indels sq.vcf.gz)
```

The `isec` counts move ~100x over the fixes where F1 moves 5 points, so use them
to tell "this changed something" from "this changed the right thing". Aardvark
needs `Number=.` in the header (see above) or it refuses the file outright.

**The five divergences, in the order they were found.** Each is worth knowing
because each is a *class* of mistake, not a typo:

1. **A read-local offset used as a genomic position.** `view.qpos()` shadowed
   the column's `pos`, so deletion REF alleles were read from the segment's
   opening bases. Guarded now by seqair's `QPos` newtype, which is why the pin
   carries it.
2. **The tract anchor.** `homopolymer_run_at`/`dinucleotide_run_at` must be
   measured one base past the pileup anchor, where a left-aligned indel starts.
   The convention now lives only inside `ref_features::indel_tract_runs_at`.
3. **Overlap dedup dropping the fragment's only indel.** The rule keeps one
   read per fragment by base agreement and template order, which can keep the
   mate that does *not* span the indel. `from_hts` never had this because it
   votes per fragment *before* deduplicating — every mate gets a turn and the
   first surviving observation is the fragment's. **The single largest one:
   +2.3 insertion / +1.5 deletion F1.**

   **The fallback is keyed on the observation, not on the indel.**
   `build_indel_observation` rejects an indel too close to a read end or on a
   read with too many non-TAPS mismatches, so a kept read can carry an indel
   and still contribute nothing, and the mate must then stand in.
   `PileupColumn::pair_indel` (seqair) cannot express that — its rule is "own
   indel wins, the mate is never consulted", which is correct for a query that
   only sees CIGARs but silently drops those fragments. **rastair therefore
   does not use `pair_indel`**; `counted_mate` plus `Option::or_else` is the
   whole mechanism, and it needs no seqair query.
4. **Two implementations of one predicate.** `has_repeat` used a 4-base window
   for the period-2 arm on one side and 6 on the other, so it fired ~8x more
   often on the seqair path. There is now one `has_terminal_repeat`.
5. **Per-read where the semantics are per-fragment.** `soft_clip_count` and the
   noisy-reference count describe a fragment; `from_hts` ORs them over both
   mates. Reading them off the surviving read alone loses what only the dropped
   mate showed. The noisy-reference count also has to exclude a fragment that
   contributed an observation — `IndelCounts::clean_depth` subtracts a fragment
   counted on both sides twice — and `AlignmentShape::noisy()` is
   `terminal_repeat || soft_clipped`, not the repeat alone.

**Result on chr12 (~26x, bundled model, `--experimental-indels=ml`), against
Platinum Genomes in `PG_ConfidentRegions`:**

| | seqair before | seqair after | htslib |
| --- | --- | --- | --- |
| SNV F1 | 0.9706 | 0.9706 | 0.9706 |
| Insertion F1 | 0.7775 | **0.8273** | 0.8274 |
| Deletion F1 | 0.8099 | **0.8506** | 0.8506 |
| disagreeing indel sites | 29,420 | **5** | — (32,326 shared) |

Deletions come out numerically identical to htslib on every column — recall,
precision, F1, `truth_fn_gt` and `query_fp_gt`. Insertions differ by two sites
in `truth_fn_gt` and nothing else.

**Judge a parity fix by convergence in every column, not by F1.** The last fix
moved deletion F1 *down* 0.0010 while moving recall, precision and both
genotype-error counts onto htslib — which is what says the semantics matched
rather than a threshold moving. A change that improves F1 while moving
`query_fp_gt` away from the reference has not fixed parity.

The residual is 5 sites of 32,331 (4 htslib-only, 1 seqair-only), all
multi-allelic columns in homopolymer or short-tandem-repeat tracts where the two
readers group one read's allele differently (`G>GTT` alongside `G>GTTT` at one
position, `GTTTT>G`, `CAAA>C`). Not worth chasing without a reason to.

### Why none of this was caught

Worth internalising, because the same blind spots are still easy to reproduce:

- **No CLI test enabled indel calling.** `grep experimental.indels
  tests/call_cli.rs` was empty and every snapshot header recorded
  `"experimental_indels":null`, so no snapshot ever held a multi-base REF.
  `indel_ref_alleles_come_from_the_deletion_site` now closes that.
- **`tests/data/test.bam` cannot produce an indel call.** It has 86
  indel-carrying reads, but no two agree on an allele, so nothing passes
  `--min-indel-ao` at any threshold. An indel test has to build its own BAM.
- **A fixture segment starting at 0 hides every position bug**, because
  `pos - segment_start` and the clamped wrong answer coincide there. Use
  `segment_at(start, seq)` with a non-zero start; the CLI fixture puts its
  deletion at 24,137 for the same reason.
- **The `from_seqair` unit tests never set `params.call_indels`**, so the whole
  indel branch was unreachable from them.


## Phasing (`--phase`, seqair only)

`src/call/phasing.rs` and its `problem` and `solver` modules. Per-column data sits
in `PileupMetrics.phase`; thresholds are constants, there are no `--phase-*` flags. `--phase` is
an error on htslib builds and under CpG-only output (`--cpgs-only`, or `--bed` alone): the
pre-filter drops the variant sites before phasing runs.

### Observations

- Captured inside the seqair column loop (`PileupMetrics` keeps no reads). A column keeps its
  whole observation list when it has `MIN_ALT_READS` mismatches on one base that TAPS cannot
  explain (`candidate_observations`). An `N` reference or an empty column is
  excluded, so `phase.is_some()` always means there is something to use.
- `accepts_read` (MAPQ and base quality ≥ 20) is applied at collection; a read can count for
  depth and genotype without carrying phase.
- **`FragmentId` is `NonZeroU64`**: seqair reports "no qname" as a missing hash, and bucketing
  those reads together would link every het in the region. They are dropped, and
  `call` warns once per run: a missing qname is a whole-file property (CRAM with `RN=false`),
  so there is no count worth reporting.

### Two TAPS gates — do not confuse them

- `taps_explains_mismatch` is *mismatch-level* (ref C + OT + T, ref G + OB + A), used when a
  column is built. It ignores sequence context: de-novo CpGs are not known yet.
- `PhaseSite::taps_confounds` is *allele-set-level* and **CpG-only**: at a C/T site whose `C`
  is in a CpG (`methylation::cpg_origin`, reference or de-novo), a `T` is evidence only on OB;
  at G/A an `A` only on OT. Outside a CpG every strand counts, since mammalian non-CpG
  methylation is rare (neurons and ES cells are the exception). Both gates use
  `methylation::CpgSide`; do not add a second copy of its bases and strand.
- Fixture trap: at a C/T site in a CpG an OT `T` is confounded, at A/G an OB `A` is; a bare
  `het()` is not in a CpG, so set `pos_metrics.cpg` or `context.after_1`. Tests of something
  other than the gate should use A/C or T/G sites. Both mates share the OT/OB assignment.

### Problem and solver

- **A site's alleles come from the stored `GenotypeTag`, never from scanning `pileup.alts`.**
  `GtAllele { Ref, Alt(NonZeroU8) }` is in *genotype space* (index into `pileup.alts`), not
  `ALT` column order; `compute_genotype` remaps both the genotype and the phase call with one
  map. Do not phase after remapping.
- Fragments reaching fewer than two sites are dropped; fragment order is first appearance, so
  the result does not depend on the hasher.
- **A zero-weight edge is not an edge**: cancelling cis and trans evidence would otherwise emit
  a coin toss as a confident phase.
- `phase_segment` runs **before** the overlap trim, because a site in the overlap links reads
  reaching into the core, and takes the segment's `is_core` predicate: each block is anchored on
  its first *emitted* site, so `PS` names a record of this file and is unique across segments,
  and a block with fewer than two emitted sites is not written. The solver itself returns raw
  orientations, defined only up to a flip per block; `phase_segment` is the one place that
  normalises, onto the anchor, so the first record is `0|1`.
- `solver::refine` is a greedy descent run to convergence, capped at `MAX_FLIPS_PER_SITE` flips
  per site of the block; hitting the cap is a `warn!`.
- `PS` is written only where the genotype is phased. Indel records and homozygotes are never
  phased.

### Testing phasing

- `tests/phasing_cli.rs` runs `call --phase` on a `Scenario` (`tests/utils/scenario.rs`): a
  pseudo-random reference with no CpG unless the test writes one (`reference_at`), het SNVs per
  haplotype (`snv`, `het`), per-haplotype methylation (every CpG on that haplotype, de-novo ones
  included, converted on the read's own strand), and read pairs per haplotype and strand
  (`pairs`, `balanced`). `check_phasing` reads the VCF back and fails on a switch error, a phased
  non-het, a `PS` not naming its block's first record, or a block not opening on `0|1`; tests
  then assert which sites share a block. Read the observations through the real seqair path
  this way rather than attaching them to a `pileups!` fixture.
- Make a scenario *decisive*: arrange the reads so the broken code gives a different block
  structure, not just a weaker edge. The TAPS gate tests do this by having converted reads
  exactly cancel the honest ones; the boundary test (`a_block_is_solved_across_…`) links two
  sites only through a third that sits in the other segment's overlap.
- `phasing_accuracy_on_an_na12878_slice` scores 200 kb of the chr12 NA12878 TAPS BAM against
  Platinum Genomes and pins phased sites, block pairs, switch errors and phased non-hets just outside what was
  measured. Its scorer agrees with `whatshap compare` on pairs and switches. It is the only
  test that sees phasing get *worse*: the pre-CpG gate loses 8 phased sites there, and no gate
  at all makes 18 switch errors. Refinement and the phase quality gate change nothing on it.
  The reads are Watchmaker data and stay out of git (Pascal, 2026-10-06): build the fixture
  with `scripts/make_phasing_slice.py tmp/taps tmp/na12878_phasing` and point
  `RASTAIR_PHASING_SLICE` at it, or the test skips. Run it before merging phasing changes.
- Check a new phasing test is not tautological by breaking the code it guards and watching it
  fail; every test added on 2026-10-06 was checked that way.
- `tests/data/test.bam` is a poor phasing fixture: its few het calls are C>T on a fully
  methylated control. Use a `Scenario` for anything that asserts phasing output.
- A run that reports zero phased sites is usually an htslib binary (see "Know which binary
  produced a result"), not a regression.

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

## Measuring accuracy

Nothing in-tree computes F1 or switch error; use the reference tools, and never judge a change
on a slice alone for genotype-error counts (they are too rare to resolve on 10 Mb).

```bash
# Genotype accuracy, haplotype-aware. Normalise the query first.
bcftools view -f PASS -i 'GT="alt"' calls.bcf -Ou | bcftools norm -f hg38.fa.gz -m -any -Oz -o q.vcf.gz
aardvark compare --reference hg38.fa.gz --truth-vcf truth.vcf.gz --truth-sample NA12878 \
  --query-vcf q.vcf.gz --query-sample sample --regions confident.bed.gz --output-dir out/
# Phasing: phased fraction, blocks, switch/flip errors (sample names must match: bcftools reheader -s)
uvx --from whatshap whatshap stats --tsv=stats.tsv calls.vcf.gz
uvx --from whatshap whatshap compare --only-snvs --names truth,rastair --tsv-pairwise p.tsv truth.vcf.gz calls.vcf.gz
```

- aardvark's `summary.tsv` splits zygosity errors out: `truth_fn_gt` (truth hom-alt, query het)
  and `query_fp_gt` (query hom, truth het). It ignores query phasing by design.
- Do not join aardvark's output VCFs on `CHROM:POS:REF:ALT`: it rewrites records. Join the
  original truth and query, both `bcftools norm -f <ref> -m -any`.
- Judge a change by every column (recall, precision, both genotype-error counts), not by F1,
  and check that an indel change leaves the SNV subset identical with `bcftools isec` of
  normalised SNVs (it compares sites, not genotypes, which also makes it the tool for
  stratifying errors).
- Do not mix depth sources when stratifying: `INFO/DP` is after overlap dedup.
- Truth data for NA12878 (Platinum Genomes phased VCF, its confident regions) lives outside git;
  GIAB HG001 has no phased records and cannot score phasing.
- A killed run can leave a smaller BCF that reads cleanly; check the last position or
  `bcftools index -n` before trusting a whole-chromosome number.

## Know which binary produced a result

`cargo build --release` without `--features experimental-seqair` replaces
`target/release/rastair` with the htslib build, and a plain `cargo test` rebuilds
`target/debug/rastair` the same way. Keep one `CARGO_TARGET_DIR` per backend when comparing,
and check the log: `grep -c "Using experimental seqair backend" run.log` is 1 on seqair.
The backends differ in mate dedup and soft-clip handling, so arms built on different
backends are not comparable. Likewise score every arm of a comparison with `--gpu` or every arm
without it: the two are not bitwise identical.

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
