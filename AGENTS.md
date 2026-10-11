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
- Run the tests with `cargo test` **and** `cargo test --features experimental-seqair`. Phasing and indel calling are seqair-only, so their tests are `#[cfg(feature = "experimental-seqair")]` and a plain `cargo test` never runs them; the htslib build instead tests that `--phase` and `--experimental-indels` are refused.
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

## Indel observations (seqair only)

Only `from_seqair.rs` collects indels; an htslib build refuses `--experimental-indels`. The rules
below came out of a parity hunt against the htslib backend's indel collection, since deleted (its
write-up is in this file's history before `feature/pair-hmm-cleanup`). Each guards a *class* of
mistake:

- **A deletion's REF comes from the deletion site.** `view.qpos()` once shadowed the column's
  `pos`, so REF alleles were read from the segment's opening bases; seqair's `QPos` newtype now
  keeps a read offset from passing as a genomic position.
- **The mate fallback.** Overlap dedup keeps one read per fragment and can keep the mate that
  does *not* span the indel (worth +2.3 insertion / +1.5 deletion F1 when fixed), so the mate
  stands in: `build_indel_observation(view).or_else(mate)`. seqair's old `pair_indel` ("own
  indel wins, the mate is never consulted") is not used.
- **The CIGAR counts only pick candidates.** `IndelCounts` holds, per allele, the fragments whose
  alignment carries it; nothing else reads it for a call. There are no read-level indel filters
  any more (end-of-read cutoff, non-TAPS mismatch limit, the noisy-fragment split behind "clean"
  support and depth, `has_terminal_repeat`): they were a second evidence pipeline beside the
  pair-HMM, which weighs each read's alignment itself. Measured on chr12 when they went: the
  mismatch filter alone moved nothing past 0.001 F1; ranking candidates by raw rather than clean
  count was +0.006 ins / +0.004 del at 26x and -0.002 / -0.004 at 73x (more scored alleles in
  repeats, more deletion FPs).

Testing blind spots, still easy to reproduce:

- **`tests/data/test.bam` cannot produce an indel call.** It has 86 indel-carrying reads, but no
  two agree on an allele, so nothing passes `--min-indel-ao`. An indel test has to build its own
  BAM (`indel_ref_alleles_come_from_the_deletion_site`) or use a `Scenario`.
- **A fixture segment starting at 0 hides every position bug**, because `pos - segment_start`
  and the clamped wrong answer coincide there. Use `segment_at(start, seq)` with a non-zero
  start; the CLI fixture puts its deletion at 24,137 for the same reason.
- **A `from_seqair` unit test must set `params.call_indels`**, or the indel branch is
  unreachable from it.

## Pair-HMM indel calling (`--experimental-indels`, seqair only)

`--experimental-indels` calls indels with the pair-HMM and nothing else: no model, no binomial
genotyper, no strand veto. It is seqair-only (htslib refuses it, like `--phase`) because it needs
the reads, which only the seqair engine still holds once a sub-segment's columns are drafted.

**Dependencies.** seqair 0.5 and compair 0.2 (both 2026-10-08) carry what this needs:
`Haplotype::with_insertion`/`with_deletion`, `qname_hash` as `Option<NonZeroU64>` (adopted
through `FragmentId::of`), `leading_soft_clip`/`trailing_soft_clip` and `fragment_idx` on
records (`PlacedRead`). Gaps are scored by compair's `PcrIndelModel` (GATK's PCR indel model,
via `Read::with_pcr_indel_model`) with the tandem-repeat gap-open floored at Q20. compair also
has `Haplotype::with_betas`/`Betas::OfHaplotype` and `ScoreMatrix`, not yet used here (see
`UNKNOWN_BETA` for why the first stays off).

- **Deferred scoring** (`process/pileups.rs`). Indel columns are not scored while they are
  built: the kept ones are collected, and after the sub-segment's column loop, columns within
  `CLUSTER_GAP` (20 bp) of the previous one form a cluster (chained, so a cluster can be longer).
  Each cluster is scored once, `from_seqair::score_indel_cluster`, through the sub-segment's
  `ClusterContext` (store, segment, params and window mask, built once), over the reads
  `guard.records_overlapping(first, last)` returns; the engine's store stays valid until it is
  reclaimed. Each column of the cluster shares the cluster's `Arc<HaplotypeScores>` (`Scored`)
  and gets its own records. Deferring alone left every record byte-identical; it is what makes
  the window mask and clusters possible.
- **`WindowMask`**: positions of the sub-segment's kept columns with at least `MIN_READS` (5)
  reads of a base TAPS cannot explain are `Base::Unknown` in the haplotype window, which compair
  matches to any base. A `T` over a `C` counts only on OB and an `A` over a `G` only on OT, so
  methylation never masks. Without it a nearby SNV absent from every haplotype is charged to the
  reads of the haplotype that carries it, and a stray indel allele that happens to absorb the
  mismatch gets the credit: it cut the 73x deletion FPs next to an SNV from 265 to 5 and cost
  nothing elsewhere. A count, not a fraction, since the FP it guards against is a handful of
  reads whatever the depth; 2, 3, 5 and 8 were tried, 5 was best at 26x and 73x.
- **Events** (`cluster_candidates`). A cluster's `Event { anchor, allele }`s are the alleles of
  its columns with at least `--min-indel-ao` reads, ranked by read count, then anchor, then
  `allele_order` (insertions first), as many as `Cluster::admits` (which `Cluster::new` enforces
  too): at most `MAX_EVENTS` (3) and **at most two per anchor**: a third at one anchor could not
  share a pair with the other two, so it would take a slot from a neighbouring anchor and never
  be called. That ranking is record order: `Cluster::events`, `EventIdx`, the calls and the
  `ALT` column all keep it. An allele not ranked in is never scored,
  so never called.
- **Haplotypes are subsets of events** (`pair_hmm.rs`). `HapIdx` is a bitmask over `EventIdx`;
  `Cluster::new` builds them once, `HapIdx::all_valid` keeping the subsets one chromosome can
  carry (`Event::compatible`: not two at one anchor, no deletion over another event's anchor),
  so at most eight haplotypes and eight alignments per read. `scoring::score` builds one compair
  `Haplotype` per subset, applying its events later anchor first so an earlier anchor's window
  index still holds; an event that does not fit the window is dropped with a `debug!`. Every CpG
  is scored at a uniform beta 0.5 (`UNKNOWN_BETA`): the columns' own level, raw or shrunk
  towards 0.5 with 2 to 16 reads of prior, measured worse at both depths. The reads are not the
  column's counted ones: one per fragment (the mate reaching further), MAPQ >= the column
  minimum, a known strand, soft clips included where the segment's reference reaches them,
  thinned evenly to `--max-coverage`.
  Stored per read on `HaplotypeScores.reads`: fragment, MAPQ, log10 likelihood per `HapIdx`
  (`-inf` for an invalid subset).
- **The verdict** is computed once, in `HaplotypeScores::new`, over pairs of valid haplotypes.
  An account (`Account`) is a pair, or a pair plus one event it lacks drawn as an artefact at
  `ARTEFACT_FRACTIONS` (10 %) less `ARTEFACT_PRIOR_LOG10` (2, chosen on chr12 26x+73x,
  docs/plans/data/artefact-sweep-2026-10-08): that is what stops a 7-of-73 allele at 73x being a
  genotype, while a 4-of-8 het survives. Margins are **per event**: log10 of the best account
  whose pair carries the event over the best whose pair does not. An event is called at
  `Margin::PRESENT` (>= 1); QUAL is 10 x the margin. `best` is the pair of the best account; `called_pair` is the
  likeliest pair carrying exactly the called events and *is* the genotype, so "present but not
  in GT" cannot happen (such a pair exists while no anchor has more than two events). A `1/2`
  whose second event is not called becomes `1/1` or `0/1` as the reads have it. `quality`
  (`GC`) is the called pair's flat-prior quality over the pure pairs, because the phasing gate
  was set on it, so it is low when the called pair is not the best one.
- **Calling** (`call_on_pair_hmm(scores, pos)`) writes the records of the events anchored at
  `pos` and no other: two called events at one anchor make one `IndelRecords::Compound` `1/2`
  record (always PASS, QUAL the weaker margin's); otherwise one `Separate` biallelic record per
  event, its genotype from `Verdict::copies_of` (0, 1 or 2 chromosomes of the called pair). A
  scored event not called stays a `0/0` call, written only under `--all` with `low_indel_qual`
  (read off the genotype, not the margin). No strand gate and no model: neither beat the margin
  alone on chr12. Emission and `phasing::problem::indel_site` match on `IndelRecords` rather
  than re-deriving it from the genotypes.
- **`AD` and `DP` are the pair-HMM's**, so every field of an indel record comes from one read
  set: a read is assigned to its likeliest haplotype when that beats the runner-up by
  `INFORMATIVE_LOG10` (0.2, GATK's informative-read threshold), once, in `HaplotypeScores::new`,
  which keeps the count per haplotype. `DP` is `HaplotypeScores::depth`, the cluster's assigned
  reads, and a record's `AD` is
  `support_over(&events)`: reads on a haplotype carrying none of the record's events (reference
  first, so a neighbouring event's carriers count there), then on one carrying each. `AD` sums to
  `DP`. The threshold is far below `Margin::PRESENT` on purpose: the pair-HMM sums a tract gap
  over every placement, so one read of a 1 bp deletion in a 13-base homopolymer is only ~0.9
  log10 likelier under the deletion, and a read-level "ten times" assigned no read at exactly
  those indels (24 % of true het deletions at 73x had `AD` 0 for the allele).

Two mechanisms explain most of what clusters changed, and are worth knowing before touching them:

- **Two-event clusters.** `AT>A` at chr12:50,733,241 and `CTT>C` at 50,733,254, both homozygous,
  13 bp apart: with one event per haplotype the true haplotype (both deletions) was not in the
  set, so one of them lost its reads to the other's mismatches; the same for an insertion and a
  deletion 9 bp apart and for poly-A tracts with three spellings plus a real event. Subsets of
  events represent these exactly, which won back recall clustering had cost (26x deletion
  recall .9085 -> .9145; 73x homozygotes written as hets 39 -> 31).
- **TAPS twin spellings.** A deletion next to a methylated CpG has two spellings whose
  haplotypes differ by exactly a `C`/`T` at the CpG (`CAT>C` at 1,376,406 against `CCA>C` at
  1,376,403). Converted OT reads match the wrong spelling for free and the right one at the
  conversion cost, so at beta 0.5 the wrong one sometimes collects them. Per-column scoring
  called both (one FP), a cluster picks one. It is TAPS-specific, and why left-aligning the
  candidates did nothing: the aligner already left-aligns, these are different events. The
  per-column beta was meant to break the tie and measured worse; open.

chr12 NA12878 (per-column = one column's two leading alleles, before clusters):

| | 26x PG Ins / Del F1 | 73x GIAB Ins / Del F1 |
| --- | --- | --- |
| forest (`main`) | .827 / .851 | .946 / .953 |
| per-column pair-HMM | .8944 / .9102 | .9525 / .9497 |
| clusters of event subsets | **.9050 / .9230** | **.9705 / .9798** |

Hold-out chr17 (26x, 73x) and chr2/chr20 (73x) agree, with 73x deletion precision back at the
forest's; tables in `docs/plans/indel-cluster-prototype-2026-10-08.md`. SNV records are
byte-identical across every arm. The constants above were each chosen on chr12 at both depths;
the reasons sit next to them in the code.

Testing: the cluster pieces have unit tests (`cluster_candidates`, `WindowMask`,
`a_masked_position_matches_any_base`, `two_deletions_on_one_chromosome_are_one_haplotype`,
`call_on_pair_hmm` with three events) and one CLI scenario,
`the_pair_hmm_calls_two_events_of_one_cluster`, which fails with one event per haplotype.

### The pair-HMM on the GPU (`--gpu`)

`--gpu` also scores the pair-HMM on the GPU (compair's `gpu` feature, part of
`experimental-seqair`): the `gpu` module of `from_seqair.rs`. The calls are identical to the
CPU path by construction — the kernel is bit-identical to compair's CPU kernels, and a pair
whose `f32` score falls under its precision floor is rescored in `f64` on the CPU
(`GpuScore::or_rescore`, ~0.8 % of pairs) — and were checked byte-identical on chr12 at 26x
and 73x. What to know when touching it:

- **`score` is split in two.** `pair_hmm::prepare` builds the cluster, haplotypes, owned
  `HaplotypeRead`s and bands; `Prepared::align` runs the CPU kernel, `Prepared::push` adds the
  cluster's pairs to a `GpuPairs` plan and `Prepared::finish_gpu` turns a launch's scores back
  into `HaplotypeScores` through the same `rows`. Keep the two arms sharing everything but the
  kernel call, or the identity check is the only thing standing between them.
- **One GPU thread for the run** (`pair-hmm-gpu`, like the ML's `InferenceStage`): a worker
  sends its plan through its `Lane` (a `thread_local` reply channel), keeps preparing clusters,
  and collects the scores later. The thread merges every plan queued while a launch ran into the
  next one (`GpuPairs::append`), so one aligner's buffers serve all workers. Being the only
  thread polling its device is also what makes `Wait::Block` safe on compair 0.2.0, whose
  blocking wait loses a map callback when several threads poll one device (~1 launch in 300 with
  the earlier one-thread-per-worker design; then every launch fails with "Buffer with 'scores'
  label is still mapped"). A failed launch rebuilds the aligner and fails every plan in it over
  to the CPU.
- **A worker's plan is at most `MAX_LAUNCH_PAIRS` (4,096) pairs**, so a sub-segment's clusters
  go up in chunks and the prepared reads a worker holds are bounded — two plans' worth, the one
  on the GPU and the one filling (a single cluster can exceed the cap; its plan holds it alone).
  The per-region `Batch` (in `ColumnScratch`) is flushed after each sub-segment's cluster loop,
  before the next sub-segment's reads are loaded: holding a launch across that load cost
  100 MB at 73x for no wall time. A batch dropped with a launch in flight drains it in `Drop`,
  or the next region would collect it.
- Every failure (no adapter, a plan too large for the device, a lost thread, a failed launch)
  falls back to the CPU kernel for the clusters concerned, with a `warn!`. Without `--gpu` the
  module is never entered, and the plain `cargo test` never exercises it.
- `RASTAIR_LOG=rastair::call::pileup::from_seqair=debug` logs one line per launch with running
  totals (launches, pairs, rescored); `grep -c "Pair-HMM GPU launch"` is the check that the GPU
  path ran (one per plan; the GPU thread logs "one launch for queued plans" with how many it
  merged). Measured only on Apple Metal (M4 Pro) so far; a Vulkan identity check is open.

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
- `PS` is written only where the genotype is phased. Homozygotes are never phased; indel records
  are phased as the next section says.

### Indel sites

- A column can hold an SNV site and an indel site (`SiteKind`); `ColumnPhase::Phased` holds a
  `RecordPhases` with one optional phase per kind. At a shared position the SNV record, pushed
  first, anchors the block. `PhaseCall` carries only the block's `set` and `first`, the allele
  haplotype 1 carries; the other haplotype's allele is the genotype's other one.
- `problem::indel_site`: the pair-HMM's called pair must be a het (`0/1` or the compound) and *be*
  its best pair. The site's two haplotypes are the called pair's own, `HapPair::split(event)`
  (carrier, other), so a read votes on them with the cluster's other events included: the
  carrier of a `0/1` may also carry a neighbouring event, and the haplotype of the event alone
  may explain no read at all.
  A read votes by the sign of its log10 ratio, cost `round(10 * |ratio|)` capped (inf is the cap,
  0 no vote); no base quality or TAPS gate, the emission has weighed both.
- **The gate**: GQ >= 100 (unclamped) and >= 5 *votes* per allele at Phred 10 or more
  (`Cost::of_margin(Margin::PRESENT)`), counted after `scored_vote` has dropped reads without a
  fragment or below the phasing MAPQ, so a site is admitted on what it feeds the solver. A
  wrong indel bridge is an SNV switch error, so judge a gate on SNV-only switches (whatshap
  `--only-snvs`) and phased non-hets, not on indel counts. Chosen on chr12 NA12878 (~26x) only:

  | gate | phased indels | SNV switches | phased non-hets |
  | --- | --- | --- | --- |
  | none (SNV-only phasing) | 0 | 66 | 612 |
  | **GQ >= 100, >= 5 reads** | 5,431 | 67 | 886 |
  | GQ >= 100 alone | 6,768 | 71 | 979 |
  | GQ >= 40 alone | 9,622 | 80 | 1,217 |

  Tried and not worth it: a margin floor, tract-length caps. The remaining bad bridges are
  tetranucleotide unit expansions and TAPS-converted twin candidates near a CpG, a candidate
  generation problem.

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
- A `Scenario` can hold an `insertion` or `deletion` per haplotype (on both for a homozygote,
  two different ones at one anchor for a `1/2`); reads crossing it carry it in their CIGAR, and
  `check_phasing` judges indel records too.
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
- Model structs: `CpgFeatures` (55), `DenovoCpgFeatures` (56), `OthersFeatures` (54). Each has
  an `extract()` returning the struct; `FeatureCalculator::calculate_*` wraps `as_row()` into an
  `Array2`. Indels have no model: the pair-HMM calls them.

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
  the recipe of the bundled model (`models/rastair_default.rff.mpk.lz4`, `include_bytes!`'d
  by `src/call/ml.rs`), so a bare command line retrains it. msgpack writes `RastairFlatModel` as
  an array, so its field order is the file layout.
- The truth VCF needs a **`.csi`** (a `.tbi` alone is refused), and `-R` must be a plain-text
  BED, not bgzipped. Without `-R`, candidates outside the truth set's confident regions are
  labelled negative.
- Train on GIAB HG001 and evaluate against Platinum Genomes, holding out the evaluation
  chromosome. **Give it enough genome**: all three forests must train or no model is written.
- Collection is bounded: rows are `f32` in one flat buffer, each model keeps per class the
  examples with the smallest uniform keys (a keyed reservoir, capped at the draw plus a
  holdout), and segments are folded with `reduce`, not collected. `KeySource` seeds one stream
  per segment *and* model from `--seed`, so the sample does not depend on scheduling, and one
  model's pool does not move another's draw.
- The draw takes at most four fifths of a class, so small models keep a holdout; the holdout
  keeps the **population's** class ratio, not the reservoir's near-balanced one, or Platt
  scaling would shift what an ML threshold means.
- Every model draws `--n-positive`/`--n-negative`. The log line `Collected training examples`
  reports pool, kept, requested and drawn per model.
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
