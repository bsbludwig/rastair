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
- **One GPU thread per worker** (`Lane`, a `thread_local`): the worker sends a plan, keeps
  preparing clusters, and collects the scores later. The aligner's buffers live on that thread
  and grow to the largest launch once, not per region. The wait is compair's `Wait::Spin`, which
  costs the helper thread a core while a launch runs: `Wait::Block` would sleep instead, but in
  compair 0.2.0 it misses the mapping callback on ~1 launch in 300 (8 % of launches with 4k-pair
  launches at 73x), and every launch after that one fails with "Buffer with 'scores' label is
  still mapped". The fix is on seqair branch `fix/compair-block-wait`; switch the one line once a
  compair with it is published. A failed launch rebuilds the lane's aligner, so one failure does
  not take the lane down for the run.
- **A launch is at most `MAX_LAUNCH_PAIRS` (4,096) pairs**, so a sub-segment's clusters go up
  in chunks and the prepared reads a worker holds are bounded — two launches' worth, the one on
  the GPU and the one filling (a single cluster can exceed the cap; it launches alone). The
  per-region `Batch` (in `ColumnScratch`) is flushed after each sub-segment's cluster loop,
  before the next sub-segment's reads are loaded: holding a launch across that load cost
  100 MB at 73x for no wall time. A batch dropped with a launch in flight drains it in `Drop`,
  or the next region would collect it.
- Every failure (no adapter, a plan too large for the device, a lost thread, a failed launch)
  falls back to the CPU kernel for the clusters concerned, with a `warn!`. Without `--gpu` the
  module is never entered, and the plain `cargo test` never exercises it.
- `RASTAIR_LOG=rastair::call::pileup::from_seqair=debug` logs one line per launch with running
  totals (launches, pairs, rescored); `grep -c "Pair-HMM GPU launch"` is the check that the GPU
  path ran. Measured only on Apple Metal (M4 Pro) so far; a Vulkan identity check is open.

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
- `PS` is written only where the genotype is phased. Indel records and homozygotes are never
  phased.
- Phase calls carry only `PS` and the allele on haplotype 1, with no phase probability or
  confidence. `GL`/`GC` remain genotype scores. Observation costs are base qualities capped at
  40; pair weights sum signed minimum costs, and every nonzero weight can link a block.
  The first emitted site's lower genotype allele anchors the orientation, so a compound-het
  anchor can open with `1|2` rather than `0|1`.
- Full reads remain accessible through seqair's pileup guard `store()` until the guard is
  dropped; `get_pileups` drops each guard at the end of its byte-budget sub-segment. Capture
  or score any full-read evidence there, before recycling, rather than trying to recover it
  from `SimpleRead` or `ColumnPhase`. Pre-ML haplotype evidence also needs its own candidate
  CpG context: `methylation::cpg_origin` requires called variants and cannot decide whether
  an uncalled candidate haplotype creates or destroys a CpG.

### Production likelihood phasing (2026-10-14)

`call --phase` now calls `haplotypes::phase::score` on **called**, genotype-conditioned SNP blocks. The older `problem`/`solver` and `legacy_phase_segment` are `cfg(test)` reference code; production does not use signed-quality weights or MEC refinement. Earlier collection/solver notes above describe that reference path.

`get_pileups_with_phase` captures full records using `haplotypes/reads.rs` before each seqair guard is recycled and returns segment-level fragments separately from `PileupMetrics`. Deduplicate alignments fetched at byte-budget boundaries by read group, name, start and flags; molecule grouping includes read group. Missing names are excluded with the existing once-per-run warning. Do not recover full sequences from `ColumnPhase`.

Genetic edges require two independent fragments and two unconfounded observations for each called allele (consensus BQ/source MQ >=20). Converted CpG T/A cannot establish an edge. The score uses all eligible fragments touching the selected region, including one-site and confounded reads. Adjacent candidate C/G partners are included as context dependencies. Apply homozygous SNPs to the reference background before constructing templates. Genotype alleles still use `GtAllele` space before VCF remapping.

Regions expand along connected sites, limited to 32 sites/2 kb. Phase-state beam 32, holdout 16, with competing site orientations protected up to 64 retained states. Only complementary chromosome pairs matching the called genotype are allowed; one supported nearby nuisance indel (<=20 bp) is marginalized. Shared CpG shrinkage defaults are 0.95/0.003, Beta(2,2), leave-one-molecule-out. Edges need retained-state phase posterior >=0.99. Orient and anchor on the first emitted genomic site even when expansion order differs from position order; core overlap semantics and PS-only-for-hets remain unchanged. These conditional posteriors are not calibrated PQ, and GL/GC remain genotype scores.

PairHMM production profiles use mate-local coordinate windows with 12 bp flanks. An unsequenced gap between mates is absent from each alignment template: likelihood-reuse keys must omit variants there too, while retaining +/-1 CpG partners. The whole-block key caused a large slowdown in dense regions because identical molecular hypotheses were rescored. Concordant overlaps are one profile; discordant CIGARs/refskips fall back to aligned scoring. Cache path normalization by row/column lengths **and gap parameters**. Scaled linear forward and the test-only log-domain oracle agree for short and 500-base profiles; both ordinary and read-local molecular probabilities normalize. The standalone diagnostic retains its existing full selected-window realignment geometry.

`tests/phasing_cli.rs` additionally tests a ten-site chain that forces pruning and a one-byte budget that repeatedly fetches the same boundary fragments. `haplotypes::phase` tests genotype-conditioned cis/trans posteriors, unlinked single-site molecules and confidence loss from disagreeing mates. Run both feature suites and clippy; a new compile warning on stderr breaks the `convert_cli::can_pipe_through` snapshot because it invokes cargo itself. Fix the warning rather than accepting that snapshot.

`scripts/benchmark_phasing.py` records matched CPU caller timing, RSS, phase statistics, binary hashes and complete non-phase VCF hashes. Caller timing excludes subsequent Python inspection. Archives and NA12878 results live in ignored `tmp/phasing-pair-hmm-benchmark/`; use the final binary/reports for conclusions, and do not treat stopped `.partial` output as complete. GIAB HG001 is unphased: unchanged unphased fields and synthetic phase oracles do not establish real-data switch-error accuracy.

`--phase-mask FILE` is repeatable and affects phasing only. `phasing::PhaseMask` reads plain/gzip/BGZF BED and the supplied UCSC `Chr Start End GapType` centromere table; intervals are zero-based, half-open and merged once before worker threads start. Each unmasked interval is an independent inference domain, with reference/context/PairHMM flanks clipped to its boundaries and nuisance indels required to fit inside. Keep variant and methylation calls in masked intervals. CLI tests assert masks split even fragment-bridged blocks and preserve every non-phase VCF field; the likelihood test checks excluded observations and indels cannot affect PairHMM or aligned fallback scores.

`scripts/assess_phasing.py` scores genotype-matched SNP phase in confident, unmasked intervals and reports genotype errors and phase coverage separately. NA12878 chr1p was assessed against official Platinum Genomes 2017 hg38 phased truth, not GIAB. The preserved full-chr1 prefix was emitted by the 9m36s progress check, before the run stalled in centromeric sequence; it was terminated on user request at 46m18s. File modification time is not time-to-first-emission. `full-chr1/chr1p_accuracy.json` and normalized WhatsHap comparisons agree: PairHMM 16/23,706 switches, legacy 59/29,681; matching phased SNPs 37,303 versus 44,282. The normalized comparison preserves called allele identity while allowing WhatsHap to assess compound-het/multi-ALT records as well. On the intersection of the two phase outputs (23,080 identical adjacent SNP pairs), PairHMM makes 15 switches versus legacy 24. Masked validation covers the former stalled centromere and three 1 Mb regions; every non-phase VCF field remains identical, and the 5 Mb centromere-spanning run completes in 40.51s with 16 threads. A complete masked-chr1 time has not yet been measured.

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
- When checking refinement fixtures, sum votes per site pair first. In
  `refinement_corrects_a_spanning_tree_mistake`, the two cis edges each total 60, ahead of the
  trans edge at 40, so the forest already gives the expected answer. Use the large-block and
  local-optimum property tests to exercise refinement. As of the 2026-10-13 checkout, the
  `phasing_accuracy_on_an_na12878_slice` test and `scripts/make_phasing_slice.py` described
  above are absent; those instructions do not provide an executable accuracy check here.
- `tests/data/test.bam` is a poor phasing fixture: its few het calls are C>T on a fully
  methylated control. Use a `Scenario` for anything that asserts phasing output.
- A run that reports zero phased sites is usually an htslib binary (see "Know which binary
  produced a result"), not a regression.

### PairHMM phasing performance (2026-10-14)

Superseded for the kernel on 2026-10-15: `forward.rs` is replaced by compair (see "compair is the
PairHMM kernel"). The `Scorer`, bounded-hypothesis skipping and parallel planning below still
apply; the bit-identity claims and `run_avx512` profiles describe the removed `f64` kernel.

`call --phase` spent ~90 % of its CPU in the scaled forward recurrence of
`pair_hmm.rs`, one closure-driven DP per (fragment, template, mapping
hypothesis), rescored at every beam step. On chr1:80000001-81000000 one 100 kb
segment with 257 het sites held 63 % of the whole run's samples while 15 workers
sat idle. Three changes, every one of them bit-identical to the old code:

- **`pair_hmm/forward.rs`** runs up to eight (profile, template) problems of
  one shape as lanes of `[f64; 8]`; every lane performs the scalar recurrence's
  operations in the scalar order, so a batch equals scoring its members one at
  a time. `run` dispatches to `#[target_feature(enable = "avx512f")]` /
  `"avx2"` copies of the same generic function (`is_x86_feature_detected!`),
  which agree bit for bit because no fused multiply-add is ever enabled and
  the row scale stays a division. The emission closure became a five-term dot
  product (`column_weights`): a column weights the A/C/G/T emissions and an
  insertion slot, each branch of the old closure maps to one or two nonzero
  weights, and adding a zero product or swapping two addends is exact.
  Batches run on the enclosing rayon pool when there is one
  (`rayon::current_thread_index()`), so the hot segment borrows the idle
  workers: `ordair`'s workers leave their loop once the segment queue is
  drained and are then stealable.
- **`Scorer`** (in `pair_hmm.rs`) takes `Plan`s — a fragment's jobs
  (profile rows × template weights) and the hypotheses that combine them —
  and batches every job of a `phase::score` matrix pass. It memoizes forward
  results by the *content* of the row and column terms (interned bit
  patterns, so a hit is an exact recomputation avoided), which is what makes
  the beam's repeated pilot passes cheap; shrinkage-pass results are looked
  up but not inserted, as the rates change every step. One `Scorer` per block
  in `phasing.rs`.
- **Bounded mapping hypotheses are skipped exactly.** Every alignment path
  emits each profile row once, through a match or an insertion, and a
  column's weights are a convex combination of the row's five terms, so a
  pass scores at most the sum of the rows' log-maxima (`ProfileTerms::bound`,
  template-independent; a mismapped mate's rows max at 0.25). A "mate
  mismapped" hypothesis (weight ~1e-6 for MAPQ 60) is skipped when
  `log_add(result_so_far, bound + 1e-3) == result_so_far`: it then cannot
  change a bit of the fragment likelihood and its DP is never run, which is
  two of the three passes an overlapping pair used to need. The all-mapped
  hypotheses are always scored first, in accumulation order, so the decision
  is made on the same running value the old loop had. The first version of
  the bound (`score ≤ -normalization`) never fired: it ignored the 0.25
  emissions and was ~100 nats too loose.

Around the kernel, the per-step work that gated wall time on the hot segment
is parallel or indexed now: `phase::score` plans fragments with `par_iter`
(each plan is independent; they run as one batch afterwards), `Linkage`
records which fragments observe each site so a block selects its fragments by
index instead of rescanning all ~6,600, `estimate_with_pairs` lists a CpG's
carrier haplotypes once (summing in haplotype order, as before), the
profile rows take one `powf` per base and derive the mismapped-mate rows from
them, and `cpgs`/`template` look alleles up in a window-filtered haplotype.

`src/call/haplotypes/pair_hmm/tests.rs` holds the oracles: the old closure recurrence
(`forward_scaled_reference`) and the old mixed emission
(`reference_score_with_model`), compared to the new path with `to_bits()`,
plus `bounded_mapping_hypotheses_are_skipped_without_changing_the_likelihood`
(a clean pair skips both mismapped-mate hypotheses, a pair whose mate carries
40 mismatches must score that mate's hypothesis). Break the kernel and these
fail.

**Measuring.** `tmp/phasing-perf/` (ignored) has the harness: `run_region.sh`
(`call --phase` with the centromere and ENCODE blacklist masks, 16 threads,
`/usr/bin/time`), `compare_vcf.sh` (byte identity of every record and every
header line but the command line) and `bench.sh`; baseline VCFs are in
`ref/`. Submit through SLURM with absolute paths and `--chdir`: `sbatch`
does not inherit the shell's `cd`, and compute nodes do not see the login
node's `/tmp` (`tabix` for the BED-reader unit tests is staged in
`tmp/phasing-perf/tools`). Profile with `perf record -e cycles:u` (user-space
sampling works with `perf_event_paranoid=2`); release binaries keep symbol
names, and the lane kernels show up as `run_avx512<8>` etc.

**Results** (`--threads 16`, baseline = commit `06691dd`, both arms on the same
SLURM node except the 30 Mb row, every output byte-identical to the baseline's):

| region | baseline wall / CPU | now wall / CPU |
| --- | --- | --- |
| chr1:80000001-81000000 (257 hets in one segment) | 33.7 s / 53.7 s | 1.7 s / 9.7 s |
| chr1:10000001-11000000 | 7.9 s / 16.0 s | 0.6 s / 4.3 s |
| chr1:121000001-126000000 (5 Mb, masked centromere) | 40.9 s / 61.5 s | 3.2 s / 23.6 s |
| chr1:1-30000000 (1,055,652 records, 11,942 phased) | 500 s / 1,858 s | 38 s / 193 s |

The unphased call is 0.24 s / 1.2 s on the first region and 1.7 s / 14 s on the
5 Mb one, so phasing itself went from ~33 s to under a second of wall time
there. CPU time on the compute nodes inflates relative to wall once the nested
batches share SMT siblings; on the 64-core login node the first region takes
1.1 s / 6.1 s. The profile after the change is ~35 % in `run_avx512<8>` and its
partial-lane variants (the row scale's divisions are most of it; a division
cannot become a reciprocal multiply without changing bits), ~10 % template and
weight construction, and 13 % in the variant-calling forest, which phasing
does not touch. The next exact wins, if ever needed, are mixed-shape lanes for
the tail batches and precomputed column weights per (window, orientation).

## ML feature layout (`src/metrics/ml/features/`)

### Experimental SNP haplotype likelihood prototype

`src/call/haplotypes.rs` exposes `infer_snp_expanded`; `likelihood.rs` owns conversion, sequencing and
mapping emissions, quadrature and phase-normalized genotype priors. `diagnostics.rs` drives
`rastair ml haplotypes` (seqair only) using full records, independently of called genotypes
and the main pileup's early filters. All focal A/C/G/T hypotheses survive. The prototype
defaults to CIGAR-conditioned scoring and shared per-CpG shrinkage with bounded fragment-connected
haplotype expansion. The integrated model and its CpG-driven context pruning are removed. It does
not modify the bundled forests or production calls. Default conversion is 0.95, background
conversion 0.003. Do not fit conversion efficiency directly to `bacteriophage_lambda_CpG`:
the spike-in may be incompletely converted and underestimate it by 2–5 percentage points.

Run `cargo test --features experimental-seqair --test haplotypes_cli` for BAM-backed tests.
`--sites-bed` takes one-base, 0-based half-open BED intervals without alleles and includes
zero-coverage targets. This is important when comparing TAPS and WGS: their naturally
discovered candidate sets differ, so calibration must use common focal positions.
`scripts/validate_haplotype_likelihoods.py --inputs <manifest.json> --rastair <seqair-binary>
--tabix <tabix> --output-dir <dir>` performs that comparison and a quadrature sensitivity
check, restricting scoring to the truth set's confident regions. `--pair-hmm` adds matched
realignment arms; `--include-indel-neighborhoods` includes and separately scores nearby indels.
The local NA12878 input manifest is `tmp/haplotypes-validation/inputs.json`; results are in
`tmp/haplotypes-validation/chr12/report.json`. Reads and generated diagnostics stay out of git.
The matching full reference, including spike-ins, is
`/mnt/lustre/shared/genomes/hg38+HPV+HBV+HVC+EBV+TAPS_SI_v2/hg38_full_gatk_HPV_HBV_HCV_spike-ins_v2.fa`.

`pair_hmm.rs` implements the opt-in `--likelihood-engine pair-hmm` local forward model and
`infer_snp_expanded`. CIGAR insertion/deletion proposals seed one nearest biallelic nuisance
indel (maximum 20 bp by default); pair posteriors marginalize its genotype and phase. Candidate
sequence adjacency, including indel edits, determines CpGs. Read sequences are cropped to the
selected loci plus flanks; masked interior bases stay in the sequence as neutral observations.
Shared rates are estimated separately for every candidate CpG, including realignment flanks.
Overlapping concordant CIGARs form a molecular profile with one consensus base and one
conversion event. Discordant CIGAR overlaps/refskips use aligned SNP likelihoods with the candidate's
edited CpG context and report `aligned_fallback_fragments`; indel alignment evidence is unscored.
Latent C/G discovery also checks neighboring SNP domains and indel-created CpG partners, so
fully converted alleles need not show an unconverted base to seed a focal hypothesis. Uniform insertion
emissions and path-mass normalization keep missing qualities uninformative across haplotype
lengths. Tests in `pair_hmm/tests.rs` compare the forward sum with exhaustive paths and sum the
observation probabilities to one; BAM tests show insertion realignment removing a false SNP.
Do not introduce a mapping event for a mate outside the scoring window: fetching that mate
depends on the input chunk and otherwise changes floating-point results. PairHMM reference
windows are focal-site-relative so normalization and CpG annotations also survive chunking.
Inactive mates must also be excluded from orientation consensus and reference-skip fallback;
the regression test combines a conflicting orientation with a skip on the unused mate.
The matched second-stage smoke report is
`tmp/haplotypes-validation/pair-hmm-final-chr12/report.json`; it includes indel neighborhoods
and PairHMM quadrature sensitivity. Slice results are not whole-genome accuracy estimates.

The earlier v1 expanded smoke suite is in `tmp/haplotypes-validation/multichrom-smoke/`: `suite.json`
records three 3-kb windows on chr1, chr12 and chr19, selected before caller results; `report.json`
combines matched aligned/PairHMM arms, probability checks and read support. Of 9,000 bp, 8,453
are GIAB-confident, containing 17 truth SNP alleles and five truth indel records. At SNP allele
presence posterior >= 0.5, both engines recover all 17 alleles; WGS has one extra allele and
TAPS has two. Each extra has only one informative supporting fragment, despite high base
quality: chr19:10002962 A>T (WGS), chr1:10018946 C>T and chr1:10019303 A>G (TAPS).
The C>T site also has two chemically confounded OT T fragments. The ordinary WGS emission
model applied to TAPS produces 90 extra alleles. PairHMM's 16/32-point quadrature comparison
changes SNP presence posteriors by at most 2e-18. These windows contain no truth SNP allele
within 10 bp of a truth indel, so they do not test near-indel SNP sensitivity. The extra calls
motivate calibration before production integration; do not infer genome-wide precision from
this smoke suite.

`scripts/benchmark_haplotype_likelihoods.py` performs the broader assessment. The earlier v1
assessment's local
manifest and results live in `tmp/haplotypes-validation/assessment/`. Three fixed 100-kb
windows contain 282,271 confident bases and 433 truth SNP alleles, including five near indels.
Discovery runs over each complete interval; inference uses all truth SNPs, all near-indel
negative candidates, known disagreements and 400 SHA-256-sampled negative candidates per
chromosome. This is 1,847 shared focal positions. Accuracy counts describe that panel;
inverse sampling weights estimate calibration on the discovery union plus forced positions,
not on all genomic bases. chr19 is recorded as held out for subsequent tuning, but these
results and the earlier chr19 regression are now inspected; they are not an untouched future
validation set. No likelihood parameter was fitted during this assessment.
The ordinary BED-reader unit tests also invoke an external `tabix`, outside the
`external-tool-tests` feature. On this checkout the validation commands used
`PATH=/tmp/rastair-tools:$PATH`; without it, two tests fail with "Failed to run tabix".

PairHMM recovers 433/433 WGS and 429/433 TAPS SNP alleles, with 5 and 13 extra panel alleles
respectively. There are 2 WGS and 8 TAPS genotype errors among the 433 truth SNP sites. All five
near-indel SNPs are recovered in both libraries. PairHMM does not uniformly improve scores:
TAPS has the same aggregate TP/FP/FN as aligned scoring, with a worse weighted log loss.
The TAPS restriction removes 1,955 neighbor-site incidences and 12,449 fragment-window links
at 1,144 panel sites. Do not confuse `omitted_sites` (all omissions, including the three-site
cap) with `cpg_omitted_sites` (only CpG-driven omissions after the cap), or these link counts
with unique molecules. They count potential coverage links, including realignment spans,
and do not guarantee informative observations at both sites. The diagnostic's genotype arrays use the ten A/C/G/T genotypes in
VCF triangular order, not REF/ALT indices. The Python scorer must reject multi-base REF
records even when a homozygous deletion's ALT is a single base; otherwise it labels indels
as SNP genotypes. Run its scoring tests with `python3 -m unittest discover -s scripts/tests`.
Of the 433 truth SNP sites, TAPS PairHMM reduces context at 169, and 94 lose every offered
neighbor link. All four missed truth alleles have reduced context; this association does
not establish that an unrestricted model would recover them.

`case_analysis.json` reproduces the three original disagreements from BAM slices and records
read-level perturbations. Each informative alternate is one fragment with two Q40 observations.
Keeping only one alternate observation changes PairHMM posteriors from 0.514, 0.996 and 0.921
to 3.52e-5, 8.88e-3 and 4.07e-4 respectively. Mates already share conversion and chromosome
assignment; the amplification follows independent sequencing errors, without a correlated
pre-sequencing substitution component. OT/OB flags agree between mates. Guessing orientation
removes the C>T disagreement, but that fragment's sole orientation motif is the focal mismatch
itself, so this does not establish a strand-assignment error. The A>G fragment has four CA
motifs per read away from the focal mismatch, supporting its existing OB assignment.

The singleton scenarios run as unit and BAM-backed CLI regressions. The earlier v1 regression
was checked by removing the second sequencing observation and observing a failure. The current
regression checks that duplication cannot amplify evidence.
`independent_two_cpg_oracle_recovers_linkage_against_a_singleton_error` is a test-only oracle
with independent per-CpG rates and negligible mapping error. With nine molecules, five
Gauss-Legendre nodes integrate the rate polynomials exactly. In the pre-consensus model, retaining
the neighbor lowered a constructed false SNP posterior from >0.99 to <0.01. The oracle now uses
the same consensus as the bounded model; retaining the neighbor still lowers that posterior
by more than 100-fold. This motivates evaluating a multi-CpG
model; it does not constitute a production implementation or a genome-wide improvement.

### Overlap consensus in haplotype likelihoods

The v2 likelihood engines resolve overlapping mates using `likelihood::consensus`, shared by
aligned SNP emissions and PairHMM profiles. Agreement takes the greater base quality without
adding qualities. A base disagreement takes the better-quality call with Phred quality equal
to the quality difference; equal-quality disagreements are Q0, with uniform base emissions.
Unknown bases provide no competing observation. Mapping quality breaks base-quality ties.
This happens after read filtering and before conversion, so TAPS C/T and G/A disagreements
cannot be interpreted as independent chemical events. Disjoint read pieces retain their
evidence and share one chromosome assignment.

Every mapping hypothesis must score the same consensus profile, replacing bases owned by
unmapped alignments with uniform emissions. Returning to a sum of individual read profiles
under partial mapping hypotheses would reintroduce overlap counting. A read owning no
consensus bases has no mapping event. CIGAR geometry remains intact, including neutral
interior bases; conflicting overlapping paths still fall back to the aligned model, which
uses the same consensus. Shared insertion bases also get a single uniform insertion emission.

`SnpEvidence.focal_observations` and TSV `observations` now count molecular bases. The raw
filtered measurements live in `focal_read_observations` / `read_observations`; overlap and
base-disagreement counts are in `overlapping_fragments` and `discordant_fragments`. These
refer to the focal SNP, not every base in the realignment window. A read spanning a deletion
at the target can supply alignment evidence without an aligned base observation. Model IDs
are `aligned_snp_v2` and `local_pair_hmm_v2`; historical v1 diagnostics are not comparable
without accounting for the observation policy change.

The v2 matched validation is in `tmp/haplotypes-validation/overlap-consensus/`, using the same
300-kb intervals and byte-identical target panels as `assessment/`. Candidate discovery is
unchanged. `case_analysis.json` reproduces the original three extra alleles with posteriors
3.52e-5, 8.88e-3 and 4.07e-4, each from one molecular observation rather than two independent
Q40 emissions. Unit and BAM tests cover duplication invariance, equal and unequal conflicts,
both TAPS strands, disjoint mates, indel realignment, fallback and chunk invariance. Deliberately
adding agreeing qualities or removing the disagreement penalty makes the consensus regression
fail; both mutations were restored.

`comparison.json` and `assessment.txt` record the matched before/after results. On the 1,847-site
panel, WGS aligned extras decrease from 6 to 4 and PairHMM extras from 5 to 3; TAPS aligned
extras decrease from 13 to 8 and PairHMM extras from 13 to 9. Recall is unchanged (433/433 WGS,
429/433 TAPS). PairHMM genotype errors decrease from 2 to 1 in WGS and 8 to 6 in TAPS, while
all five truth SNPs near indels remain recovered. These panel counts are not whole-genome
accuracy estimates. The one-CpG context restriction remains unchanged: 169/433 TAPS truth
sites lose context and 94 lose every offered neighbor link. Verification: 651 default tests,
680 seqair tests, four Python scoring tests, both Clippy feature sets, documentation and
formatting checks pass; `cargo xtask insta` passes 651 tests with no snapshots to review.

### Shared CpG shrinkage and neutral conversion evidence

`ml haplotypes` uses shared CpG shrinkage exclusively; `--methylation-model` is removed.
Neutral emissions remain internal to the pilot. Rates use Beta(2,2), molecular CpG-carrier
probability >=0.9 and leave-one-fragment-out conversion observations. These are approximate
plug-in scores calibrated with CIGAR coordinates. WGS bypasses conversion and rate fitting.

`expansion.rs` grows focal-rooted contexts through molecular connections. Defaults are
initial radius 150 bp, expansion radius 1000 bp, 64 retained haplotypes, 32 holdouts and
32 SNP sites as a safety cap. Neighbors require two conversion-unconfounded alternate
molecules and links from two independent fragments, with consensus BQ/source MAPQ >=20.
Conversion-only CpGs cannot form genetic bridges, but their observations remain available
for methylation fitting and scoring. Candidate CpG partners are inferred jointly, including
latent C/G neighbors, with two qualified molecules covering either partner. This context
dependency prevents fixing a variant partner to reference and falsely creating CpGs at
uncovered targets; it does not establish a distant genetic connection. Haplotypes rank by diploid posterior occupancy;
reference, every focal genotype's best pair and each active allele's carrier are protected.
Holdouts are reconsidered at the next extension. Evidence and rates are recomputed, never
multiplied across overlapping stages. Genotype priors are shared among surviving phase
assignments. Pruning mass is conditional on the current proposals, not a global error bound.
TSV v3 reports disconnected/resource-limited sites, pruning/rescue counts and context bounds.
These contexts do not emit global phase blocks or change production ML/calling/phasing.

PairHMM's `Prepared::interval()` is the authoritative scoring and calibration window,
including nuisance indels. Do not reconstruct its bounds from SNP loci: distant indels can
extend the actual window. CpG annotation includes one partner base outside either edge,
but observations outside the scored interval cannot fit rates or contribute mate evidence.

The matched assessment is `tmp/haplotypes-validation/cpg-models-large-final/`: six predefined
1-Mb windows on chr2–7, 5,172,807 confident bases, 20,976 targets and 6,161 truth SNP alleles.
At allele-presence posterior >=0.5, integrated/neutral/shrinkage recover 6,026/5,965/6,021
truth alleles, with 50/45/46 extra panel alleles and 178/240/183 genotype errors. Shrinkage
recovers 58 truth alleles lost by neutral and loses two, a net +0.909 percentage points of
recall; neutral has lower weighted Brier error on this panel. Both alternatives remove
CpG-driven context loss. These panel results are not whole-genome precision estimates.
The matched 500-kb PairHMM confirmation is `cpg-models-pair-hmm-final/`: integrated and
shrinkage recover 566/571 truth alleles with no extra panel alleles; neutral recovers 561
with one extra. `cpg-models-sensitivity-final/` checks prior strengths 2/4/8 and 16/32-node
quadrature on the same large-panel targets. Verification passes 657 default and 687 seqair
tests, five Python tests, both Clippy feature sets, and snapshot review with no updates.

The expansion assessment is `tmp/haplotypes-validation/expansion-smoke/`, with the same
1,847 targets and 433 truth SNP alleles as the previous 300-kb assessment on chr1/12/19.
Compared with the archived three-site shared-shrinkage model, WGS keeps 433 TP/4 FP/0 FN
and genotype errors change 2 to 3. TAPS keeps 429 TP/4 FN, extra panel alleles decrease
6 to 5, and genotype errors stay 7. WGS/TAPS have 528/620 contexts exceeding three sites,
maximum 21/32 sites, 343/416 contexts with pruning and 93/44 with rescued candidates.
These are panel results, not whole-genome precision. PairHMM's scalar forward model remains
expensive for dense wide contexts; the completed chr1/chr19 smoke in `expansion-pair-hmm/`
uses explicit eight-site/32-haplotype/16-holdout budgets and the same growth radii, recovering
all six truth alleles in both libraries with no extra panel alleles. The dense chr12 run
was stopped for runtime and remains explicitly incomplete. Superseded
wide/default-budget runs are marked incomplete and must not be used as completed assessments.
Verification passes 666 default and 698 seqair tests, five Python tests and both Clippy sets.
Snapshot review passes with no updates. Removing CpG context-support edges makes the
uncovered-focal regression fail; the exact source was restored and all 67 likelihood unit
tests pass afterward.

### Existing feature layouts

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

## Molecular ML feature experiment (2026-10-14)

`ml train --ml-features haplotypes` appended, at the time, 13 features to all five frozen Standard layouts (68/69/67/30/30); since the merge of 2026-10-16 (see "Merged base" below) it appends 16 features to the three SNV layouts (71/72/70) and indels have no forest. Model metadata automatically enables full read capture and pre-ML PairHMM scoring in both backends; final `--phase` calls are never training inputs. `canonical` is the original-layout ablation with the same canonical indel preprocessing. The chr19 training and full chr12 validation in `tmp/haplotype-ml-chr19/slurm-20261014-optimized-64/` completed successfully. The bundled model remains unchanged pending investigation of the canonical indel genotype regression; improved F1 alone does not satisfy that validation.

`src/call/haplotypes/features.rs` owns shared likelihood parameters and preprocessing for training and calling. Simple indels are normalized before filtering/label lookup; reference tracts, per-molecule depth and query boundary metrics use that anchor. Reads retain their original CIGAR profiles. A boundary outside the preceding matching run or read mask is unresolved. Keep already-accepted indel allele observations even when the canonical boundary is lower-quality/missing: encode that uncertainty rather than imposing a second, different BQ cutoff. A supported SNP within an observed indel shift tract makes the complex chromosome hypothesis unsafe without assembly, so its molecular features are unavailable. Local ML alignment geometry includes the full reference-equivalent placement span; the diagnostic retains its existing geometry. Indels over 100 bp and wholly aligned-fallback indels have unavailable molecular features.

Training and calling use the same pre-ML eligibility predicate, de-novo adjacency setup and shared CpG shrinkage. Use normalized, split GIAB truth and remove unused ALT rows after splitting (`GT="alt"`). A deletion's full reference span must be inside the confidence BED. New molecular feature extraction fails the training run on any failed segment. Relative-weight normalization in `posterior_pairs` prevents large common negative log offsets from producing aggregate probabilities above one.

`tests/call_cli.rs::model_metadata_activates_pair_hmm_features_without_a_phase_flag` uses a tiny forest that responds to the availability feature; masking makes that feature unavailable. It runs on both readers and proves the model selected the actual capture/scoring path. `features.rs` and `reads.rs` tests cover canonical boundaries/qualities, shifted CIGAR placements, read-group identity, molecular deduplication, aligned fallback and CpG conversion guards.

`scripts/benchmark_haplotype_ml.py` provides separate train/call/compare stages. Run long experiments via SLURM `short` (24 h maximum), supplying realistic `-t`, `-c`, and `--mem`. Current jobs request 64 CPUs/96 GB/6 h for chr19 training, 64 CPUs/96 GB/12 h for full chr12 calling, and 8 CPUs/16 GB/1 h for comparison; record actual resources to refine these estimates. Keep binaries, source patches, masks, commands, sample names and complete-status markers in the run directory. Do not install a model from incomplete data or a failed job.

Compute nodes do not share the login/session `/tmp`: put test helper executables in shared storage and add that path to the SLURM job's PATH. GIAB's sample name is `HG001`, not the Platinum Genomes `NA12878`; read it from the VCF header. The augmented reference has IUPAC `Y` bases that Rastair renders as `N`; these can fail bcftools normalization outside GIAB confidence. Subset query calls to the confidence BED before normalization for benchmarking, retaining strict reference validation inside that BED. Do not use a permissive REF-repair flag to hide a coordinate error.

### Comparing GATK physical phasing

Check `PGT`/`PID` as well as `GT`/`PS` in HaplotypeCaller/GenotypeGVCFs output. GenotypeGVCFs can reorder the called genotype while retaining the physical phase in `PGT`: the supplied GATK 4.5.0.0 chr12 callset has 1,049 PASS SNP heterozygotes with `GT=0|1` but `PGT=1|0`. Scoring its final GT orientation alone therefore measures a representation problem as phase errors. [GATK's explanation](https://gatk.broadinstitute.org/hc/en-us/community/posts/28270342277147-Remove-the-ALT-Allele-Reporting-From-a-g-vcf-Before-or-After-Genotyping?sort_by=votes).

Use PGT/PID for the physical-phase comparison only when it unambiguously describes the final called alleles: the assessment restricts that interpretation to biallelic SNPs with final GT 0/1 and PGT 0|1 or 1|0. PGT is not a replacement for allele dosage or a multiallelic genotype. Preserve the original file and score GT/PS separately as a diagnostic. Variant F1 continues to use the original genotype, with phasing ignored by Aardvark. `gatk-comparison/` under the completed chr19/chr12 run preserves the original-GT comparison and the PGT/PID assessment, with WhatsHap crosschecks. Compare SNP-only blocks on identical confidence/mask domains and report common-link errors alongside phase coverage; neither PS labels nor a whole-block orientation flip is an error.

### Indel genotype regression audit (2026-10-15)

The following describes the pre-fix behavior preserved in the audit artifacts. Indel GT was assigned by `indel_calling::binomial_genotype` from clean molecule counts and fitted repeat-dependent allele fractions before ML runs. The forest predicts variant presence, not dosage. With ML enabled, indel VCF QUAL comes from that variant-presence score; indels emit no genotype GQ/GL/GC. A high QUAL cannot be interpreted as high confidence in 1/1. The count model always emits its winning genotype without a minimum confidence requirement.

`normalize_indels` previously computed reference noise from **all** reads in a fragment. A noisy non-covering mate can exclude a clean REF observation while the corresponding unrelated mate does not exclude ALT evidence. This is a measured regression, not a proposed explanation: in the completed chr12 GIAB-confident audit, 393 previously correct simple heterozygous indels become homozygous ALT, all with denominator loss and no increased target ALT count, changed anchor or changed tract. Restoring only reference observations excluded by remote mate noise reverses 367 of those decisions. The canonical and full haplotype models have identical genotypes/counts at all 30,700 shared simple indels. Correct this scope and validate genotype counts before retraining or changing the bundled default.

The other 26 audited changes involve local handling of rejected indel reads (15 sites, confirmed in BAM with >5 non-TAPS mismatches) or clean competing-allele observations lost during mate aggregation (13 sites; two overlap). Do not automatically turn rejected indels into clean reference evidence. The compound-het resolver also skips a locus whose first allele is already hom-alt, so biased single-allele decisions can obstruct joint interpretation. `indel-genotype-investigation/` under the completed experiment preserves before/after counts, unchanged-output hashes, BAM evidence and an isolated failing invariance reproducer for both backend builds. Custom tracing targets require `RASTAIR_LOG=target=info`; `RUST_LOG` is not Rastair's logging variable.

Pre-ML molecular scoring already computes unrestricted diploid genotype evidence: `SnpEvidence.genotype_posteriors` contains all ten unordered A/C/G/T genotypes, marginalized over retained haplotype pairs; `features::indel_likelihood` computes zero/one/two copies of one focal indel. Previously feature extraction retained only per-allele summaries, and `process_collected_pileups` cleared `ColumnPhase::Likelihood` before count-based genotype assignment. Reusing molecular genotype evidence required retaining the full distribution through this boundary. The final production phaser instead conditions on called genotypes and cannot supply an independent hom/het decision. Single-indel presence/copy marginals also do not provide a joint distribution over competing indel alleles or complex SNP/indel genotypes.

### Molecular genotype assignment (2026-10-15)

Models with `Haplotypes` metadata now retain the unrestricted pre-ML PairHMM genotype posterior through filtering. `haplotypes::genotype` assigns SNP GT before methylation and final phasing, marginalizing all ten DNA genotypes over retained chromosome pairs. Indels use the same pre-ML likelihoods for zero/one/two copies; competing alleles at one canonical anchor are scored jointly (up to seven alternatives), including latent anchor SNP genotypes, so a `1/2` is an actual joint hypothesis. Per-allele ML features retain their existing layout and scoring to allow evaluation with existing forests. All molecule likelihoods are used once; ML gates alleles without multiplying its score into the genotype posterior or renormalizing away excluded allele mass.

`ColumnPhase::Genotyped` preserves genotype provenance through final phase assignment without enlarging `PileupMetrics`. FORMAT `GTM` reports `pairhmm` or `counts`, and `GQ` reports model-derived genotype confidence capped at 99. Molecular GQ uses retained-state posterior error, limited by maximum observed pruning mass: it is conditional, not calibrated or a bound on all omitted hypotheses. Legacy scalar GL/GC are omitted on molecular genotypes; ML-based indel QUAL remains variant-presence confidence, including compound records. `ml-rescue` cannot override a molecular hom-ref genotype to force a het. Missing/unsupported molecular evidence, masked domains, unsafe complex repeat alleles, >100 bp events and over-budget competing-allele sets retain the count fallback. Models with Standard/Canonical metadata remain count-genotyped.

Canonical REF noise now considers only mates with usable observations at the actual column; non-covering noisy mates cannot remove clean reference support. Raw and normalized indel anchors survive read-level rejection, so a rejected indel cannot silently become clean reference evidence at either placement. Competing indel identities must appear in priors, PairHMM reuse keys and protected search states: a boolean present/absent key aliases distinct sequences and invalidates joint genotyping.

`tmp/pairhmm-genotyping/` preserves test logs, the frozen caller, commands, model hashes and chr12 validation outputs. The targeted panel combines 40 previously wrong homozygous calls with three 500 kb intervals; its error-enriched aggregate is not an unbiased chromosome estimate. Full chr12 evaluation compares the original forests without retraining, including a fixed canonical preprocessing arm. The bundled default remains unchanged until full genotype-aware validation has been reviewed.

The cargo-managed `mdbook-mermaid` shim currently attempts installation with `--version 0.16`, which cargo-binstall rejects. `cargo xtask docs` still generates the CLI and field references; verify the final book with the installed mdbook binary and `/users/bschuster/.cargo/bin` first on PATH to use the working Mermaid 0.16.0 binary. Restore host-dependent generated thread defaults to avoid unrelated documentation churn.


Full chr12 genotype validation finished successfully (jobs 194776/194779). The caller took 7h23m05s versus 7h22m27s previously; the finalized, indexed BCF has 2,795,772 records and passed a complete read-through plus indexed tail query. Against GIAB, changing from the same molecular forest with count genotypes to PairHMM genotypes moves SNP/insertion/deletion F1 from 0.978441/0.837510/0.850791 to 0.979569/0.812642/0.869645. **Do not promote this model yet:** insertion recall falls from 0.765208 to 0.700315 despite improved precision. Overcalled heterozygous insertion copies (`query_fp_gt`, EA=1/OA=2) fall 487→90, but undercalled homozygous insertions (`truth_fn_gt`, EA=2/OA=1) rise 477→1036. Deletion error directions similarly shift 897→159 versus 312→662. Reporting only query FP_GT as "genotype errors" conceals the opposite error direction. `full/allele-dosage-decisions.json` preserves both truth and query allele-copy decisions. Raw switch counts improve 16/31,554→13/31,648 against Platinum Genomes; on the identical 30,176-link intersection both make 10 switches, so this is not established improved common-link phase accuracy. Reports and original forests are preserved in `tmp/pairhmm-genotyping/full/`; the bundled default remains unchanged.

Live progress caveat: `runtime::segments::process_in_order` uses ordair's bounded look-ahead (four segments per worker). A slow early segment can therefore hold back both publication and admission of more regions after that window fills; idle workers do not establish that the rest of the chromosome is complete. `ProgressTracker` counts written segments, and its linear ETA can be misleading during this stall. In chr12 genotype-validation job 194776, at ~2h31m the output frontier was 8,293,365 with 246,027 records, while a live sample showed one busy core. A 4,258-sample profile confirmed ongoing AVX-512 PairHMM/template/prior computation with no lost samples. The profiler step was stopped separately; the caller continued. Preserve `tmp/pairhmm-genotyping/status-194776.perf` for follow-up.

**Important:** Whenever you learned something new about how to develop features, find code, or how to debug issues, you **must** add it to this document.
This is the single source of truth for how to work on this codebase, and it must be kept up-to-date with any new insights or changes.
If you find yourself asking "How do I do X?" and you figure it out, add a section here so that the next person doesn't have to ask the same question.

### compair is the PairHMM kernel (2026-10-15)

`pair_hmm.rs` scores every (profile, template) forward pass with seqair's `compair` crate; the
lane-batched `f64` kernel (`pair_hmm/forward.rs`) is gone. compair runs GATK's recurrence
(free start `D[0][j] = 1/h`, total over the last row's M+I) in banded `f32` SIMD kernels with an
`f64` rescue, and returns log10 (`align` converts to ln).

**Vendored, patched compair.** compair 0.2.0 emits an inserted base at one, as GATK does; rastair
emits it at 0.25, which keeps the observation model a distribution (at one, each inserted base
is 4x cheaper, ~10^6 for 10 bp, which moves indel genotypes). `vendor/compair` is the published
0.2.0 plus `Read::with_insertion_emissions`, wired in via `[patch.crates-io]` and excluded from
the workspace. The emission folds into the two transitions into the insertion state
(`Transition::insertion_continuation` beside `gap_continuation`, a ninth row track; in the pairs
kernel it rides beside the 8x8 `RowEntry` transpose), so emissions of one are bit-identical to
upstream. The GPU path refuses such reads (`GpuError::InsertionEmissions`). The upstream patch,
with tests, is `tmp/compair-patch/compair-insertion-emissions.patch`. The patched clone with the
full compair suite is `tmp/compair-patch/seqair` (all compair tests, the `gpu` feature and clippy
pass there). Drop the vendor copy once a compair release carries the change.

**How rastair's model maps onto compair.**
- A profile becomes one `compair::Read` per mapping hypothesis (`ProfileTerms`). Bases owned by a
  mismapped read get Q0 (eps capped at 0.75, a quarter under every column). Unknown bases get
  Q254 (emission one) and insertion emission one. Other bases insert at `INSERTION_EMISSION`.
  Gap open/extend must be whole Phred qualities (`phred` errors otherwise; defaults are Q45/Q10).
- A template becomes a `TemplateTrack` per read orientation: a `compair::Haplotype` whose betas
  carry each column's conversion rate. A stateless `TemplateEmission` maps a level onto the C/T
  or G/A pair, so neutral columns (0.5 at a template T or A too), shrinkage rates and window-edge
  CpG roles are exact, and compair's own CpG roles and the read strand are unused.
- **Accepted difference:** an unknown (`N`, incl. IUPAC) template base emits `1 - eps` (GATK's "N
  matches"), not the old 0.25; `SiteWeights` cannot express a uniform column, and all compared
  haplotypes share the column.
- Normalization (path mass) is computed with compair too, as unknown read vs unknown template
  through the *same band*. A band cannot unnormalize the model.

**Bands.** Each job's band is anchored at the template index of the profile's first base
(`ProfileSpan.first`; `partition_point` over template positions). Its half-width is
max(23, profile CIGAR indel bases + max(template indel bases, `Prepared::indel_reach`) + 16),
rounded to compair's `16n - 2` widths. **`indel_reach` is the trap**: a read can carry an indel
its CIGAR aligns straight through, and against the reference template neither side shows it.
The first band rule missed exactly that, scoring such a read 29 nats too low against the
reference. `Prepared::new` seeds the reach from its indel, and `likelihood_matrix_geometry`
widens it with `.reaching(...)` over every haplotype's indel (competing indel alleles).
`the_band_holds_a_path_across_a_long_indel` (60 bp deletion, CIGAR-carrying and straight reads,
both templates, against the full `f64` oracle) fails with a default-width band.

**Scorer.** The API (`new` is now fallible, `plan`, `run`) and the memo/bound/skip logic are
unchanged. Memo keys are interned read bytes (base, quality) and template (base, level bits), plus
band offset and width; normalizers are keyed by `Shape` (rows, columns, offset, width). `align`
sorts pairs by (band width, read, template), because a compair pairs-kernel group shares one
width and caches derived tracks by address. It scores chunks of 64 pairs on a thread-local
`Workspace`, on the enclosing rayon pool, and sends chunks shorter than `PAIRS_BREAK_EVEN`
through `align_candidates` (the strip kernel). `ProfileTerms::bound` stays valid: a band only
removes paths.

**Precision and tests.** Scores agree with the `f64` recurrence to ~1e-6 relative.
`pair_hmm/tests.rs` compares against the old closure oracle with `COMPAIR_PRECISION`, not
`to_bits()` (`compair_reproduces_the_mixed_emission_through_the_full_matrix`; its oracle treats
an unknown template base as compair does). Outputs are no longer byte-identical to pre-compair
runs, so `tmp/phasing-perf/ref/` VCFs are pre-compair references.

**Validation** (`tmp/compair-integration/`, ignored; `report.json`, `panel/phase-accuracy.json`).
The baseline is rebuilt from HEAD plus the staged diff, both arms in one SLURM job on node01 at
64 threads, with the haplotype-metadata model on the `tmp/pairhmm-genotyping/panel` intervals:

| | baseline f64 | compair |
| --- | --- | --- |
| panel calling wall / CPU | 908 s / 15,407 s | 671 s / 6,273 s |
| chr1:80000001-81000000 `--phase`, bundled model | 1.79 s / 25.4 s | 1.28 s / 17.8 s |

Aardvark GT metrics are identical in every column for SNV, insertion and deletion (recall,
precision, F1, `truth_fn_gt`, `query_fp_gt`), and match the frozen `pairhmm-genotypes` arm. All
35,196 records agree on GT, FILTER, QUAL, phase and PS; 12 differ in GQ by 1-2. Platinum phase:
514 phased SNPs, 0 switch errors in both. On the chr1 region every non-phase field is identical,
with 524 phased records in both. CPU drops 2.5x but panel wall only 1.35x: wall is now bound by
something serial (~9 of 64 cores busy on average), likely the slow-segment look-ahead stall
described under "Live progress caveat". Profile that next, not the kernel.

**Building a baseline tree under `tmp/`:** `git apply` run inside the repository's `tmp/`
silently applies nothing, because it ignores patch paths outside the current directory. Use
`git archive HEAD | tar -x -C <dir>` and then `patch -p1 < staged.patch`, and diff a few files
against `git show :<path>` to confirm.

### SNP-dense regions: pre-ML scoring cost (2026-10-15)

With compair in place, the haplotype-metadata model (`--experimental-indels=ml --phase`, 64
threads) spent ~99 % of its CPU in `add_likelihoods` → `infer_snp_expanded`. Dense regions such
as chr12:1-12M hold columns with 24-32 context sites (31 expansion steps, 64 haplotypes, 128
proposals, 8,256 diploid pairs), which took 3-3.5 s each. Most of that cost was not the
alignment kernel. The changes below are all exact, and every run's VCF was identical to the
previous version's:

- **One `Scorer` per task, not one per matrix pass.** `likelihood_matrix_geometry` used to build
  a fresh `Scorer` for every pass, so the content-keyed memo never hit (0 hits measured), and
  every expansion step realigned every (profile, template) pair. The `Scorer` is now threaded
  through `infer_snp_expanded`/`expand_indels`/`score_indel`/`joint_indel_genotypes`. Memo keys
  are exact content, so sharing a `Scorer` across passes and columns cannot change a result.
- **Columns run in parallel.** `add_likelihoods` scores eligible columns with
  `par_chunks(COLUMNS_PER_TASK).with_max_len(1)`, with one `Scorer` per chunk, and writes the
  results back in column order. **`with_max_len(1)` matters:** without it rayon's adaptive
  splitting handed one thread ~22 consecutive heavy columns, which ran serially at the end of
  the run. On chr12:1-12M that tail cost 250 s of a 870 s wall.
- Template windows are index ranges into one base vector. `TemplateTrack` stores only its memo
  key, and the `compair::Haplotype` is decoded only for memo misses.
- `retained_pair_priors` precomputes allele and indel codes and per-site log priors, and keys
  haplotypes with an `FxHashMap` (`retained_pair_priors_reference` is the bit-identical test
  oracle). `methylation::estimate` reuses the per-(pair, fragment) mixtures it already computed.
  `log_add` has an equal-argument fast path.

**Results** (`tmp/compair-profiling/`, ignored; one SLURM job per arm on `short`, 64 threads;
"compair" is commit `8d3e99b`):

| region | compair wall / CPU | now wall / CPU |
| --- | --- | --- |
| chr12:8.3-8.4M (densest segment, alone) | 4,154 s / — | 65.5 s / 2,692 s |
| chr12:1-12M (321,198 records) | 4,430 s / 87,421 s | 620 s / 37,232 s |
| chr19:1-6M (384,293 records, 3,009 phased) | >2,051 s, cancelled, 58/60 segments | 436 s / 24,886 s |
| chr12:60-72M control (195,165 records, 4,786 phased) | not run | 226 s / 13,071 s |

The runs now average 57-60 of 64 cores, so wall time follows CPU. The cancelled compair chr19 run had
written output only up to chr19:297,515: the remaining two segments, chr19:299802-400201 and
599802-700201, blocked publication through the look-ahead stall (see "Live progress caveat").

**What is left is structural.** On the dense segment the profile is flat: compair kernel ~13 %,
step plans ~11 %, pair priors and `posterior_pairs` ~11 %, `exp` ~6 %, memo bookkeeping ~7 %.
The cost is O(steps × haplotypes² × fragments) per candidate column, repeated for every
candidate, and ~97 % of SNP-scored columns there end up non-variant. Further gains need
model-changing work, so they need the full genotype validation and the user's approval:
region-based joint inference (as in GATK and Octopus), a cheaper candidate pre-gate, or
dense-region haplotype caps.

**Profiling harness** (`tmp/compair-profiling/`):
- `instrument.py <tree>` patches a *copy* of the source with span-close timing and stage spans.
  In parallel-column mode it adds per-column `COLUMN snp` logs, which carry `pos`, `ms`, sites,
  steps and haplotypes.
- `run.sh <region> <name> [threads]` (the binary comes from `BIN`) runs the call with a CPU
  sampler.
- `analyze.py runs/<name>` reports per-segment stage times and the busy-core timeline;
  `columns.py`/`tail_columns.py` analyse the per-column logs.

Rayon worker threads do not inherit the segment span, so attribute a column by `pos` rather
than by span. Under parallel columns, segment busy times overlap and stop being a cost measure;
use the per-column logs and `perf` instead. The CPU sampler must follow the `rastair` process
(`pgrep -x rastair`), not `/usr/bin/time`.

### Investigating molecular indel undercalls (2026-10-15)

`tmp/indel-dosage-mechanism/findings.json` preserves a diagnostic investigation at `8a6cc42`,
including the actual optimized scorer, isolated instrumented source, scalar oracle and ten
5 kb real-data reruns. No production likelihood change was made. Distinguish three mechanisms:

- `Prepared::new` crops molecular sequences to the same reference-coordinate interval used
  for templates. A spanning insertion-bearing profile can therefore be longer than REF.
  `Scorer` subtracts a separate path-mass normalizer for each profile/template shape. When REF
  requires insertions just to accommodate that length, its normalizer is tiny, cancelling much
  of the gap penalty. Deletion-bearing profiles do not suffer that cancellation against REF.
  With 20 ALT molecules and one REF outlier, a mirrored 5 bp test through current compair scores
  insertion as 99.8015% het and deletion as 99.2166% hom-alt, with methylation disabled. This is
  a property of the cropped, length-conditioned model, not floating-point drift. Do not simply
  remove normalization: its neutral-observation invariants still matter. Separate evidence
  cropping from padded template geometry when testing a replacement.
- The genotype mixture gives each decisive ALT molecule at most ln(2) evidence for hom-alt
  over REF/ALT; a REF-favoring outlier can contribute much more in the opposite direction.
  PairHMM uses fixed Q45 gap-open / Q10 extension, without the old count model's repeat-dependent
  hom-alt fraction (down to 0.683) or an explicit local-artifact/stutter model.
- `normalize_indels` removes noisy REF from count evidence, but the full molecular likelihood
  does not use `ReadAnnotations.noisy` as an artifact probability. At chr12:10,905,014 T>TA,
  AD=0,9 accompanies GT=0/1, GQ=13: one terminal-repeat-marked molecule contributes log(ALT/REF)
  = -8.5209. Removing just that molecule from the saved, two-haplotype matrix changes the
  conditional-on-variant result from 95.245% het to 99.208% hom-alt. This diagnostic removal
  establishes influence; it is not a validated proposal to discard all noisy reads.

In the direct normalized-allele intersection of both chr12 query callsets, true hom-alt
insertions are undercalled at 1,022/4,141 sites, deletions at 647/4,199. Of these, 983 and 614
have at least six reference-equivalent right shifts; count shifts *after the deleted sequence*
for deletions, rather than counting the deletion itself as repeat length. Insertions with only
A/T inserted bases account for 806/1,022 errors. These subsets are not the full Aardvark totals
(1,036/662); complex/absent alleles are excluded. Trace unsplit calls too: normalization turns a
wrong 1/2 into one-copy records. Compound genotypes explain 29 insertion / 47 deletion errors
here, only 7/723 and 14/477 newly regressed calls, so they are a secondary contribution.

The mechanisms have controlled and real-locus evidence, but their separate whole-chromosome
accuracy contributions have not been measured. Fit any repeat/artifact parameters on chr19,
keep chr12 held out, and assess both dosage error directions. The test and tracing worktree
remains isolated under `tmp/indel-dosage-mechanism/tree`; production source is unchanged.

Upstream comparison: `tmp/indel-dosage-mechanism/upstream/` preserves source from GATK 4.5.0.0
(`8317d8b`) and Octopus `ed012a6`. GATK's `LoglessPairHMM` returns its forward sum with a uniform
haplotype-start prior and numerical scaling; it does not divide by Rastair's length-conditioned
path mass. Its insertion-state emission convention is one, whereas Rastair uses 0.25. Do not
confuse `AlleleLikelihoods.normalizeLikelihoods` with our path normalizer: it floors per-read
scores relative to the best haplotype (Q45 maximum separation by default in 4.5), symmetrically
including REF. GATK also filters poorly modeled reads and adjusts indel qualities with its PCR
model, or supplied DRAGstr context parameters. Its documentation assumes haplotypes long enough
to contain reads; assembly has padded regions and read clipping, not a claim of never clipping.

Octopus explicitly checks read length plus HMM padding and throws `ShortHaplotypeError` when
it cannot fit. The inspected `pair_hmm.hpp`/`simd_pair_hmm.hpp` scoring path is a banded best-path
penalty approximation, not GATK's forward sum over paths. There is no equivalent length-conditioned
path-mass division. Its repeat model assigns both gap-open and gap-extension penalties by the
candidate haplotype's repeat motif and length, with library/technology-specific or custom parameters.
With mapping quality enabled, its relative likelihood is mixed as `(1-e)*L + e`, assuming an
off-locus score of one. Rastair's off-locus sequence probability instead decays as `4^-n`; that
component can be negligible for a locally wrong read whose flanks match, so it is not equivalent
to Octopus's floor or GATK's likelihood-ratio cap. These are modeling differences, not a guarantee
that either upstream caller avoids all hom/het errors. The supplied filtered GATK VCF retains a
GenotypeGVCFs 4.5.0.0 command but no HaplotypeCaller command, so its exact likelihood options cannot
be established from that header.

### Regional forward likelihoods (2026-10-15)

Production molecular features now use `haplotypes::regional`, with one diploid posterior
per connected region (32 sites / 2 kb, 64 retained and 32 buffered chromosomes). Compatible
indels at multiple anchors can coexist. Pruning protects the best pair for every local SNP
and indel genotype plus the all-reference chromosome. Insufficient protection capacity or
a deletion spanning the next candidate splits the region before that candidate; never erase
the earlier deletion chromosome to accommodate an interior anchor. Site marginals supply
both forest features and molecular genotyping. The focal-rooted diagnostic and called-SNP
phase solver remain separate inference searches sharing the likelihood kernel.

This supersedes the length-conditioned scorer described above. Compair now uses insertion
emission one, with no separate path-mass normalizer. Full molecular profiles within each
unmasked reference domain have independently padded templates: 12 bp plus the maximum sum
of indel lengths on a candidate chromosome. A symmetric per-molecule cap, bounded by source
MAPQ and the library profile's maximum (default Q45), applies to both neutral and final
matrices. Repeat errors are evaluated on candidate DNA before conversion, and projected into
compair read transition qualities using original coordinates; this is not a per-path,
per-template-column transition model. Cache keys include the resulting gap tracks.

`haplotypes::ErrorProfile` stores platform/library names, separate insertion/deletion opening
qualities and extension qualities by repeat period/copy count. `ml train
--haplotype-error-profile FILE` embeds it in the forest; call and phase use that profile.
Older Haplotypes forests lacking this metadata are rejected at caller initialization, since
their feature definitions differ. Standard/Canonical forests still load. Default profiles
use an explicitly conservative GATK-style repeat curve. Calibration is available through
`scripts/calibrate_haplotype_errors.py`: confident hom-reference repeat tracts with matching
flanks, known-variant exclusion, molecule deduplication and shrinkage of sparse contexts.

Regression coverage includes mirrored insertion/deletion outlier cases, multiple indels
plus SNP inference, overlap deduplication, local genotype protection, spanning-deletion
splitting, raw-forward oracles, profile serialization and invalid-profile rejection.
Run manifests, calibration counts, frozen old/new binaries and SLURM validation artifacts
live in `tmp/regional-pairhmm/`. Chr19 training and held-out chr12 evaluation must complete
before any bundled forest promotion; submission alone is not validation.

### Regional deletion sensitivity audit (2026-10-16)

The completed 2 Mb chr19 / 2 Mb chr12 pilot lives in
`tmp/regional-pairhmm/pilot-2mb/`. Relative to the equally trained Canonical
control, deletion truth-side recovery loses 19 alleles and gains 3: 11 losses
are ML-passing PairHMM `0/0` calls, 2 are hom-to-het, and 6 are filtered out.
Use `call --all` to retain rejected indels for diagnosis (`--vcf-all-fields`
only expands fields). A PASS indel can still have GT=0/0; filter both fields
when making a benchmark query. Aardvark truth-side TP alone does not establish
correct query dosage: excess copies are penalized on its query side.

The regional model inherited `indel_prior=0.0001` from the nuisance-indel model.
At a biallelic locus its 0/1:0/0 prior odds are 1:19998. In the calibrated
repeat model, chr12:10270785 AT>A has AD=9,6 and ML=0.78947, but GT=0/0:
the marginalized read likelihood favors het over reference by about 85:1,
insufficient to overcome this prior. The resulting reference posterior is
0.99577. This example reproduces the original pilot's call in an isolated
trace. The same issue appears in training: 37/214 true deletions with available
molecular features have P(variant)<0.5; 33 of those have positive log Bayes
factors. More forest examples cannot directly change a later molecular veto.

`tmp/regional-pairhmm/deletion-mechanism/` preserves the isolated instrumented
source, matrices, profile ablations, training-feature audit and `summary.json`.
Its 7 kb windows are error-enriched; some change context relative to the
original 100 kb windows, so check exact call fields before claiming identity.
Replacing only the error profile with the conservative preset restores
non-reference calls at all 16 traced hom-ref losses, but does not establish
overall accuracy or resolve the two original hom-to-het errors. Forests must
be retrained when profiles change; diagnostic overrides are confined to the
isolated tree. Local prior reweighting conditions on retained states and fixed
methylation estimates, not a fresh complete search.

The follow-up `tmp/regional-pairhmm/pilot-4mb/` doubles each original interval
symmetrically, expanding adjoining original windows together to retain exactly
4 Mb per chromosome without overlap. It retains all original positions and
selects no regions based on errors. Calibrated and conservative profiles each
receive matched forest retraining and calling; `jobs.json` records dependencies
and the combined report job. Keep these pilot results separate from a full
chromosome validation and do not promote the bundled default from them.

The 4 Mb runs completed (194839–194850; combined report at 09:31:58 BST on
2026-10-16). Calibrated training/calling took 36m51s/27m28s; conservative
36m55s/28m11s. Genotype-aware SNP/insertion/deletion F1 is respectively
0.974176/0.830063/0.857295 for the equally trained Canonical control,
0.976374/0.831112/0.838387 for calibrated PairHMM, and
0.975779/0.827098/0.857426 for conservative PairHMM. Conservative improves
deletion recall (0.740047→0.779859) but lowers precision (0.966867→0.952128).
The calibrated arm loses 39 truth-side deletion recoveries relative to
Canonical and gains 10; 31 losses are PASS `0/0`, 3 hom-to-het and 5 absent.

The conservative profile reintroduces dosage undercalls even on identical
emitted alleles: among 142 shared true homozygous insertions, hom-to-het errors
are 9 calibrated versus 25 conservative; among 134 shared homozygous deletions,
7 versus 21. Shared heterozygote overcalls are insertion 7 versus 5 (200 truth
hets), deletion 2 versus 2 (189 truth hets). See `shared-allele-dosage.json`:
comparing only total errors without this intersection confounds dosage with
missing variants. This supports investigating the variant-presence decision
while retaining calibrated dosage evidence, not adopting the conservative
profile on deletion F1 alone. No genotype-policy change or default promotion
was made. Phase switches are 1/1205 for calibrated and 0/1198 conservative;
all four arms have 0/1144 switches on their shared links, so this pilot does
not establish a difference in shared-link phase accuracy.

### Review of the regional PairHMM: why it does not beat the count model (2026-10-16)

The evidence is the 4 Mb pilot (`tmp/regional-pairhmm/pilot-4mb/`, forests trained on 4 Mb of
chr19: 224 deletion / 201 insertion positives) and a per-record audit of it in
`tmp/regional-pairhmm/veto-audit/` (`audit.py`, `results.txt`; every confident truth variant
classified as correct / absent / filtered / PASS-but-0/0 / wrong dosage, per arm, with each arm's
own aardvark annotations). The full-chr19 training (job 194825) timed out at 6 h, so no
whole-chromosome result exists.

**Root cause of the indel losses: the molecular hom-ref veto.** `genotype::assign_snp` /
`assign_indels` take the argmax of the full diploid posterior with hom-ref always eligible, so an
allele the forest accepts is emitted as PASS `0/0` whenever the prior-weighted posterior favours
reference. On the pilot the regional arm does that at 22 SNV / 24 insertion / 38 deletion truth
variants against 4 / 7 / 8 non-variants. With the consistent per-record metric, lifting the veto
(GT from the non-reference genotypes once ML passes) bounds F1 at SNV 0.979 / ins 0.864 /
del 0.887 against canonical 0.976 / 0.843 / 0.877 and the current 0.977 / 0.837 / 0.840. Nothing
else in the audit is large: absent truth indels are equal between arms (91/66 vs 95/74), and the
calibrated arm already makes fewer dosage errors than counts (hom→het 7/5 vs 16/8).

**Why the posterior favours reference at real hets.** The per-locus prior is 1e-4 for an indel
(het:ref odds 1:20,000) and 1e-3 for a SNP, and the calibrated TAPS profile is pessimistic
because the library is: `taps-error-profile.counts.json` shows 19 % of molecules at
hom-reference 14-mer homopolymers carrying a deletion (39 % at 20-mers), i.e. a gap-open of
Q19 per position for a forward-sum HMM. That calibration is self-consistent (per-position rate ×
L equivalent placements = per-molecule rate), and it is 60× GATK's CONSERVATIVE curve (Q37 at
14 copies ≈ 0.3 % per molecule) and ~10× Octopus's PCR-free table (Q17, but Octopus's SIMD HMM
is best-path, `_mm_min_epi16`, so no sum over placements ≈ 2 % per molecule). At such stutter
rates an AD 9,6 het is only ~85:1 likelihood evidence, so no prior-weighted presence decision
keeps it; the forest keeps it because it is trained on the candidate population. One gap-open
per context cannot serve both decisions: the conservative profile restored presence at 16/16
traced losses but doubled hom→het errors (26/25 vs 7/5), because a REF-looking stutter read at
a true hom-alt then looks like a second haplotype.

**What GATK and Octopus do instead** (sources: `tmp/indel-dosage-mechanism/upstream/`, and
sparse clones of GATK 4.5.0.0 and Octopus made during the review).
- GATK: the heterozygosity prior (1e-3 SNP, 1.25e-4 indel) enters only QUAL = P(AF>0),
  thresholded at 30; the genotype is the argmax of the PLs with no prior
  (`GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN`; posteriors are used only in `--dragen-mode`,
  which also drops the call confidence to 3). Per-read allele likelihood = max over the
  haplotypes carrying the allele (`AlleleLikelihoods.marginalize`), then each site's genotype
  likelihoods are computed independently (`GenotypeLikelihoodCalculator`), reads capped at Q45
  below the best haplotype (`normalizeLikelihoods`, symmetric to reference in 4.5) and
  poorly-modelled reads dropped. Haplotypes come from assembly: only observed alleles.
- Octopus: joint genotype posterior over haplotype pairs with a coalescent prior (θ 1e-3 /
  1e-4), variants reported at marginal posterior ≥ Phred 0.1 and then filtered at QUAL ≥ 10 or
  by its random forest; GT is the MAP genotype, and candidates above threshold that the MAP
  genotype lacks are force-called into viable genotypes. Only observed alleles form haplotypes
  (200, with lagging and holdout); repeat penalties are set per haplotype position.
  Rastair's regional model is Octopus's structure with GATK-magnitude priors and no separate
  presence filter, which is the combination that vetoes.

**Other findings in the regional code.**
- `features::add_likelihoods_with_profile` proposes all four bases at every SNP candidate
  (`Base::KNOWN`) and `regional::infer` protects the best pair of every genotype of every site,
  so unobserved alleles take ~3 of the 64 beam slots per site and regions split at ~21
  candidates (feature histogram peaks at 11 and 21; `max_sites` 32 is never reached). Later
  regions see earlier sites as reference.
- Once the beam is full every candidate runs a full matrix pass over all local fragments and
  all proposals, and `context_profile` builds one `compair::Read` per (read, template) with the
  template's gap tracks in the memo key: 24,500 CPU-s/Mb on the pilot (98,075 CPU-s for 4 Mb)
  against 1,100-3,100 CPU-s/Mb for the per-site PairHMM after the dense-region work and
  15 CPU-s/Mb for counts. This is why chr19 training timed out.
- The regional SNV forest filters more true SNVs than the count forest (174 absent vs 129; 50
  are canonical PASS with median ML 0.99), on 4 Mb of training. The previous per-site PairHMM
  arm, trained on full chr19, has SNV recall 0.972 against 0.961; the pilot cannot separate the
  forest size from the model.
- Everything structural is sound: GATK forward convention, symmetric Q45/MAPQ cap, independent
  padded templates, multi-indel haplotypes, profile serialised with the forest, tests pass.

**Plan, in order.** (A) Decide presence with the forest and dosage with the likelihood: when an
allele passes ML, choose the genotype among those carrying accepted alleles (flat prior, or at
most het:hom 2:1), keep 0/0 only when nothing passes; report the conditional posterior as GQ.
Rerun the 4 Mb pilot; expect every class above canonical. (B) Propose observed alleles only and
protect observed genotypes; rescore incrementally so a candidate's pass touches only pairs that
include it; drop the template-dependent read objects from the memo key. Target the per-site
cost so full-chr19 training fits in `short`. (C) Retrain on full chr19 and validate full chr12
with both dosage directions, SNV isec and Platinum phase. (D) For homopolymers, replace the
per-position gap-open with a haplotype-level stutter mixture (length-specific, −1 ≫ −2) so
presence and dosage share one calibrated model; calibrate on chr19 hom-ref tracts as now.
(E) Give the forest prior-free molecular features (conditional het/hom posterior, per-strand
and per-mate allele consistency) rather than the prior-weighted `haplotype_alt_probability`.

### Regional PairHMM v2: the review's five changes (2026-10-16)

All five recommendations above are implemented on Codex's regional tree; artifacts, job ids,
the frozen binary and the source patch are in `tmp/regional-pairhmm/v2/` (`README.md`).

- **A. Presence is the forest's, dosage is the likelihood's.** `DiploidPosterior::select`
  chooses among genotypes carrying ML-admissible alleles whenever there is one; `0/0` is the
  molecular genotype only when no allele passes. GQ is the dosage error conditional on that
  choice (still floored at the pruning mass). The het:hom prior inside the posterior is 1:1
  per locus, so this is a flat prior over the accepted genotypes.
- **B. Search and cost.** `add_likelihoods_with_profile` proposes observed alleles only (alt
  depth ≥ 2 or passing `pre_ml_filter`); observed sites that are not candidates stay as
  context within `window_radius` (150 bp) of a candidate and emit nothing; fragments are
  borrowed as the index run of their buckets, as the per-site code did. In `pair_hmm`, gap
  qualities come from the reference repeat context at each read's aligned positions
  (`RepeatTracks`, GATK's convention), so a profile is one `compair::Read` under every
  template and the memo serves every template again; the memo key carries only the
  positions whose gaps differ from the plain ones. `TemplateTrack` is back to its key.
- **D. Stutter mixture.** `ErrorProfile` gained a `stutter` table (`StutterRate`: per
  period and minimum copies, per-molecule probabilities of −2, −1, +1, +2 units).
  `error_model::stutter_mixture` builds, for every compared haplotype and every tract one
  of the compared indels sits in (`tract_of`: normalized start, primitive unit,
  `equivalent_end` copies), its stuttered copies at *that chromosome's* copy count,
  deduplicated against the compared haplotypes; `likelihood_matrix_geometry` scores the
  extended set and folds each row as `log Σ w_k exp(L_k)` before the cap and the shrinkage
  estimate (`base_context` keeps `estimate` on the compared haplotypes). With a stutter
  table the per-position repeat qualities are plain. `calibrate_haplotype_errors.py` counts
  whole-unit steps per molecule (`step_minus_2` … `step_plus_2` in `.counts.json`) and
  shrinks them toward the conservative curve's per-molecule rate (`fit_stutter`). On the TAPS
  NA12878 chr19 calibration a 14-mer homopolymer reads −1 in 14 %, −2 in 4.5 %, +1 in 7.6 %
  and +2 in 0.8 % of molecules.
- **E. Prior-free features**: `haplotype_conditional_homozygous`,
  `haplotype_molecular_alt_molecules`, `haplotype_molecular_alt_fraction` (GATK-style allele
  depth from the best allele-carrying pair, counting molecules whose two chromosome
  likelihoods differ by more than log10 0.2). The layout is 16 molecular features; the
  `feature_name_layout` snapshot covers the Standard structs only and did not change.

**Cost** (4 Mb pilot training, 64 threads): Codex's regional 2,202 s wall / 133,369 CPU-s;
v2 per-position profile 343 s / 20,763 CPU-s; v2 stutter profile 445 s / 26,606 CPU-s.
Tests: 715 default + 751 seqair, clippy under both feature sets and the Python tests pass
(`check-*.log`); the full-chr19 training that timed out at 6 h now fits in an hour or two.

**4 Mb pilot** (`tmp/regional-pairhmm/v2/pilot-4mb-report/report.json`; same 4 Mb chr19 forests
and 4 Mb chr12 panel as Codex's pilot; aardvark GT-aware F1 / dosage errors hom→het + het→hom
on exact shared alleles):

| arm | SNV F1 | ins F1 | del F1 | SNV dosage | ins dosage | del dosage | call CPU-s |
| --- | --- | --- | --- | --- | --- | --- | --- |
| canonical (counts) | 0.9742 | 0.8301 | 0.8573 | 22+20 | 16+14 | 8+17 | 46 |
| Codex regional | 0.9764 | 0.8311 | 0.8384 | 14+12 | 9+7 | 7+2 | 98,075 |
| v2, per-position profile (A+B+E) | 0.9767 | 0.8498 | 0.8757 | 14+13 | 11+8 | 12+4 | 14,804 |
| v2, stutter profile (A+B+D+E) | 0.9769 | **0.8566** | **0.8763** | 14+14 | 9+5 | 11+4 | 18,568 |
| previous per-site PairHMM, full-chr19 forest | 0.9826 | 0.8011 | 0.8532 | 20+12 | 38+7 | 27+7 | — |

Both v2 arms beat the count model in every class, and the veto is gone (0-1 PASS `0/0` at
truth variants per class against 22/24/38). The stutter mixture adds a little on insertions
and on dosage for 25 % more CPU; on 4 Mb that is within noise. The per-site arm's SNV recall
(0.972 vs 0.963) remains the forest-size effect: it is the only arm trained on full chr19.
Phase: 0 switches on the 1,141 shared links in every arm. The full-chr19 training and
full-chr12 validation of both v2 profiles run as jobs 195272-195277
(`tmp/regional-pairhmm/v2/full/{perpos,stutter}/benchmark`).

### Merged base: `feature/haplotypings` plus the SNV features, genotype rule and likelihood phasing (2026-10-16)

Branch `feature/haplotypings-merge` (worktree `tmp/merge-trial`) merges `feature/pairhmm_bsb`
into `feature/haplotypings` at `3429c89`. What is kept from each side:

- **Indels are the cluster pair-HMM caller's alone** (`--experimental-indels`, the "Pair-HMM
  indel calling" section): no indel forests, no `normalize_indels`, no canonical anchor
  preprocessing, no `assign_indels`, no count genotyper. `MlFeatureSet` is
  `Standard | Simple | Haplotypes`; `Canonical` and the `InsertionFeatures`/`DeletionFeatures`
  layouts are gone, and a model holds three forests. compair 0.2.0 comes from crates.io: the
  vendored insertion-emission patch was unused once the regional model emitted inserted bases
  at one.
- **SNVs keep the regional molecular features** (`ml train --ml-features haplotypes
  --haplotype-error-profile`; 13 depth-free features appended to the CpG/de-novo/other layouts,
  68/69/67, see "Depth-free molecular features" below),
  the genotype rule (the forest decides presence, the likelihood decides dosage; `GTM`, `GQ`)
  and the **likelihood phaser** (`--phase`, repeatable `--phase-mask`). `ColumnPhase` is
  `Observed | Phased(RecordPhases) | Likelihood(ColumnEvidence) | Genotyped(RecordPhases)`;
  `call_of` reads `Phased` and `Genotyped`. Indel records are not phased: the MEC solver and
  `problem::indel_site` survive as `cfg(test)` reference code (`legacy_phase_segment`), so the
  "Indel sites" gate table above is historical.
- **`--indel-error-profile FILE`** (experimental, INDELS group) is the per-read repeat-quality
  experiment: the indel pair-HMM's reads take insertion-open, deletion-open and extension
  qualities from the calibrated `ErrorProfile.repeats` evaluated on the read's own bases
  (`ErrorProfile::tracks` → `HaplotypeRead::new(.., gaps)` → `compair::Read::new`), in place of
  compair's `PcrIndelModel` (GATK's conservative curve, Q45 open floored at Q20 on the read's
  tandem repeats). `PileupMappingParams.indel_gap_profile` carries the `Arc<ErrorProfile>` to
  `prepare_cluster`; without the flag the output is HEAD's. Benchmark arms are `gatk`
  (default) and `calibrated` in `tmp/regional-pairhmm/merge/` (`checks.sbatch`,
  `build.sbatch`, `bench.py`, `train.sbatch`, `call-*.sbatch`, `compare-*.sbatch`,
  `submit.sh`); the baseline is `tmp/regional-pairhmm/competing/benchmark/chr12/haplotypings`
  and the SNV/phasing reference `tmp/regional-pairhmm/v2/full/`.

Two traps met while checking the worktree:

- **A worktree under the repository's `tmp/` inherits the parent checkout's
  `.cargo/config.toml`**, and cargo concatenates alias arrays across config files, so `cargo
  insta`, `cargo xtask` and `cargo bin` expand to doubled commands (`run --package cargo-bin --
  run --package cargo-bin -- cargo-insta bin cargo-insta …`, failing with "No package found for
  binary run"). Call the installed binary directly
  (`.bin/rust-1.98.1/cargo-insta/1.46.1/bin/cargo-insta test --accept`) or put the worktree
  outside the repository.
- `clippy.toml` sets `allow-indexing-slicing-in-tests = true`, so an
  `#[expect(clippy::indexing_slicing)]` inside a `cfg(test)` module is an "unfulfilled lint
  expectation"; use `#[allow(.., reason = ..)]` there. The `len_without_is_empty` warning on
  `IndelAllele::len`, `regions[0]` in `call.rs` and `usize::from(n.get())` in `emit.rs` are
  HEAD's own and were left alone.

### Depth-free molecular features (2026-10-16)

The forests must work at depths other than the training BAM's (~26x), which is why the
Standard layout holds only ratios, log-ratios and per-read means. The first merged layout
broke that: measured on the chr19 training export of the merged forests
(`tmp/regional-pairhmm/merge/benchmark-16feat/haplotypes/features/`, medians at true
variants across quartiles of informative molecules at the site, Spearman against that count
at negatives), `haplotype_log_bayes_factor` rose 41 → 56 → 65 → 64 (rho −0.53; it is a sum
over molecules), the raw counts `molecular_alt_molecules` 6 → 17, `fragments` 40 → 58 and
`linking_fragments` 12 → 18 (rho +0.4 to +0.7), and `context_sites` sat at the 32-site cap
in every quartile. The bounded posteriors saturate at 1 or 0 and only sharpen with depth;
the fractions (`molecular_alt_fraction`, `expected_fraction`, `realigned_fraction`,
`pruned_mass`) and the Standard `alt_score`/`depth_ratio_alt` were flat (|rho| < 0.1). By
importance the Bayes factor was fourth or fifth in every forest, the counts under 0.007.

The layout is now 13 features (68/69/67 with the Standard parts): the Bayes factor is a
per-molecule mean over the fragments covering the site
(`haplotype_log_bayes_factor_per_molecule`, clamped ±50 after dividing), the counts are
`haplotype_linking_fraction` (fragments covering the site that also reach another site of
the region, over those covering it) and `haplotype_informative_fraction` (fragments whose
likelihoods under the best carrying pair's two chromosomes differ by more than 0.2 decades,
over the region's fragments); `context_sites`, `fragments`, `linking_fragments`,
`molecular_alt_molecules` and the always-one `boundary_resolved` are gone. The layout change
retrains every Haplotypes forest; nothing bundled uses it. The CLI test
`cpg_only_output_keeps_non_cpg_snps_in_the_likelihood_context` probes
`haplotype_linking_fraction` instead of the site count; its read pairs already span both SNVs
(a `Scenario` pair's second read starts at `right`, so `balanced(100, 350, _)` covers 350-450,
and moving `right` to 450 breaks the link). Depth-free features are necessary, not sufficient: the check that the model
transfers is training at 26x and scoring a BAM of another depth against a forest trained at
that depth, which has not been run yet.

### Merged branch: full chr19 training, full chr12 validation (2026-10-17)

Jobs 202903-202909 on commit `25096b9` (depth-free layout), artifacts in
`tmp/regional-pairhmm/merge/` (`benchmark/`, `report/report.json`, `report/indels-*.txt`; the
superseded 16-feature training is `benchmark-16feat/`). Same protocol as the v2 full run: GIAB
HG001 chr19 training with the chr19-calibrated SNV error profile, chr12 calls with
`--experimental-indels --phase` and the centromere/ENCODE masks, aardvark GT-aware on the GIAB
confident BED, Platinum phase assessment. `merge_gatk` is the merged caller as committed;
`merge_calibrated` adds `--indel-error-profile` with the same calibrated profile.

| arm | SNV R / P / F1 | ins R / P / F1 | del R / P / F1 | dosage SNV · ins · del (hom→het + het→hom) | chr12 call |
| --- | --- | --- | --- | --- | --- |
| `feature/haplotypings` | 0.9621 / 0.9855 / 0.9737 | 0.8705 / 0.9414 / 0.9045 | 0.9114 / 0.8983 / 0.9048 | 621+689 · 269+188 · 208+233 | 7 s / 214 CPU-s |
| merged, PCR gap model | 0.9721 / 0.9879 / **0.9799** | 0.8705 / 0.9414 / 0.9045 | 0.9112 / 0.8982 / 0.9047 | 321+535 · 269+188 · 208+233 | 7,601 s / 465,389 CPU-s |
| merged, calibrated per-read gaps | 0.9721 / 0.9879 / 0.9799 | 0.8884 / 0.9016 / 0.8949 | 0.9057 / 0.9075 / 0.9066 | 321+535 · 254+197 · 275+159 | 7,708 s / 471,465 CPU-s |
| v2 per-position (16 features, for reference) | 0.9726 / 0.9873 / 0.9799 | 0.7936 / 0.9598 / 0.8688 | 0.8239 / 0.9558 / 0.8849 | 322+530 · 360+176 · 351+192 | 7,571 s / 463,611 CPU-s |

Phasing against Platinum (masks applied by the assessment): `feature/haplotypings` 59,421
phased SNPs, 87 switches / 39,338 links, 734 phased non-hets; merged 49,245, 8 / 31,496, 327
(both gap models identical, as SNVs are untouched); v2 49,349, 7 / 31,545. On the 29,834 links
all four arms share: 14 / 7 / 7 / 7.

- **The merge is exact on indels**: `compare_indels.py` finds the merged PCR-model arm's
  40,179 indel records identical to `feature/haplotypings` in POS/REF/ALT/QUAL/FILTER/GT/AD/DP
  (phase aside). The one-record differences in the deletion row (correct 11,596 vs 11,599) are
  SNV-side changes to multi-allelic records, not indel calls.
- **SNVs carry over in full**: +0.6 F1 points over the Standard forest, SNV dosage errors
  halved, and the depth-free layout matches the 16-feature v2 arm (0.9799 both; v2 had
  slightly more recall, 0.9726 vs 0.9721, and less precision). The top molecular features are
  now `expected_fraction`, `alt_probability`, `heterozygous_probability` and the per-molecule
  Bayes factor (0.07-0.08); `linking_fraction` and `informative_fraction` are below 0.01.
- **Calibrated per-read repeat qualities do not help the indel caller**: insertions gain recall
  (absent 1,294 → 1,093) but lose more precision (false passes 492 → 999, F1 −1.0 point);
  deletions trade the other way (F1 +0.2) with deletion hom→het errors up 208 → 275; joint
  indel F1 0.9047 → 0.9009. Same conclusion as the SNV-side profile experiments: the TAPS
  library's stutter rates are too high for a per-position gap-open to serve the presence
  margin, and the PCR model stays the default. `--indel-error-profile` remains available.
- **Phasing** keeps the likelihood phaser's accuracy (8 vs 87 switches; 7 vs 14 on shared
  links) and loses ~17 % of phased SNPs against the MEC solver, which bridged blocks through
  indel sites and had no posterior gate.
- **Cost** is the SNV regional model: 2 h 07 m / 465k CPU-s for chr12 against 7 s for the
  baseline, the same as v2. This is the next thing to work on if the merged caller is to be
  the default; see "SNP-dense regions: pre-ML scoring cost".
- Still open: the depth-transfer check (26x forest on another depth) described under
  "Depth-free molecular features".

Operational: a cancelled call keeps executing its binary for a while, so a build job that
`cp`s over `rastair-merge` fails with "Text file busy"; `build.sbatch` now copies to a
temporary name and renames.
