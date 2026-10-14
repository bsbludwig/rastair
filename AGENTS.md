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

## Molecular ML feature experiment (2026-10-14)

`ml train --ml-features haplotypes` appends 13 features to all five frozen Standard layouts (68/69/67/30/30). Model metadata automatically enables full read capture and pre-ML PairHMM scoring in both backends; final `--phase` calls are never training inputs. `canonical` is the original-layout ablation with the same canonical indel preprocessing. The bundle is replaced only after held-out validation; the current experiment's full PairHMM results are pending.

`src/call/haplotypes/features.rs` owns shared likelihood parameters and preprocessing for training and calling. Simple indels are normalized before filtering/label lookup; reference tracts, per-molecule depth and query boundary metrics use that anchor. Reads retain their original CIGAR profiles. A boundary outside the preceding matching run or read mask is unresolved. Keep already-accepted indel allele observations even when the canonical boundary is lower-quality/missing: encode that uncertainty rather than imposing a second, different BQ cutoff. A supported SNP within an observed indel shift tract makes the complex chromosome hypothesis unsafe without assembly, so its molecular features are unavailable. Local ML alignment geometry includes the full reference-equivalent placement span; the diagnostic retains its existing geometry. Indels over 100 bp and wholly aligned-fallback indels have unavailable molecular features.

Training and calling use the same pre-ML eligibility predicate, de-novo adjacency setup and shared CpG shrinkage. Use normalized, split GIAB truth and remove unused ALT rows after splitting (`GT="alt"`). A deletion's full reference span must be inside the confidence BED. New molecular feature extraction fails the training run on any failed segment. Relative-weight normalization in `posterior_pairs` prevents large common negative log offsets from producing aggregate probabilities above one.

`tests/call_cli.rs::model_metadata_activates_pair_hmm_features_without_a_phase_flag` uses a tiny forest that responds to the availability feature; masking makes that feature unavailable. It runs on both readers and proves the model selected the actual capture/scoring path. `features.rs` and `reads.rs` tests cover canonical boundaries/qualities, shifted CIGAR placements, read-group identity, molecular deduplication, aligned fallback and CpG conversion guards.

`scripts/benchmark_haplotype_ml.py` provides separate train/call/compare stages. Run long experiments via SLURM `short` (24 h maximum), supplying realistic `-t`, `-c`, and `--mem`. Current jobs request 64 CPUs/96 GB/6 h for chr19 training, 64 CPUs/96 GB/12 h for full chr12 calling, and 8 CPUs/16 GB/1 h for comparison; record actual resources to refine these estimates. Keep binaries, source patches, masks, commands, sample names and complete-status markers in the run directory. Do not install a model from incomplete data or a failed job.

Compute nodes do not share the login/session `/tmp`: put test helper executables in shared storage and add that path to the SLURM job's PATH. GIAB's sample name is `HG001`, not the Platinum Genomes `NA12878`; read it from the VCF header. The augmented reference has IUPAC `Y` bases that Rastair renders as `N`; these can fail bcftools normalization outside GIAB confidence. Subset query calls to the confidence BED before normalization for benchmarking, retaining strict reference validation inside that BED. Do not use a permissive REF-repair flag to hide a coordinate error.

**Important:** Whenever you learned something new about how to develop features, find code, or how to debug issues, you **must** add it to this document.
This is the single source of truth for how to work on this codebase, and it must be kept up-to-date with any new insights or changes.
If you find yourself asking "How do I do X?" and you figure it out, add a section here so that the next person doesn't have to ask the same question.
