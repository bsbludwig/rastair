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

`Pileup` objects are constructed in `src/call/pileup/from_hts.rs` via `Pileup::from_hts()`.
The `reference_base` comes from the FASTA sequence via `Base::from(u8)`, which maps any non-ACGT byte (e.g. `N`) to `Base::Unknown`.
There is **no filtering** to skip pileups with zero reads or Unknown reference bases before `PileupMetrics::new()` is called.

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

## Generated docs are partly machine-dependent

`cargo run -- internal cli-docs docs/src/cli.md` bakes `available_parallelism()` into the
`--threads` default, so regenerating on a machine with a different core count produces spurious
diffs. Check `git diff docs/src/cli.md` and revert those hunks. `internal vcf-docs` currently
also renders `Phred` as `Phred>`, which the committed file does not have.

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

## Indel calling: two pathways, and the invariants they share

`--experimental-indels` runs a hard-filter chain (`src/call/variant_calling/indel_calling/hard_filters.rs`);
`--experimental-indels-ml` runs the ML model instead, degrading to the hard-filter chain under `--no-ml`.
Both consume the same `IndelCounts` from `aggregate_indels()` in `src/metrics/pileup_metrics.rs`,
so the invariants below apply to both.

**Everything is per fragment, not per alignment.** Indel bookkeeping (`IndelReadData`) is built in
lock-step with `raw_reads` in `src/call/pileup/from_hts.rs` and the same overlap deduplication
`swap_remove`s both. `ref_count = reads.len() - total_indel_reads` is only a real partition because of
this. Any per-alignment property that differs between mates — reverse flag, position in read, which
mate carried the base — is *not* recoverable after that collapse.

**Strand is OT/OB, never the reverse flag.** Both mates of a fragment share an OT/OB assignment but
have opposite reverse flags, so OT/OB is the only strand notion invariant under deduplication. A
reverse-flag split collapses to whichever mate survived (always the leftmost, hence forward-flagged).

**The strand-bias test is the binding constraint on indel recall, and it does not pay for itself.**
Measured on chr12 (BED-restricted, normalized), turning it off with `--indel-strand-bias-alpha 0`
moves the hard-filter pathway from **P 0.9819 / R 0.7485 / F1 0.8495** to
**P 0.9747 / R 0.9214 / F1 0.9473**: it was rejecting **4,288 true indels to remove 252 false ones**,
17:1 against. That single flag beats the pre-fix `d4244fcc` (F1 0.9326) at far higher precision.

The test is well constructed; the hypothesis it tests is false. Its null is "the supporting fragments
are drawn from the locus' own strand mix", and for TAPS indels they are not — OT and OB reads present
different sequence after C→T conversion, so genuine indel support is strand-asymmetric for reasons
that have nothing to do with artifacts. The p-value is still *informative*, so the intended home for
it is an ML feature, not a hard gate.

Beware that this filter masks everything downstream of it. The `--indel-noise-exclusion`,
`--indel-error-rate` and `--indel-het-vaf` sweeps below all saturate at F1 ≈ 0.85 **because they were
measured with the strand test on**, throwing the recall away before those knobs could act. Re-measure
anything in that table before trusting it at a lower alpha.

**Strand bias is a test, not a rule.** `IndelAlleleCounts::strand_bias_p_value` is a two-sided exact
binomial test against `IndelCounts::null_ot_fraction` — the strand mix of the *non-supporting*
fragments at that locus, smoothed with one fragment of prior mass so an all-one-strand background
(routine at low depth) does not produce a degenerate null. That prior is centred on the **locus'** own
OT share, not on 0.5: a homozygous indel has no non-supporting fragments at all, and a flat 0.5 prior
then judges it against balanced coverage it never had, rejecting every hom-alt call at a one-strand
locus for a skew that belongs to the coverage. Rejection is at `--indel-strand-bias-alpha`
(default 0.05).
Do **not** replace this with an `ot > 0 && ob > 0` rule: at the default `min_indel_ao` of 2 that rejects
a 2/0 split, which is a coin flip, and measured on real data it tracked the chance rate almost exactly
below AO≈4 — rejecting ~26% of all indel candidates, ~85% of them by chance.

**"Repeat" means a real repeat — and that, not the noise sidedness, is the lever.** `terminal_repeat`
once used a single 3 bp window for both periods, which made the period-2 arm
`seq[0]==seq[2] || seq[len-3]==seq[len-1]` — true for **43.75% of random reads**.
`TerminalRepeatLimits` now counts whole repeat units per period (homopolymer ≥4 bp, dinucleotide
≥3 units), ~3-4% on random sequence; `has_repeat_seq` guards `units < 2` because a sub-two-unit
window makes the periodicity `all()` vacuously true. This changes `has_repeat`, an ML feature, so it
deepens the retrain need.

Which side of the ratio noisy fragments come off is `--indel-noise-exclusion`
(`symmetric` default / `ratio-only` / `depth-only`). The theory says `symmetric`: noise is a property
of the read, so it lands on supporting and non-supporting fragments alike, and a one-sided haircut
inflates VAF by the noise rate with the binomial flipping hom-ref→het at VAF ≈ 0.218. **Measured, the
choice is nearly inert** — chr12, GIAB-BED-restricted and normalized: symmetric F1 0.8509,
ratio-only 0.8510, depth-only 0.8561. Do not spend time here.

**The hom-ref/het boundary is not where the recall is. Stop tuning it.** The pre-fix `d4244fcc`
scores **P 0.9194 / R 0.9462 / F1 0.9326** against HEAD's **P 0.9823 / R 0.7505 / F1 0.8509** (chr12,
BED-restricted, normalized), and essentially all of the gap is recall. There are three knobs that
move that boundary, and *all three saturate around F1 0.85*:

| knob | best | at |
| --- | --- | --- |
| `--indel-noise-exclusion` | 0.8561 | `depth-only` |
| `--indel-error-rate` | 0.8545 | 0.03 (default 0.05 → 0.8509) |
| `--indel-het-vaf` | 0.8527 | 0.40 (default 0.5 → 0.8495) |

Each buys recall by giving up more precision than it gains, and none gets near 0.93. `--indel-het-vaf`
is the most defensible of the three — reads carrying an indel are harder to place than
reference-matching reads, so a genuine het indel *is* observed below 0.5, and modelling that is
better than reaching the same boundary via `--indel-error-rate`, which also weakens the hom-ref
hypothesis against genuine sequencing noise. But defensible is not the same as effective: it is worth
+0.003 F1. Leave it at 0.5 unless you have a reason beyond F1.

The conclusion this forces is that HEAD's *candidate set* is smaller, not merely genotyped more
strictly — the missing true indels are being lost before genotyping, so look at the depth gate, the
min-AO gate, the strand test, and the read-level indel filters in `from_hts.rs`, not at the binomial.

The ML pathway deliberately keeps its own one-sided `depth_offset` so those distributions do not
move further.

**Known gap: mate discordance.** `resolve_pair` picks the surviving mate from the anchor base and
`second` flag, blind to indels, so when the two mates' CIGARs disagree the fragment's vote is decided
arbitrarily and the loser's evidence is dropped. Measured at ~1,562 lost votes per 5 Mb, **97% of it
genuine mate-level alignment disagreement** and only 3% asymmetric read filters. Do not "fix" this by
OR-ing the mates' observations — that promotes alignment ambiguity to positive support in exactly the
repeat contexts where indel artifacts live. The intended fix is to make discordance explicit (abstain,
and carry it as a feature).

**Indels obey the same emission contract as SNVs.** Every genotyped allele is *built* with a FILTER
(`indel_strand` / `indel_hom_ref` / `indel_no_ml`) rather than dropped in the caller, but `to_vec`
emits only PASS by default, all of them under `--all`, and none under `--cpgs-only`. Emitting them
unconditionally — the earlier behaviour — put tens of thousands of `indel_hom_ref` lines into a plain
WGS run. `indel_no_ml` exists because the ML path's `FORMAT/ML` is empty both when the model declined
to score and when it passed, so PASS alone could not distinguish them.

`Pileup` carries **no `#[serde(default)]`**: `.mpk` encodes structs as positional arrays, so a field
added anywhere but the end shifts everything after it and a default decodes neighbouring data instead
of failing. `MPK_FORMAT_VERSION` (`src/io/mpk/format.rs`) is the compatibility mechanism — bump it
when `Pileup`/`PileupMetrics` change.

Measurements behind all of the above, and how to reproduce them, are in
`.claude/notes/indel-strand-concordance-vs-fragment-dedup.md`. Validate changes with
`rastair verify --experimental-indels --truth <GIAB HG001>` — the sample in `tmp/taps/` is NA12878.
End to end, with the truth set and reference already on disk:

```sh
rastair call tmp/taps/NA12878_aa_chr12.bam \
  -r tmp/na12878/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -l chr12 --experimental-indels --no-ml -o out.bcf
bcftools index -f out.bcf
rastair verify --truth tmp/1000genomes/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  -l chr12 -R <GIAB high-confidence BED> out.bcf --experimental-indels
```

Tune on one chromosome, confirm on another before believing a threshold — chr20 is the conventional
GIAB validation chromosome.

`verify` loads **only PASS records**, so emission-policy changes do not perturb the comparison. It
reports per-category FN/recall/F1 for the indel categories (they are classified by allele length, so
the truth VCF bins them identically; the CpG/DeNovo/Other categories depend on rastair INFO flags the
truth set does not carry, so those stay precision-only). Baseline a change by building the comparison
commit in a throwaway worktree (`git worktree add --detach /tmp/rastair-baseline HEAD`) rather than
stashing — the working tree stays intact and both binaries can run back to back.

**Always pass `-R/--regions-file <GIAB high-confidence BED>`.** Without it `verify` scores calls in
regions the truth set makes no claim about, and every one of them lands as a false positive: 40.6% of
our PASS indels on chr12 fall in the 7.7% of the chromosome GIAB excludes. Unrestricted scoring
reported P 0.618 for a caller that is actually at P 0.982, which is what made an internally-reported
hap.py F1 of 0.814 look irreconcilable with our own numbers. `verify` also normalizes indels to
minimal representation before matching (`minimal_representation` in `src/verify.rs`), so
`100 TC>TCAA` and `99 C>CAA` are the same variant. Quote both scorings if you quote either — the
whole-chromosome numbers are still a useful relative signal between two builds, they are just not
comparable to anything external.

Matching is on `(chrom, pos, ref, alt)` only — `verify` is **genotype-blind**, so it is more lenient
than hap.py, which counts a GT mismatch as both an FP and an FN. That bias is currently unquantified.

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
  `InsertionFeatures` (23), `DeletionFeatures` (23). Each has an `extract()` returning the
  struct; `FeatureCalculator::calculate_*` wraps `as_row()` into an `Array2`.

**Feature order is frozen by every trained model.** Reordering a field silently corrupts
predictions. Two tests guard this in `features.rs`: `feature_counts_are_stable` (pins the
counts) and `feature_name_layout` (an insta snapshot of every `name→index` mapping —
this replaces the old "never change the order" comments; update it via `cargo xtask insta`
only after verifying a layout change is intentional).

A third guard is at runtime: `check_feature_widths` (`src/call/ml.rs`) refuses a model whose
`ForestMeta::n_features` disagrees with what the build extracts. This has to be an error, not a
warning — `FlatForest::predict` does no width check, so a stale forest fed a wider row reads
whichever columns its split indices name and predicts confidently from the wrong features. It is the
only way a model file can fail while still producing plausible output.

**Always train with `-R <high-confidence BED>`.** Training labels every candidate not in the truth
VCF as negative, and outside a truth set's high-confidence regions the VCF asserts nothing — so a
real variant there is absent from it and gets taught to the model as a *negative*. Those mislabelled
examples are not spread evenly: high-confidence BEDs exclude repetitive, low-mappability sequence,
which is where indels concentrate. Restricting to the GIAB HG001 BED over chr1/chr6/chr11 drops
**~40% of indel training examples, ~89% of them negatives**, against ~10% for the SNV sets (in line
with the 7.7% of sequence excluded). Candidates outside the regions are dropped, not relabelled.

**But apply it to the indel models only.** The same restriction measurably *hurts* cpg/denovo/others:
chr12 overall F1 .9884 -> .9877, chr20 .9896 -> .9882, with false positives roughly doubling in every
SNV category on both chromosomes (chr12 Other 842 -> 1,550 FP; DeNovo 151 -> 368). Recall improves,
precision drops more. The asymmetry is in what gets dropped: of the out-of-BED examples removed,
**11.7% of the indel ones were positives against 0.4% of the SNV ones** (cpg: 712,067 removed, 3,037
positive). Out-of-BED SNV candidates are overwhelmingly genuine sequencing error — correctly labelled
negatives, and valuable hard ones from the worst sequence, so dropping them starves the model and it
turns permissive. Out-of-BED indel candidates are dominated by real variants GIAB could not confirm.
Same operation, opposite sign. Train indels on the BED, keep the SNV models genome-wide, and let the
splice below combine them — that is what the splice is *for*, not just a blast-radius convenience.

The BED restriction was worth more than everything else tried on the indel ML path combined:

| chr12 / chr20 | P | R | F1 |
| --- | --- | --- | --- |
| ML, trained genome-wide | .9948 / .9945 | .8066 / .8045 | .8909 / .8895 |
| ML, trained on the BED | .9843 / .9820 | .9222 / .9214 | **.9523 / .9508** |
| hard filters | .9747 / .9736 | .9214 / .9183 | .9473 / .9452 |

Holdout Brier roughly halves for both indel models (.1119→.0556, .1034→.0576) and the forests get
*smaller* (max_tree_size 6509→3501), because they stop spending depth memorising labels that were
guesses. With this, **the ML indel path overtakes the hard-filter path** — same recall, better
precision. Two caveats: a BED-trained model is only validated inside those regions but still scores
everywhere at call time, so out-of-BED behaviour is unmeasured (and unmeasurable against GIAB, for
the same reason the labels were unusable); and these figures are in-BED performance specifically.

**Retraining for an indel change: splice, do not ship the whole retrain.** `ml train` retrains all
five forests, so a change confined to the indel feature vector would otherwise also replace the
methylation and SNV models — much wider blast radius, and unmeasured unless you separately validate
methylation. `examples/splice_indel_models.rs` takes cpg/denovo/others from `--old` and
insertion/deletion from `--new`; `--old` accepts either a current (11-field) or a pre-indel
(7-field) model file. Confirm the result with `examples/model_stats.rs`: the kept and spliced
forests should show the tree sizes of their respective source files.

A symptom to recognise: replacing the cpg/denovo forests fails
`test_adjacent_denovo_cpgs_dual_role_middle_position`, which used to run the bundled model live over
a five-read synthetic pileup. It now forces its calls with `set_pass`/`reprocess`, so it pins the
de-novo logic rather than the current model's opinion of made-up data. If a test like that starts
failing after a retrain, decouple it rather than editing the expected values to match — the numbers
it asserts are the behaviour under test.

Feature names flow to training output via `FeatureCalculator::feature_names() -> FeatureNames`.
`train.rs` uses them for the `--feature-analytics` importance CSVs (`index\tfeature\timportance`)
and the `--export-features` TSV headers, so both exports agree by construction.

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
