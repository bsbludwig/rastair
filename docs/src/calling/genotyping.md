# Genotyping

Rastair uses @ML scores to determine which @altAllele:pl represent true genetic @variant:pl, then applies statistical methods to estimate zygosity.
Rastair only supports diploid genotyping (two chromosome copies) at this time.

## Multi-Allelic Sites

When multiple alternate alleles are present at the same genomic position, they are combined into a single VCF record with comma-separated ALT values.
For example, with `REF=A` and `ALT=T,G`, the genotype represents the complete set of observed alleles:

- A genotype of `1/2` means "on one chromosome there was a T (allele 1), on the other was a G (allele 2)"
- A genotype of `0/1` means "on one chromosome there was the reference A (allele 0), on the other was a T (allele 1)" (G alt did not pass filters)
- A genotype of `1/1` means "both chromosomes had a T (allele 1)" (only T alt passed filters)

This representation follows VCF specification standards where allele indices start at 0 for the reference, and 1, 2, 3... for the alternate alleles in order.

## Genotype Calls

- Homozygous Reference (`0/0`): no alternate alleles passed filters
- Heterozygous (`0/1`): One alternate allele passed filters but read counts show a mix of reference and alternate reads
- Homozygous Alternate (`1/1`): an alternate allele passed filters and read counts show predominantly alternate reads, consistent with two copies of the variant
- Compound Heterozygous (`1/2`, `2/3`, etc.): Multiple alternates passed filters and read counts support different variants on each chromosome copy

## Confidence Scoring

Confidence values reflect call certainty:

- For `0/0` calls: based on the margin between the ML threshold and the highest-scoring alternate
- For variant calls: based on how well read count ratios match the expected distribution

## Biological Interpretation

A genotype represents the complete set of observed alleles in an individual.
All alternate alleles listed in a VCF record's ALT field are simultaneously present at that position in at least one chromosome copy.
In most cases:

- If only one variant passes filters, the genotype indicates whether it's present on one chromosome copy (heterozygous `0/1`) or both copies (homozygous alternate `1/1`)
- If multiple variants pass filters, they represent compound heterozygosity where different variants are present on different chromosome copies (e.g., `1/2`)

## Strand-Specific Counting

For `C→T` and `G→A` variants, only one strand is used to avoid confounding with methylation.
For all other variant types, both strands contribute to genotyping.

## Phasing

`--phase` can be used to group heterozygous @SNV:pl that share a read pair
into phase blocks and written with a `|` separator and a `PS` field:

```
chr12  20028875  .  A  G  ...  GT:PS  0|1:20028875
chr12  20028894  .  C  T  ...  GT:PS  1|0:20028875
```

Both records carry `PS=20028875`, the position of the block's first heterozygous site, so they belong to one block.
Within a block, the alleles written first all sit on the same chromosome copy: above, the reference `A` and the alternate `T` are on one copy, the alternate `G` and the reference `C` on the other.
Each block opens with its lower genotype allele first (`0|1` for a reference/alternate genotype, or e.g. `1|2` for two alternates); nothing relates one block to another, and heterozygous sites no read pair links stay unphased without `PS`.

Phasing now uses fragment likelihoods and a local PairHMM, conditioned on the two alleles of each **called** heterozygous SNP. It changes chromosome assignments only; variant calls, unphased genotypes, genotype scores and methylation calls retain their existing calculation. The previous signed-quality graph solver is used only as a reference in unit tests.

Full sequences, qualities, CIGAR coordinates and OT/OB orientation are captured before seqair recycles its read store. Read name plus read group identifies a molecule, and boundary alignments fetched again in byte-budget subdivisions are captured once. Overlapping mates form one consensus observation: agreement retains the better quality, while disagreement reduces it by the competing quality. Disagreeing CIGAR paths or reference skips use the aligned likelihood fallback.

A site needs at least two independent, quality-qualified molecules supporting **each** called allele before it can provide a genetic link. Consensus base and source mapping quality must be at least 20, and two molecules must link an edge. Converted-strand `T` at a C/T CpG and `A` at a G/A CpG cannot establish these links. Once a block is selected, its likelihood uses all eligible captured fragments covering its sites, including conversion-confounded observations and fragments reaching only one site. Adjacent candidate C/G partners are included together for CpG context even when they lack a genetic link.

Connected chromosome hypotheses expand progressively, with a beam of 32 phase states and 16 held-out states reconsidered on the next extension. Competing orientations at every selected site are protected; protection can enlarge the beam to 64 states. Each state pairs a chromosome with its complementary called alleles, so no alternate unphased genotypes are introduced. Retained states have equal phase priors. A nearby, independently supported CIGAR indel of at most 20 bp can be included as one nuisance allele, marginalizing its genotype and phase without emitting phased indels.

The neutral pilot and final PairHMM scores use the shared CpG shrinkage model described below: default conversion is 0.95, background conversion 0.003, with Beta(2,2) shrinkage and leave-one-fragment-out rate estimates. Homozygous SNP calls are applied to both chromosome backgrounds. Candidate sequence adjacency after SNP and indel edits determines reference and de-novo CpGs.

Each mate is aligned within its own CIGAR-coordinate footprint plus a 12 bp realignment flank; concordant overlaps use the joint molecular footprint. An unsequenced gap between mates provides no alignment observations. Scores are reused for chromosome hypotheses with identical DNA over these footprints, retaining an additional base for CpG partners outside a scoring window. Alignment path sums use row scaling to prevent underflow; their normalization is cached by sequence lengths and gap parameters.

A genetic edge is emitted only if the retained-state posterior probability of its chosen relative phase is at least 0.99. Ambiguous edges can split a connected inference region into separate output phase sets. These probabilities are conditional on called genotypes, retained hypotheses and estimated shared rates, and are not externally calibrated phase quality scores. `GL` and `GC` remain genotype scores; only `GT` orientation and `PS` are changed.

Use repeated `--phase-mask FILE` arguments to exclude centromeres and poorly mappable intervals from **phasing only**. Inputs can be plain, gzip or BGZF BED files, including the UCSC centromere table with its `Chr Start End GapType` header. Coordinates are zero-based and half-open; overlapping intervals from all files are merged. Masked positions retain their variant and methylation calls but have no phased genotype or `PS`. Masks are barriers even when a read pair bridges them: inference, CpG context, nearby indel proposals and PairHMM flanks are restricted to each unmasked interval. Paths must be supplied explicitly, and contig names must match the BAM/reference.

```bash
rastair call --phase --phase-mask centromeres.txt --phase-mask blacklist.bed.gz \
  --fasta-file reference.fa reads.bam --vcf calls.vcf.gz
```

```admonish warning
Phasing requires the `experimental-seqair` backend and cannot be combined with `--cpgs-only` or BED output alone. Inference is bounded to 32 heterozygous sites and a 2 kb span. Blocks also stop at processing-segment boundaries; overlap sites contribute evidence, but each phase set is anchored on its first emitted site. The production phasing model currently assumes TAPS chemistry with the defaults above.
```

Runtime can vary substantially with local variant density, coverage and nearby indels. Use `scripts/benchmark_phasing.py` to compare unphased, PairHMM and an optional archived legacy binary on the same BAM, reference and regions. It records caller wall time, CPU time, peak memory and phase-block statistics, and verifies that every non-phase VCF field is identical. The caller timing excludes the script's subsequent VCF inspection. The HG001 GIAB benchmark is unphased and therefore cannot measure switch-error accuracy.

The benchmark accepts repeated `--phase-mask` arguments or a `phase_masks` list in its input JSON. These apply to the new phasing arm; archived legacy binaries do not support the mask flag. Accuracy comparisons should exclude masked intervals from both outputs.

`scripts/assess_phasing.py` compares called SNP phase with independent phased truth, such as the [Illumina Platinum Genomes NA12878 calls](https://github.com/illumina/platinumgenomes). Supply its confident-region BED and any masks with `--exclude-bed`. Genotype errors are reported separately; switch errors are orientation changes between adjacent genotype-matched phased SNPs sharing both query and truth phase sets. A chromosome-wide orientation flip is irrelevant, while an isolated erroneous SNP orientation contributes two switches. Singleton intersections cannot assess phase. Report phase coverage alongside error rates, since conservative thresholds can improve the rate by leaving difficult links unphased.

## Experimental haplotype likelihoods

The `ml haplotypes` diagnostic calculates SNP likelihoods before final variant or genotype calls. It requires the `experimental-seqair` backend:

```bash
cargo run --features experimental-seqair -- ml haplotypes \
  --fasta-file reference.fa reads.bam --region chr12:20028995-20030086 \
  --output snps.tsv --pair-output pairs.tsv
```

Use `--non-taps` for ordinary WGS reads. For TAPS, `--conversion-rate` defaults to `0.95` and `--false-conversion-rate` to `0.003`. These describe conversion of methylated and unmethylated cytosines respectively; they are separate from base-quality sequencing errors. The lambda spike-in is not used to fit these parameters automatically: incomplete lambda conversion can underestimate conversion efficiency by 2–5 percentage points.

Every focal position is evaluated against all four genetic bases. Shared per-CpG shrinkage is the sole methylation-rate model. Both mates share a chromosome assignment, and overlapping bases form one molecular consensus observation and share a conversion event. There is no depth cap; configured flag, quality and read-end filters still apply.

The default `--likelihood-engine aligned` conditions on the input CIGAR alignment. A conversion-neutral pilot gives C/T identical emissions on OT at potential CpG C positions, and G/A identical emissions on OB at potential CpG G positions. Potential positions are the union across candidate sequences, including de-novo and indel-created CpGs. Other genetic alleles and opposite-strand observations remain informative. The pilot is internal; `--methylation-model` and the integrated-model path have been removed.

The pilot supplies genotype and molecular chromosome probabilities for estimating one methylation rate per candidate CpG, shared across haplotypes that contain it. A symmetric Beta prior shrinks low-information rates toward 0.5 (`--methylation-prior-strength 4`, corresponding to Beta(2,2)). Base and mapping quality enter a separate one-dimensional rate estimate. Only molecules with pilot CpG-carrier probability at least `--methylation-min-cpg-probability 0.9` can fit the rate. Each fragment is scored with its own conversion observation's contribution removed. Estimates use CIGAR coordinates, including insertion coordinates; ambiguous orientations do not fit rates. Candidate CpGs without usable observations retain the prior mean. These are approximate plug-in scores and do not model allele-specific methylation. `--quadrature-order` controls the one-dimensional rate estimate, not a joint haplotype-rate integral.

No CpG-count limit removes neighboring SNPs. `--non-taps` bypasses conversion and rate estimation. Conversion defaults remain 0.95 and 0.003; the spike-in is not used to estimate them. No benchmark truth alleles or WGS genotypes enter TAPS inference.

### Bounded haplotype expansion

The initial candidate radius is 150 bases in either direction (`--window-radius 150`). Fragment connections can extend the context to `--expansion-radius 1000` bases in either direction. A 150-base read cannot necessarily reach 150 bases from a SNP it covers: reach depends on where the SNP lies within that read. Paired reads and chains of overlapping fragments can link positions farther apart than one read's length.

A neighbor can extend the region when at least two independent molecules provide quality-qualified, conversion-unconfounded alternate evidence there, and at least `--min-linking-fragments 2` molecules link it to an active site. Linking bases require consensus base quality and source mapping quality at least 20. Converted observations at biallelic C/T or G/A CpGs on the confounded strand cannot establish a phase connection. Reference-only and conversion-only candidate positions do not bridge otherwise disconnected genetic regions. Candidate SNPs that jointly determine one CpG are also included together, including latent converted C/G candidates and partners adjacent after a deletion. At least two quality-qualified molecules covering either partner suffice for this context dependency; it does not supply a genetic bridge to distant variants. Joint partner hypotheses prevent a genetic A/T at one base from falsely implying a de-novo CpG at an uncovered neighbor. All passing fragments covering the selected loci remain available for final scoring, including observations unsuitable for establishing links.

The tree starts with every focal A/C/G/T hypothesis and, for PairHMM, the eligible nuisance indel. Sites within the initial radius are considered first, followed by connected extensions. Each extension branches the retained and buffered haplotypes over the next site's alleles. When proposals exceed `--max-haplotypes 64`, their fragment likelihoods, shared rates and diploid-pair posteriors are recomputed. Haplotype rank is its posterior chromosome occupancy summed across diploid pairs. The reference haplotype, the best pair for every focal genotype, and a carrier of every retained site's allele are protected. If protection cannot fit the budget, extension stops and reports the excluded sites.

Up to `--haplotype-holdout 32` discarded candidates are buffered and reconsidered at the next extension. This provides bounded recovery, not exhaustive backtracking. `--max-sites 32` is a safety limit rather than the former three-SNP enumeration limit; it accepts 1–64. The scoring pool is bounded by four times the retained-plus-holdout budgets. The single nuisance-indel limit remains.

Each iteration scores the current molecular observations once; it does not multiply overlapping-window likelihoods or reuse the old posterior as a prior. A genotype's prior is divided among its surviving phase assignments after pruning. Posteriors are conditional on the retained candidate set. `maximum_pruned_mass` reports the largest pair-posterior mass discarded at one expansion step, conditional on that step's proposals; it is not an upper bound on all omitted genomic hypotheses.

`SnpEvidence::genotypes_at` projects any selected site from the same joint diploid posterior. These are focal-rooted inference contexts for SNP diagnostics. Overlapping focal contexts are not reconciled into a single chromosome-wide posterior or emitted as VCF phase sets. Production calling, ML forests and the separate `--phase` solver remain unchanged.

### Local PairHMM realignment

Use `--likelihood-engine pair-hmm` to sum over alternative alignments:

```bash
cargo run --release --features experimental-seqair -- ml haplotypes \
  --fasta-file reference.fa reads.bam --region chr12:20028995-20030086 \
  --likelihood-engine pair-hmm --output snps.tsv --pair-output pairs.tsv
```

This scalar forward algorithm has match, insertion and deletion states, following the [GATK PairHMM model](https://github.com/broadinstitute/gatk/blob/master/docs/pair_hmm.tex). `--gap-open` defaults to Phred 45 (`10^-4.5`) for each gap type and `--gap-extend` to `0.1`. Match emissions combine conversion and sequencing error; overlapping mates contribute one consensus base to a single conversion draw. Insertion emissions are uniform over consensus bases, including shared insertions counted once. The alignment sum is normalized by its path mass, conditioning on the profile's length and containment in the local window. Unlike GATK's omitted insertion emission, this gives a normalized observation model: uninformative bases cannot favor a longer haplotype.

The realignment interval spans the selected SNPs and nuisance indel, padded by `--realignment-flank` (12 bases by default). Read sequence is cropped using its CIGAR coordinates; insertion bases are retained, soft clips are excluded, and masked or low-quality interior bases become uninformative observations. This is a local realignment prototype, not assembly of full active regions. Reads spanning a deletion at the target still contribute to the alignment sum.

Shared CpG shrinkage uses this same interval for candidate CpG context and rate estimation, including when a nuisance indel extends beyond the SNP window. One additional partner base at either edge can establish CpG context, but observations outside the scoring interval cannot fit a methylation rate.

CIGAR insertions and deletions passing the read filters seed competing genetic sequences, without using called genotypes or truth alleles. Proposals are left-normalized against the local reference. The nearest eligible simple indel is retained as a biallelic nuisance site. Indels outside the initial radius require at least the configured number of fragments spanning both the focal position and the indel anchor, with a default maximum length of 20 bases (`--max-indel-length`). Deletions covering a selected SNP and other indel proposals are omitted; `omitted_indels` reports their count. Its genotype prior is `1 - --indel-prior` for reference/reference and `--indel-prior / 2` for reference/indel and indel/indel (default prior `0.0001`). Compatible phase assignments share the joint SNP/indel genotype prior. The SNP probabilities marginalize over the indel's presence and phase, rather than committing to an indel call first. Candidate sequence adjacency determines CpGs after either SNPs or indels create or destroy them.

Mapping hypotheses reuse each profile's score when its contributing alignments have the same mapping state. Fully unmapped profiles have uniform emissions, so their alignment sum cancels against the path-mass normalization and is evaluated directly.

For concordant overlapping CIGAR paths, mates form one molecular profile with a consensus observation at each shared base. Disjoint mates remain separate profiles assigned to the same chromosome. Discordant overlapping CIGAR paths and reference skips use aligned SNP emissions, with the same base-consensus policy and the candidate's edited CpG context; deleted partner positions are skipped, and indel alignment evidence is not scored. `realigned_fragments` and `aligned_fallback_fragments` make that fallback visible. A distant mate outside the scoring interval contributes no mapping event, sequence likelihood, orientation evidence or reference-skip fallback, keeping results independent of input chunks.

Shared shrinkage also estimates rates for CpGs in realignment flanks, using the same candidate sequence adjacency and authoritative scoring interval. Richer indel candidates, allele-specific methylation and full-read assembly remain future steps. These bounds and fallback behavior need to be considered when interpreting diagnostic probabilities.

`snps.tsv` contains alternate-presence, heterozygous and homozygous-alternate posteriors, expected allele fraction, genotype entropy and natural-log Bayes factors. The window and site-only Bayes factors compare the expanded context, including the nuisance indel, with the focal-site calculation. `pairs.tsv` contains every candidate diploid haplotype pair's natural-log likelihood, prior and posterior. Its `haplotype_a/b` strings give SNP bases at the listed 1-based positions; `indel_a/b` give an unanchored `0-based-start:deleted>inserted` replacement, or `.` for reference. The SNP prior assigns `1 - --snp-prior` to the reference genotype and distributes the remaining mass uniformly across the other genotypes at each site; compatible phase assignments share that mass equally. These model probabilities are not empirically calibrated ML call scores.

By default, focal candidates come from observed alternate bases and potential de-novo CpG alleles whose converted base resembles the reference. Potential CpG partners include the reference, neighboring observed SNP alleles and, for PairHMM, indel proposals. Thus an indel-created CpG can seed a C/G hypothesis even when conversion leaves only reference-looking T/A observations. `--sites-bed` instead scores supplied positions regardless of alternate evidence or coverage; it accepts plain BED with one-base, 0-based half-open intervals and no allele information. This allows matched libraries to be compared at identical positions. Input chunk size does not change each SNP's window.

`scripts/validate_haplotype_likelihoods.py` compares matched WGS/TAPS diagnostics with a truth VCF and GATK calls in confident regions. Its input JSON names `reference`, `truth`, `confidence_bed`, `gatk`, `taps`, `wgs` and `region`. Truth alleles are used only for scoring; common focal positions include the discovery union and truth SNP positions. It also evaluates TAPS reads under the ordinary model and compares 16- versus 32-node quadrature. Add `--pair-hmm` for matched realignment arms, and `--include-indel-neighborhoods` to include and separately score positions within `--indel-padding` bases of truth indels. The diagnostic currently leaves the trained forests and `call` output unchanged; integrating its features requires a new feature set and retraining.

### Assessing model limitations

`--candidates-only` writes `chrom`, 1-based `pos` and `ref` without evaluating likelihoods. It cannot be combined with supplied sites or pair output. This separates candidate discovery sensitivity from inference at forced truth positions.

The diagnostic also reports `genotype_posteriors` and `site_genotype_posteriors`: ten unphased base-pair probabilities in the order A/A, A/C, C/C, A/G, C/G, G/G, A/T, C/T, G/T, T/T. These use genetic bases rather than the input VCF's REF/ALT indices. `site_posterior` scores the focal SNP without neighboring SNPs or the nuisance indel.

`disconnected_sites` counts candidates excluded by the fragment-link requirements. `resource_limited_sites` counts reachable candidates excluded by site or haplotype-protection budgets; their sum is `omitted_sites`. `pruned_haplotypes`, `rescued_haplotypes`, `expansion_steps`, `maximum_pruned_mass` and the 1-based `context_start`/`context_end` expose the actual search. `offered_linking_fragments` counts fragments covering the focal site and any offered candidate; `retained_linking_fragments` counts those reaching selected neighbors. These support counts include passing observations that cannot establish a genetic link. A fragment can appear in several focal windows, so summed counts are fragment-to-window incidences, not unique molecules.

`scripts/benchmark_haplotype_likelihoods.py` assesses predefined windows with a shared target panel. Its manifest contains `inputs`, a `regions` list of `region`/`split` objects, `negative_sites_per_region`, `seed`, and optional `known_sites`. It retains all confident truth SNPs, near-indel candidate positions and specified disagreements; other negative candidate positions are sampled by SHA-256. Inverse sampling weights estimate calibration over that target population. Reported TP/FP/FN counts refer to the evaluated panel; separately labelled estimates account for unsampled negative candidates. Genotype scoring excludes indels and compares unphased base pairs at truth SNP sites. The script records per-arm commands, runtime and peak memory, along with context strata, reliability bins, genotype errors and the focal-site ablation.

An optional `arms` list selects custom comparisons. Each entry names `name`, `library` (`wgs` or `taps`), `engine` (`aligned` or `pair-hmm`), and optional `arguments`, such as `["--max-haplotypes", "32"]`. Without it, the original four library/engine arms run. `--panel-source <earlier-suite>` reuses byte-identical targets for matching inputs, regions and sampling choices, recording the discovery provenance. This supports parameter sensitivity checks without changing the evaluated positions.

```bash
python3 scripts/benchmark_haplotype_likelihoods.py \
  --manifest assessment.json --rastair target/release/rastair \
  --output-dir tmp/haplotypes-assessment --tabix tabix
```

Both likelihood engines resolve overlapping measurements before scoring haplotypes. Agreement retains the better base quality, without adding qualities. Disagreement retains the better-quality base and reduces its Phred quality by the competing quality: Q40 versus Q30 becomes Q10, while equal-quality conflicts become Q0 and carry no allele preference. Unknown bases carry no competing evidence. The consensus base contributes one sequencing emission and one conversion draw. Mapping uncertainty uses the selected observation's alignment; a mate that supplies no consensus base adds no mapping event. Mapping quality breaks ties between equal base qualities.

`observations` now counts molecular base observations at the focal position. `read_observations` retains the number of filtered read measurements; `overlapping_fragments` counts fragments contributing two measurements there, and `discordant_fragments` counts their base disagreements. A fragment spanning a deletion may contribute realignment evidence without an aligned focal base, so molecular observations can be fewer than fragments. The TSV model identifiers are `aligned_snp_v3` and `local_pair_hmm_v3`; v2 reports used the bounded three-SNP search and optional methylation models, while v1 also used independent overlapping sequencing measurements. Existing read filters are applied before consensus, and a quality reduced to Q0 remains an uninformative observation.
