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

Since Rastair operates on @TAPS and @5Base data,
it has to take methylation artifacts into account.
In a @CpG site with alleles `C` and `T`,
the `T` on the converted strand may be a methylated `C`,
so it says nothing about which copy the read came from
(same for `G>A`).
Such reads are not counted.
Outside a CpG, cytosines are almost never methylated, so these reads count like any other.

```admonish warning
Phasing only works with the `experimental-seqair` backend.
It can not be combined with `--cpgs-only` or BED output alone.
Blocks do not cross the boundary between two processing segments,
so a block is at most one segment long (`--segment-max-length`, 100 kb by default).
```

## Experimental haplotype likelihoods

The `ml haplotypes` diagnostic calculates SNP likelihoods before final variant or genotype calls. It requires the `experimental-seqair` backend:

```bash
cargo run --features experimental-seqair -- ml haplotypes \
  --fasta-file reference.fa reads.bam --region chr12:20028995-20030086 \
  --output snps.tsv --pair-output pairs.tsv
```

Use `--non-taps` for ordinary WGS reads. For TAPS, `--conversion-rate` defaults to `0.95` and `--false-conversion-rate` to `0.003`. These describe conversion of methylated and unmethylated cytosines respectively; they are separate from base-quality sequencing errors. The lambda spike-in is not used to fit these parameters automatically: incomplete lambda conversion can underestimate conversion efficiency by 2–5 percentage points.

Every focal position is evaluated against all four genetic bases. Candidate haplotypes combine up to three nearby SNP sites, with CpG context calculated from each resulting sequence. Reference and de-novo CpGs, and loss of a reference CpG, therefore receive the corresponding conversion model. Both mates share a chromosome assignment. Overlapping bases form one molecular consensus observation and share a conversion event. There is no depth cap; configured flag, quality and read-end filters still apply.

The default `--likelihood-engine aligned` conditions on the input CIGAR alignment. With the default `--methylation-model integrated`, both engines marginalize a uniform methylation prior independently for each haplotype's scored CpG using Gauss-Legendre quadrature. To keep the integral at most two-dimensional, they reduce neighboring context whenever a haplotype would contain more than one scored CpG; `cpg_context_limited` and `omitted_sites` report this. `alignment_indel_reads` reports alignments with indels or reference skips overlapping the SNP window.

Two alternative models retain multiple CpGs without a joint methylation integral:

- `--methylation-model neutral` gives C and T hypotheses identical emissions on OT at potential CpG C positions, and G and A hypotheses identical emissions on OB at potential CpG G positions. Potential positions are the union across candidate sequences, including de-novo and indel-created CpGs. Both confounded alleles receive the same normalized mixture of sequencing-error emissions. Other alleles and opposite-strand observations remain informative; removing only the C/G hypothesis's evidence would bias the comparison toward T/A. PairHMM applies this rule at template coordinates across alignment paths.
- `--methylation-model shrinkage` first uses neutral likelihoods to obtain genotype and molecular chromosome probabilities. It estimates a methylation rate for each candidate CpG, shared across haplotypes that contain that CpG. A symmetric Beta prior shrinks low-information rates toward 0.5 (`--methylation-prior-strength 4`, corresponding to Beta(2,2)). Base and mapping quality enter a separate one-dimensional rate estimate. Only molecules with pilot CpG-carrier probability at least `--methylation-min-cpg-probability 0.9` can fit the rate, preventing ambiguous genetic T/A observations from being treated as conversions. Each fragment is scored with its own rate contribution removed. Estimates use CIGAR coordinates, including insertion coordinates; ambiguous orientations do not fit rates. Candidate CpGs without usable observations retain the prior mean. These are approximate plug-in scores, not an exact joint marginal likelihood, and they do not model allele-specific methylation.

The three-SNP cap and single nuisance indel remain in all models. `--non-taps` bypasses conversion and rate estimation in every model. The conversion defaults remain 0.95 and 0.003; the spike-in is not used to estimate them. No benchmark truth alleles or WGS genotypes enter TAPS likelihood inference.

### Local PairHMM realignment

Use `--likelihood-engine pair-hmm` to sum over alternative alignments:

```bash
cargo run --release --features experimental-seqair -- ml haplotypes \
  --fasta-file reference.fa reads.bam --region chr12:20028995-20030086 \
  --likelihood-engine pair-hmm --output snps.tsv --pair-output pairs.tsv
```

This scalar forward algorithm has match, insertion and deletion states, following the [GATK PairHMM model](https://github.com/broadinstitute/gatk/blob/master/docs/pair_hmm.tex). `--gap-open` defaults to Phred 45 (`10^-4.5`) for each gap type and `--gap-extend` to `0.1`. Match emissions combine conversion and sequencing error; overlapping mates contribute one consensus base to a single conversion draw. Insertion emissions are uniform over consensus bases, including shared insertions counted once. The alignment sum is normalized by its path mass, conditioning on the profile's length and containment in the local window. Unlike GATK's omitted insertion emission, this gives a normalized observation model: uninformative bases cannot favor a longer haplotype.

The realignment interval spans the selected SNPs and nuisance indel, padded by `--realignment-flank` (12 bases by default). Read sequence is cropped using its CIGAR coordinates; insertion bases are retained, soft clips are excluded, and masked or low-quality interior bases become uninformative observations. This is a local realignment prototype, not assembly of full active regions. Reads spanning a deletion at the target still contribute to the alignment sum.

Alternative methylation models use this same interval for candidate CpG context and rate estimation, including when a nuisance indel extends beyond the SNP window. One additional partner base at either edge can establish CpG context, but observations outside the scoring interval cannot fit a methylation rate.

CIGAR insertions and deletions passing the read filters seed competing genetic sequences, without using called genotypes or truth alleles. Proposals are left-normalized against the local reference. The nearest eligible simple indel is retained as a biallelic nuisance site, with a default maximum length of 20 bases (`--max-indel-length`). Deletions covering a selected SNP and other indel proposals are omitted; `omitted_indels` reports their count. Its genotype prior is `1 - --indel-prior` for reference/reference and `--indel-prior / 2` for reference/indel and indel/indel (default prior `0.0001`). Compatible phase assignments share the joint SNP/indel genotype prior. The SNP probabilities marginalize over the indel's presence and phase, rather than committing to an indel call first. Candidate sequence adjacency determines CpGs after either SNPs or indels create or destroy them.

For concordant overlapping CIGAR paths, mates form one molecular profile with a consensus observation at each shared base. Disjoint mates remain separate profiles assigned to the same chromosome. Discordant overlapping CIGAR paths and reference skips use aligned SNP emissions, with the same base-consensus policy and the candidate's edited CpG context; deleted partner positions are skipped, and indel alignment evidence is not scored. `realigned_fragments` and `aligned_fallback_fragments` make that fallback visible. A distant mate outside the scoring interval contributes no mapping event, sequence likelihood, orientation evidence or reference-skip fallback, keeping results independent of input chunks.

In the integrated model, at most one CpG rate per haplotype is fitted by the numerical integral. Other CpGs in its alignment flanks use a fixed nuisance methylation fraction of `0.5`; their conversion probability still depends on the supplied chemistry rates and each candidate's sequence. The alternative models apply their neutral or estimated-rate treatment in the flanks. Richer indel candidates, allele-specific multi-CpG inference and full-read assembly remain future steps. These bounds and fallback behavior need to be considered when interpreting diagnostic probabilities.

`snps.tsv` contains alternate-presence, heterozygous and homozygous-alternate posteriors, expected allele fraction, genotype entropy and natural-log Bayes factors. The window and site-only Bayes factors compare the expanded context, including the nuisance indel, with the focal-site calculation. `pairs.tsv` contains every candidate diploid haplotype pair's natural-log likelihood, prior and posterior. Its `haplotype_a/b` strings give SNP bases at the listed 1-based positions; `indel_a/b` give an unanchored `0-based-start:deleted>inserted` replacement, or `.` for reference. The SNP prior assigns `1 - --snp-prior` to the reference genotype and distributes the remaining mass uniformly across the other genotypes at each site; compatible phase assignments share that mass equally. These model probabilities are not empirically calibrated ML call scores.

By default, focal candidates come from observed alternate bases and potential de-novo CpG alleles whose converted base resembles the reference. Potential CpG partners include the reference, neighboring observed SNP alleles and, for PairHMM, indel proposals. Thus an indel-created CpG can seed a C/G hypothesis even when conversion leaves only reference-looking T/A observations. `--sites-bed` instead scores supplied positions regardless of alternate evidence or coverage; it accepts plain BED with one-base, 0-based half-open intervals and no allele information. This allows matched libraries to be compared at identical positions. Input chunk size does not change each SNP's window.

`scripts/validate_haplotype_likelihoods.py` compares matched WGS/TAPS diagnostics with a truth VCF and GATK calls in confident regions. Its input JSON names `reference`, `truth`, `confidence_bed`, `gatk`, `taps`, `wgs` and `region`. Truth alleles are used only for scoring; common focal positions include the discovery union and truth SNP positions. It also evaluates TAPS reads under the ordinary model and compares 16- versus 32-node quadrature. Add `--pair-hmm` for matched realignment arms, and `--include-indel-neighborhoods` to include and separately score positions within `--indel-padding` bases of truth indels. The diagnostic currently leaves the trained forests and `call` output unchanged; integrating its features requires a new feature set and retraining.

### Assessing model limitations

`--candidates-only` writes `chrom`, 1-based `pos` and `ref` without evaluating likelihoods. It cannot be combined with supplied sites or pair output. This separates candidate discovery sensitivity from inference at forced truth positions.

The diagnostic also reports `genotype_posteriors` and `site_genotype_posteriors`: ten unphased base-pair probabilities in the order A/A, A/C, C/C, A/G, C/G, G/G, A/T, C/T, G/T, T/T. These use genetic bases rather than the input VCF's REF/ALT indices. `site_posterior` scores the focal SNP without neighboring SNPs or the nuisance indel.

`cpg_omitted_sites` counts neighboring SNP candidates removed by the CpG restriction after the ordinary site cap; `cpg_omitted_indels` counts eligible nuisance indels removed for that reason. `offered_linking_fragments` counts fragments covering the focal site and any neighbor retained by the site cap, and `retained_linking_fragments` counts those still reaching a neighbor after CpG reduction. A fragment can appear in several focal windows, so summing these fields counts fragment-to-window incidences, not unique molecules.

`scripts/benchmark_haplotype_likelihoods.py` assesses predefined windows with a shared target panel. Its manifest contains `inputs`, a `regions` list of `region`/`split` objects, `negative_sites_per_region`, `seed`, and optional `known_sites`. It retains all confident truth SNPs, near-indel candidate positions and specified disagreements; other negative candidate positions are sampled by SHA-256. Inverse sampling weights estimate calibration over that target population. Reported TP/FP/FN counts refer to the evaluated panel; separately labelled estimates account for unsampled negative candidates. Genotype scoring excludes indels and compares unphased base pairs at truth SNP sites. The script records per-arm commands, runtime and peak memory, along with context strata, reliability bins, genotype errors and the focal-site ablation.

An optional `arms` list selects custom comparisons. Each entry names `name`, `library` (`wgs` or `taps`), `engine` (`aligned` or `pair-hmm`), and optional `arguments`, such as `["--methylation-model", "shrinkage"]`. Without it, the original four library/engine arms run. `--panel-source <earlier-suite>` reuses byte-identical targets for matching inputs, regions and sampling choices, recording the discovery provenance. This supports parameter sensitivity checks without changing the evaluated positions.

```bash
python3 scripts/benchmark_haplotype_likelihoods.py \
  --manifest assessment.json --rastair target/release/rastair \
  --output-dir tmp/haplotypes-assessment --tabix tabix
```

Both likelihood engines resolve overlapping measurements before scoring haplotypes. Agreement retains the better base quality, without adding qualities. Disagreement retains the better-quality base and reduces its Phred quality by the competing quality: Q40 versus Q30 becomes Q10, while equal-quality conflicts become Q0 and carry no allele preference. Unknown bases carry no competing evidence. The consensus base contributes one sequencing emission and one conversion draw. Mapping uncertainty uses the selected observation's alignment; a mate that supplies no consensus base adds no mapping event. Mapping quality breaks ties between equal base qualities.

`observations` now counts molecular base observations at the focal position. `read_observations` retains the number of filtered read measurements; `overlapping_fragments` counts fragments contributing two measurements there, and `discordant_fragments` counts their base disagreements. A fragment spanning a deletion may contribute realignment evidence without an aligned focal base, so molecular observations can be fewer than fragments. The TSV model identifiers are `aligned_snp_v2` and `local_pair_hmm_v2`; earlier v1 reports used independent overlapping sequencing measurements. Existing read filters are applied before consensus, and a quality reduced to Q0 remains an uninformative observation.

The restriction to one fitted CpG can also discard useful genetic linkage. A test-only two-CpG oracle integrates independent rates exactly on a small constructed case, retaining evidence against a singleton error that the bounded model favors. The oracle applies the same overlap consensus as the production likelihood calculation. This establishes a possible failure mode, not the accuracy of a future multi-CpG implementation. Integrated-model flank CpGs still use fixed methylation 0.5; the alternatives retain multiple CpGs without this context restriction. The ordinary site cap remains three SNPs.
