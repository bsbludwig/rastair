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

Every focal position is evaluated against all four genetic bases. Candidate haplotypes combine up to three nearby SNP sites, with CpG context calculated from each resulting sequence. Reference and de-novo CpGs, and loss of a reference CpG, therefore receive the corresponding conversion model. Both mates share a chromosome assignment, and overlapping observations share a conversion event. There is no depth cap; configured flag, quality and read-end filters still apply.

The first prototype conditions on the input CIGAR alignment. It marginalizes a uniform methylation prior independently for each haplotype's CpG using Gauss-Legendre quadrature. To keep the integral at most two-dimensional, it reduces neighboring context whenever a haplotype would contain more than one scored CpG; `cpg_context_limited` and `omitted_sites` report this. `alignment_indel_reads` reports alignments with indels or reference skips overlapping the SNP window. These limitations matter when interpreting the probabilities: the prototype does not perform PairHMM realignment or infer multiple CpG rates per haplotype window.

`snps.tsv` contains alternate-presence, heterozygous and homozygous-alternate posteriors, expected allele fraction, genotype entropy and natural-log Bayes factors. The window and site-only Bayes factors distinguish neighboring evidence from the focal-site calculation. `pairs.tsv` contains every candidate diploid haplotype pair's natural-log likelihood, prior and posterior. The prior assigns `1 - --snp-prior` to the reference genotype and distributes the remaining mass uniformly across the other genotypes at each site; compatible phase assignments share that mass equally. These model probabilities are not empirically calibrated ML call scores.

By default, focal candidates come from observed alternate bases and potential de-novo CpG alleles whose converted base resembles the reference. `--sites-bed` instead scores supplied positions regardless of alternate evidence or coverage; it accepts plain BED with one-base, 0-based half-open intervals and no allele information. This allows matched libraries to be compared at identical positions. Input chunk size does not change each SNP's window.

`scripts/validate_haplotype_likelihoods.py` compares matched WGS/TAPS diagnostics with a truth VCF and GATK calls in confident regions. Its input JSON names `reference`, `truth`, `confidence_bed`, `gatk`, `taps`, `wgs` and `region`. Truth alleles are used only for scoring; common focal positions include the discovery union and truth SNP positions. It also evaluates TAPS reads under the ordinary model and compares 16- versus 32-node quadrature. The diagnostic currently leaves the trained forests and `call` output unchanged; integrating its features requires a new feature set and retraining.
