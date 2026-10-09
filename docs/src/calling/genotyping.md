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
Each block opens with `0|1`; nothing relates one block to another, and heterozygous sites no read pair links stay `0/1` without `PS`.

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
