# Genotyping

Rastair uses @ML scores to determine which @altAllele:pl represent true genetic @variant:pl.
Their zygosity comes from the pair-HMM that also calls indels: it scores the reads against the local haplotypes, TAPS conversion included, and picks the likeliest pair of haplotypes that carries the called alleles.
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
- Heterozygous (`0/1`): the likeliest pair carries one copy of the alternate allele
- Homozygous Alternate (`1/1`): the likeliest pair carries two copies
- Compound Heterozygous (`1/2`): the likeliest pair carries two alternate alleles, one on each copy
- Missing (`./.`): the alternate allele passed filters but the pair-HMM could not score it, at an `N` in the reference

## Confidence Scoring

`GC` is the Phred-scaled chance that the genotype is wrong:

- For `0/0` calls with no passing alternate: a binomial test of the reads that show a SNP against the error rate
- For variant calls: the posterior of the other copy count over the pairs of haplotypes consistent with the called alleles; `GL` carries the same value

With `--no-ml` the pair-HMM also decides which alleles are present, and a passing alternate it finds absent is written `0/0`.

## Biological Interpretation

A genotype represents the complete set of observed alleles in an individual.
All alternate alleles listed in a VCF record's ALT field are simultaneously present at that position in at least one chromosome copy.
In most cases:

- If only one variant passes filters, the genotype indicates whether it's present on one chromosome copy (heterozygous `0/1`) or both copies (homozygous alternate `1/1`)
- If multiple variants pass filters, they represent compound heterozygosity where different variants are present on different chromosome copies (e.g., `1/2`)

## Strand-Specific Counting

At a @CpG, a `T` over a `C` on the converted strand may be a methylated `C` (`A` over `G` likewise).
The pair-HMM scores such reads against both alleles at the strand's methylation level, so they count only as far as they tell the alleles apart.
For all other variant types, both strands contribute to genotyping.

## Phasing

Rastair groups heterozygous variants that share a read pair
into phase blocks, written with a `|` separator and a `PS` field:

```
chr12  20028875  .  A  G  ...  GT:PS  0|1:20028875
chr12  20028894  .  C  T  ...  GT:PS  1|0:20028875
```

Both records carry `PS=20028875`, the position of the block's first heterozygous site, so they belong to one block.
Within a block, the alleles written first all sit on the same chromosome copy: above, the reference `A` and the alternate `T` are on one copy, the alternate `G` and the reference `C` on the other.
Each block opens with `0|1`; nothing relates one block to another, and heterozygous sites no read pair links stay `0/1` without `PS`.

Each read votes with the same pair-HMM that called the genotype:
for the allele whose local haplotype explains the read better, weighted by how much better.
A read needs a clear preference (ten times likelier under one allele) to vote,
and a read neither haplotype explains (one with a soft clip or an indel the pair-HMM does not hold) does not vote.

Since Rastair operates on @TAPS and @5Base data,
it has to take methylation artifacts into account.
In a @CpG site with alleles `C` and `T`,
the `T` on the converted strand may be a methylated `C`.
The pair-HMM scores such a read as about equally likely under both alleles,
so it casts no vote (same for `G>A`).

Heterozygous indels called by the pair-HMM (see [indels](indels.md)) take part the same way.
An indel only becomes part of a block when the pair-HMM is confident of its heterozygous genotype: a doubtful one would mostly join two blocks the wrong way round.

```admonish warning
Phasing is not done with `--cpgs-only` or BED output alone, which drop the variants it links.
Blocks do not cross the boundary between two processing segments,
so a block is at most one segment long (`--segment-max-length`, 100 kb by default).
```
