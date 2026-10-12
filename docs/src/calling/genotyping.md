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

`--phase` can be used to group heterozygous variants that share a read pair
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

Heterozygous indels called by the pair-HMM (`--experimental-indels`, see [indels](indels.md)) take part too.
At an indel a read has no single base to vote with, so it votes for the allele whose local haplotype explains it better, weighted by how much better.
An indel only becomes part of a block when the pair-HMM is confident of its heterozygous genotype: a doubtful one would mostly join two blocks the wrong way round.

### Methylation by haplotype

`--phase` also splits each CpG's reads by the haplotype of the read pair they came from:

```
chr12  20028930  .  C  .  ...  GT:M5mC:HPS:HPC:HPX:HPF  0/0:0.48:20028875:9,0,1,11:1,1:31,2,3,29
```

`HPC` counts the reads on the strand where TAPS conversion shows (the top strand at a `C`, the bottom strand at a `G`): unconverted then converted on haplotype 1, then the same on haplotype 2.
`HPS` names the phase set whose haplotypes these are: a `PS` from the same file, or the position of a nearby heterozygous SNV or indel that could not be linked to another, in which case haplotype 1 carries the first allele of that site's `GT`.
Where reads reach several sets, the one that placed the most of them is reported.
A column is never split by its own alleles, only by the heterozygous sites around it.
Above, haplotype 2 is methylated and haplotype 1 is not: allele-specific methylation.

`HPX` says whether each haplotype has a CpG here at all.
At a heterozygous `C>T` the chromosome carrying the `T` has no cytosine, and every one of its reads converts by construction; `HPX` marks that haplotype `0`, so its share of `HPC` is read as the genotype rather than as methylation.
Where the column cannot tell, `HPX` is omitted.

At such a CpG, `M5mC` without `--phase` assumes half the converted reads are the variant's.
With `--phase`, where `HPX` is known and the haplotype carrying the cytosine has reads, `M5mC` is that haplotype's methylation instead.
The CpG's other base needs no such correction: its `M5mC` only counts reads that also show the cytosine next to it, which the variant's chromosome never does.

`HPF` counts the same evidence per molecule: each fragment votes once for the phase set, by the majority of its reads at the set's CpGs, counting only reads on a haplotype that carries the CpG at that column.
A molecule crossing several CpGs is one observation, not several, so `HPF` is what a statistical test of allele-specific methylation should use; it is the same on every record of a set.

```admonish warning
Phasing only works with the `experimental-seqair` backend.
It can not be combined with `--cpgs-only` or BED output alone.
Blocks do not cross the boundary between two processing segments,
so a block is at most one segment long (`--segment-max-length`, 100 kb by default).
```
