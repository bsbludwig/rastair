# Measuring accuracy

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

## Comparing binaries

Score every arm of a comparison with `--gpu` or every arm without it: the two are not bitwise
identical. A binary built before seqair became the only pileup backend read pileups with htslib
unless built with `--features experimental-seqair` (its log then said "Using experimental seqair
backend" once); the backends differ in mate dedup and soft-clip handling, so do not compare arms
built on each.
