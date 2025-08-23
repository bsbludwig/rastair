# Examples

rastair has two main modes:
1. Call methylation per genomic position
2. Annotate individual reads with the methylation states of the CpGs they contain

In addition, there are a number of convenience commands, e.g. to convert between different output file formats, and a small set of utility scripts written in R to produce quality-control metrics.

Below we will give some examples that cover common use-cases. You can find an in-depth documentation of the complete command line syntax [here](cli.md).

## 0. Parameters shared between various sub-commands
Nearly all rastair sub-commands are capable of multi-threading. If your CPU has multiple cores, you can e.g. use `-@ 4` to parallelise operations across 4 cores.

Similarly, you can restrict processing to only one genomic interval using the `-l` parameter. To only process a specific chromosome, you can do `-l chr19`. You can also chose an interval within the chromosome: `-l chr19:6103156-6143156`.

## 1. Call genomic variants and methylated positions

By default, rastair will use a built-in machine-learning model to classify variants as real or false. The output should be a `.vcf.gz` or a `.bcf` file:

```bash
rastair2 call -r reference.fa.gz -o test.bcf test.bam
```

The resulting @VCF output still contains candidate variants that did not pass all filters. If you just want to store high-confidence calls, you can write uncompressed vcf output to `STDOUT` and filter it on the fly with e.g. [bcftools](https://samtools.github.io/bcftools/bcftools.html):

```bash
rastair2 call -r reference.fa.gz --vcf - test.bam | bcftools view -f PASS -o test.bcf
```

You can find a description of all custom VCF fields used by rastair [here](formats/vcf-fields.md).

## 2. Only call methylation, do not report genetic variants (very fast)

In cases where all you need is a table of methylation counts in genomic - and putatively _de-novo_ - CpG positions, you can use the `--cpgs-only` parameter and use `--bed` output, which will greatly speed up the processing compared to calling all putative genetic variants:

```bash
rastair2 call -r reference.fa.gz --cpgs-only --bed test.bed.gz test.bam
```

The reference for the meaning of the different columns in the bed output format can be found [here](formats/bed.md). Rastair will automatically produce [bgzip compressed](https://www.htslib.org/doc/bgzip.html) files if the output file name ends in `.gz`. You can then [index these with tabix](https://www.htslib.org/doc/tabix.html) for rapid access to specific genomic ranges:

```bash
tabix -p bed test.bed.gz
# Fetch calls in a random genomic region:
tabix test.bed.gz chr19:6103156-6143156
```

In some cases, you might prefer to write the bed output to `STDOUT` and pipe it into another unix tool, e.g. to only report positions that are CpG in the references (ie exclude @denovo):

```bash
rastair2 call -r reference.fa.gz --cpgs-only --bed - test.bam | grep -Fw REF
```

Sometimes, it might be desirable to ignore a certain number of bases at the beginning or end of a read when counting methylated positions. This can e.g. account for loss of methylation due to sonication damage. The command-line argument for this was inspired by [MethylDackel](https://github.com/dpryan79/MethylDackel). However, we decided to only have one set of parameters: `--nOT` and `--nOB`. Each of them takes a comma-separated list of 4 integers: `[r1_5',r1_3',r2_5',r2_3']`, denoting the number of bases from the start/end of read 1/2 that should be ignored. **Unlike MethylDackel, rastair's command line arguments always refer to the start/end position of the read in 5' -> 3' direction, not the position in the reference after alignment**. To give an example: imagine the following read pair

```text
    000000000111111111122
    123456789012345678901
R1: CG--------TG--------->
R2:                         <AC-------------GC-
                             876543210987654321
                             111111111000000000
```

This "F1R2" read represents the OT (ie R1 is the OT, and R2 is the reverse complement of the OT). A parameter of `--nOT 0,5,0,5` will exclude the `A` at position 18 in R2, because it ocurrs within 5 bases from the end of R2 _in read coordinates_, not in reference coordinates.

## 3. Report methylation per-read in bed format

