# Examples

rastair has two main modes:
1. Call @methylation per genomic position
2. Annotate individual @read:pl with the methylation states of the @CpG:pl they contain

In addition, there are a number of convenience commands, e.g. to convert between different output file formats, and a small set of utility scripts written in R to produce quality-control metrics.

Below we will give some examples that cover common use-cases. You can find an in-depth documentation of the complete command line syntax [here](cli.md).

```admonish tip
Nearly all rastair sub-commands are capable of multi-threading.
You can e.g. use `-@ 4` to parallelise operations across 4 cores.
By default, rastair will try to use all available cores on your system.
```

## 1. Call genomic variants and methylated positions

By default, rastair will use a built-in @ML model to classify @variant:pl as real or false.

If writing to a @VCF (give a `.vcf.gz` file name) or @BCF file (give a filename ending in `.bcf`),
all high-confidence calls will be included:

```bash
rastair call -r reference.fa.gz -o test.vcf.gz test.bam
```

```admonish tip
You can restrict processing to only a specific genomic region. For example `-l chr19` will process the entire chromosome 19.
You can also chose an interval within the chromosome: `-l chr19:6103156-6143156`.
```

Rastair will only report variants that pass its internal filters. If you want to report absolutely all candidate variants, regardless of their score, you can use the `--all` flag:

```bash
rastair call -r reference.fa.gz --all -o test.vcf.gz test.bam
```
```admonish warning
Beware that vcf files generated with `--all` can get very large!
```

Rastair will annotate variants with only a subset of all possible `INFO` and `FORMAT` tags. You can choose a custom set of annotations with e.g.

```bash
rastair call -r reference.fa.gz --vcf-info-fields AD,BQ,MQ,MQ0,ENT100,SC5,CPG,CPGnovo --vcf-format-fields GT,M5mC,ML -o test.bcf test.bam
```

You can find a description of all custom VCF fields used by rastair [here](formats/vcf-fields.md).

```admonish tip
rastair can produce uncompressed VCF, bgzip-compressed `vcf.gz` and `bcf` files.
The format will be guessed from the file ending.
If no output name is given, rastair writes VCF to @STDOUT.
```

## 2. Only call methylation, do not report genetic variants (very fast)

@BED output will only report CpG -- and de-novo CpG -- positions, which is a lot faster than processing all potential variants:

```bash
rastair call -r reference.fa.gz --bed test.bed.gz test.bam
```

```admonish tip
Rastair will automatically produce @bgzip compressed files if the output file name ends in `.gz`. They are also automaticall indexed with @tabix for rapid access to specific genomic ranges:

> ```bash
> tabix test.bed.gz chr19:6103156-6143156
> ```

This is usually orders of magnitude faster than using e.g. @bedtools. You can also stream the @BED file to @STDOUT by setting `--bed -`.
```

In some cases, you might prefer to write the bed output to @STDOUT and pipe it into another unix tool, e.g. to only report positions that are CpG in the references (ie exclude @denovo):

```bash
rastair call -r reference.fa.gz --cpgs-only --bed - test.bam | grep -Fw REF
```

### Ignore positions at the edges of reads

Sometimes, it might be desirable to ignore a certain number of bases at the beginning or end of a read when counting methylated positions.
This can e.g. account for loss of methylation due to sonication damage.
The command-line argument for this was inspired by [MethylDackel](https://github.com/dpryan79/MethylDackel).
However, we decided to only have one set of parameters: `--nOT` and `--nOB`.
Each of them takes a comma-separated list of 4 integers: `[r1_5',r1_3',r2_5',r2_3']`, denoting the number of bases from the start/end of read 1/2 that should be ignored.

```admonish warning
Unlike MethylDackel, **rastair's command line arguments always refer to the start/end position of the read in 5'&rarr;3' direction**, not the position in the reference after alignment!
```

To give an example: imagine the following read pair

```text
    000000000111111111122
    123456789012345678901
R1: CG--------TG--------->
R2:                         <AC-------------GC-
                             876543210987654321
                             111111111000000000
```

This "[F1R2](https://gatk.broadinstitute.org/hc/en-us/community/posts/360075017111-strand-bias-and-orientation-bias)" read represents the @OT (ie R1 is the OT, and R2 is the reverse complement of the OT).
A parameter of `--nOT 0,5,0,5` will exclude the `A` at position 18 in R2, because it ocurrs within 5 bases from the end of R2 _in read coordinates_, not in reference coordinates.

## 3. Report methylation per-read

Many applications require information about the methylation of individual bases in each sequenced read. Rastair can report this information in bed format:

```bash
rastair per-read -r reference.fa.gz --bed test_per-read.bed.gz test.bam
```

For a description of the per-read bed format, see the [BED format](formats/bed.md#per-read-methylation) section.

```admonish info
The per-read output is also the input for your QC pipeline. For more details, read the [qc section](src/qc.md) of this manual.
```
