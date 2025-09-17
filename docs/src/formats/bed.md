# BED Format

Rastair can output @BED files of two different kinds:

1. **CpG sites**:
   A file containing all @CpG sites with their @methylation status.
   Generated using the [`call`] command and specifying `--bed`
   (or using `convert`).
2. **Per-read methylation**:
   A file containing the methylation status of each CpG site for each @read.
   Generated using the [`per-read`] command.

[`call`]: ../cli.md#rastair-call
[`per-read`]: ../cli.md#rastair-per-read

## CpG Sites

The BED file for CpG sites contains the following columns:

| Column          | Description                                                                 |
| --------------- | --------------------------------------------------------------------------- |
| `chrom`         | Chromosome name                                                             |
| `start`         | Start position of the CpG site (0-based)                                    |
| `end`           | End position of the CpG site (1-based)                                      |
| `name`          | Name of the CpG site (e.g., "CpG1")                                         |
| `beta_est`      | Estimated beta value for methylation (empty string if not present)          |
| `strand`        | Strand information (e.g., "+", "-")                                         |
| `unmod`         | Number of unmethylated reads                                                |
| `mod`           | Number of methylated reads                                                  |
| `no_snp`        | Number of reads not counting as @SNP:pl                                     |
| `snp`           | Number of reads counting as @SNP:pl                                         |
| `coverage`      | Total coverage at the CpG site                                              |
| `genotype`      | `C/C`, `C/T`, `G/G`, `G/A`, `T/T`, or `A/A`                                 |
| `gt_p_score`    | P-value for the genotype call                                               |
| `gt_conf_score` | Confidence score for the genotype call                                      |
| `cpg`           | `REF` if CpG site occurs in @referenceGenome, `NEW` if it is a @denovo site |

## Per-Read Methylation

The BED file for per-read methylation contains the following columns:

| Column          | Description                                   |
| --------------- | --------------------------------------------- |
| `chr`           | Chromosome name                               |
| `start`         | Start position                                |
| `end`           | End position                                  |
| `read_id`       | Name of read                                  |
| `mapq`          | Mapq of read                                  |
| `orientation`   | Orientation of read, either `+` or `-`        |
| `insert_size`   | Absolute fragment length (non-directional)    |
| `read_length`   | Read length                                   |
| `flag`          | Flag of read (decimal, same as in @BAM)       |
| `num_cpg`       | Number of CpGs in a read                      |
| `num_mod`       | Number of modified CpGs                       |
| `mod_cpgs`      | Positions in read of modified CpGs            |
| `unmod_cpgs`    | Positions in read of unmodified CpGs          |
| `snp_cpgs`      | Positions in read that are @SNP:pl (mutated)  |
| `mod_denovos`   | Positions in read of @denovo that are mutated |
| `unmod_denovos` | Positions in read of @denovo that are mutated |

**Note:** The positions in reads take @indel:pl into account,
meaning that the positions are relative to the read, not the reference genome.
If `--count-clipped` is set, it will also include leading clipped bases.
