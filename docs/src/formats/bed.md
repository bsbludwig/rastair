# BED Format

Rastair2 can output @BED files of two different kinds:

1. **CpG sites**: A file containing all @CpG sites with their @methylation status.
2. **Per-read methylation**: A file containing the methylation status of each CpG site for each @read.

## CpG Sites

The BED file for CpG sites contains the following columns:

| Column          | Description                                                                     |
| --------------- | ------------------------------------------------------------------------------- |
| `chrom`         | Chromosome name                                                                 |
| `start`         | Start position of the CpG site (0-based)                                        |
| `end`           | End position of the CpG site (1-based)                                          |
| `name`          | Name of the CpG site (e.g., "CpG1")                                             |
| `beta_est`      | Estimated beta value for methylation                                            |
| `strand`        | Strand information (e.g., "+", "-")                                             |
| `unmod`         | Number of unmethylated reads                                                    |
| `mod`           | Number of methylated reads                                                      |
| `no_snp`        | Number of reads not counting as @SNP:pl                                         |
| `snp`           | Number of reads counting as @SNP:pl                                             |
| `coverage`      | Total coverage at the CpG site                                                  |
| `genotype`      | `C/C`, `C/T`, `G/G`, `G/A`, `T/T`, or `A/A`                                     |
| `gt_p_score`    | P-value for the genotype call                                                   |
| `gt_conf_score` | Confidence score for the genotype call                                          |
| `cpg`           | `REF` if CpG site occurs in @referenceGenome, `NEW` if it is a @denovo site |

## Per-Read Methylation

t.b.d.
