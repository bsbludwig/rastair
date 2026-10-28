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
| `chr`           | Chromosome name (the header line starts with `#chr`)                        |
| `start`         | Start position of the CpG site (0-based)                                    |
| `end`           | End position of the CpG site (1-based)                                      |
| `name`          | Always `.`                                                                  |
| `beta_est`      | Estimated beta value for methylation (empty if no read shows the CpG)       |
| `strand`        | The strand `beta_est` is read on: `+` (OT) or `-` (OB), `.` if none        |
| `unmod`         | Number of unmethylated reads, the ones `beta_est` counts                    |
| `mod`           | Number of methylated reads, the ones `beta_est` counts                      |
| `no_snp`        | Number of reads not counting as @SNP:pl                                     |
| `snp`           | Number of reads counting as @SNP:pl                                         |
| `coverage`      | Total coverage at the CpG site                                              |
| `genotype`      | The genotype as two bases, e.g. `C/C`, `C/T` or `T/T`; the reference twice where the VCF has `./.` |
| `gt_p_score`    | The same value as `gt_conf_score`                                           |
| `gt_conf_score` | The genotype's `GQ`: Phred-scaled chance that it is wrong; 0 where the VCF has no `GQ` |
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
| `num_cpg`       | Number of CpGs in a read: the reference CpGs it covers and the @denovo ones it shows (in `mod_denovos` or `unmod_denovos`) |
| `num_mod`       | Number of modified CpGs: `mod_cpgs` and `mod_denovos` |
| `mod_cpgs`      | Positions in read of modified CpGs            |
| `unmod_cpgs`    | Positions in read of unmodified CpGs          |
| `snp_cpgs`      | Positions in read that are @SNP:pl (mutated)  |
| `mod_denovos`   | Positions in read of @denovo sites that are methylated (needs `--calls`) |
| `unmod_denovos` | Positions in read of @denovo sites that are unmethylated (needs `--calls`) |

Each read's own base decides whether a CpG is modified, at reference and @denovo CpGs alike,
not the site's β, as in the tags [`bam`](../cli.md#rastair-bam) writes:
an OT read shows a methylated C as `T` and an unmethylated one as `C`,
an OB read a methylated G as `A` and an unmethylated one as `G`.
A read showing another base is listed in `snp_cpgs` at a reference CpG and in neither column at a @denovo one.

**Note:** The positions in reads take @indel:pl into account,
meaning that the positions are relative to the read, not the reference genome.
If `--count-clipped` is set, it will also include leading clipped bases.
