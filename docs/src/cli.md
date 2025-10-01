# Command-Line Help for `rastair`

This document contains the help content for the `rastair` command-line program.

**Command Overview:**

* [`rastair`↴](#rastair)
* [`rastair call`↴](#rastair-call)
* [`rastair per-read`↴](#rastair-per-read)
* [`rastair bam`↴](#rastair-bam)
* [`rastair convert`↴](#rastair-convert)
* [`rastair view`↴](#rastair-view)
* [`rastair mbias`↴](#rastair-mbias)

## `rastair`

Rastair -- detect genetic variants and methylated positions from short-read sequencing data created using TET-Assisted Pyridine-Borane Sequencing.

See <https://docs.rastair.com/> for more information.

**Usage:** `rastair [OPTIONS] <COMMAND>`

###### **Subcommands:**

* `call` — Call methylated positions
* `per-read` — Call methylation per-read
* `bam` — Add methylation information to BAM files
* `convert` — Convert between different file formats
* `view` — View internal format as JSON lines
* `mbias` — Calculate conversion per base position in read

###### **Options:**

* `-v`, `--verbose` — Enable more logging

   You can also use the `RASTAIR_LOG` environment variable to configure logging in a more precise way. See the documentation of the `tracing-subscriber` library to learn more.



## `rastair call`

Call methylated positions

Process TAPS-sequenced BAM files and call methylated positions.

If no output file is specified, the output is written to stdout. You can use `--vcf` and `--bed` to write to files instead.

If using `-c` (`--cpgs-only`), all CpG positions in the reference as well as de-novo CpGs are written. Stdout will default to BED.

Only variants that pass all filters are written by default. Use `--all` to get a full VCF file.

**Usage:** `rastair call [OPTIONS] --fasta-file <FASTA_FILE> <BAM_FILE>`

###### **Arguments:**

* `<BAM_FILE>` — Path to sorted and indexed BAM file

###### **Filter Options:**

* `--keep-overlapping-reads` — Whether to keep overlapping reads

  Default value: `false`
* `-c`, `--cpgs-only` — Report CpGs only and default to BED output

   Only report positions that are CpGs in the reference or variants that would result in a de-novo CpG.

  Default value: `false`
* `-q`, `--min-mapq <MIN_MAPQ>` — Minimum mapping quality to consider a read

  Default value: `1`
* `-Q`, `--min-baseq <MIN_BASEQ>` — Minimum base quality to consider a base

  Default value: `10`
* `--nOT <N_OT>` — For OT reads, exclude `[r1_start, r1_end, r2_start, r2_end]` bases from counting.

   The coordinates are relative to the read, so start is the distance from the 5' of the read, the end is the distance to the 3', irrespective of which way around the read aligns to the reference.

   Also note that the distance is relative to read length, not alignment length, so soft-clipped bases count, too!

  Default value: `0,0,0,0`
* `--nOB <N_OB>` — For OB reads, exclude `[r1_start, r1_end, r2_start, r2_end]` bases from counting.

   The coordinates are relative to the read, so start is the distance from the 5' of the read, the end is the distance to the 3', irrespective of which way around the read aligns to the reference.

   Also note that the distance is relative to read length, not alignment length, so soft-clipped bases count, too!

  Default value: `0,0,0,0`
* `-f`, `--include-flags <INCLUDE_FLAGS>` — Include reads that match all of these bit-flags

  Default value: `3`
* `-F`, `--exclude-flags <EXCLUDE_FLAGS>` — Exclude reads that match any of these bit-flags

  Default value: `3852`
* `--cpg-novo-min-depth <CPG_NOVO_MIN_DEPTH>` — Minimum reads needed in support of de-novo CpG

  Default value: `2`
* `--cpg-novo-min-baseq <CPG_NOVO_MIN_BASEQ>` — Minimum base quality for de-novo CpGs

  Default value: `15`
* `--cpg-novo-min-mapq <CPG_NOVO_MIN_MAPQ>` — Minimum mapping quality for de-novo CpGs

  Default value: `50`
* `--cpg-novo-min-vaf <CPG_NOVO_MIN_VAF>` — Minimum variant allele frequency for de-novo CpGs

  Default value: `0.2`
* `--m-vaf-min <M_VAF_MIN>` — The minimum variant allele frequency

  Default value: `0.2`
* `--m-min-depth <M_MIN_DEPTH>` — The minimum number of reads to call a position as methylated

  Default value: `3`
* `--m-min-denovo-depth <M_MIN_DENOVO_DEPTH>` — The minimum number of reads required as evidence for a de novo CpG

  Default value: `2`
* `--m-bq-ratio-min <M_BQ_RATIO_MIN>` — The minimum quality ratio `(ad_alt*bq_alt + 1) / (ad_ref*bq_ref + 1)`

  Default value: `0.27`
* `--m-read-position-min <M_READ_POSITION_MIN>` — The minimum relative position in read for alt allele evidence

  Default value: `0.2`
* `--m-read-position-max <M_READ_POSITION_MAX>` — The maximum relative position in read for alt allele evidence

  Default value: `0.8`
* `--m-max-coverage <M_MAX_COVERAGE>` — The maximum coverage depth for methylation calling

  Default value: `1000`
* `--thresholds` — Only use hard thresholds to call variants and methylation events.

   This disables using the machine learning models. This will make rastair much faster, but at the cost of accuracy.
* `--ml <ML>` — Use machine learning model with this threshold value to call variants and methylation events

   When specified, a ML model will classify positions with a prediction score. Anything above this threshold is considered PASS.

   For consistency with `--thresholds`, this option can be also be specified as `--ml` without a value, which will use the default threshold.

  Default value: `0.80`
* `--model-cpg <MODEL_CPG>` — Path to the model for CpG positions

   Default is the bundled model in the Rastair binary.
* `--model-denovo-cpg <MODEL_DENOVO_CPG>` — Path to the model for de novo CpG positions

   Default is the bundled model in the Rastair binary.
* `--model-others <MODEL_OTHERS>` — Path to the model for other positions

   Default is the bundled model in the Rastair binary.
* `--bed-include-empty` — Include CpG positions with zero coverage

   This can be useful to get a complete list of CpG positions in the output BED file. Note that this requires the input data to contain a complete list of CpG positions, e.g. by using the `--cpgs-only` option when calling methylation.

###### **Input Options:**

* `-r`, `--fasta-file <FASTA_FILE>` — Path to sorted and indexed (via samtools faidx) FASTA file. Can be bgzip compressed, but requires both a gzi index and a fai index
* `-l`, `--region <REGION>` — Restrict to a specific chromosome or region of a chromosome. Format is "chr", "chr:start" or "chr:start-end", where start is 1-based and end is inclusive

###### **Output Options:**

* `-o`, `--vcf <VCF>` — VCF/BCF output file path (use - to write to stdout)

   Format is guessed based on the file extension: `.vcf` for VCF (uncompressed), `.vcf.gz` for VCF (compressed), `.bcf` for BCF (compressed) `.mpk.lz4` for internal format (Message Pack, LZ4-compressed)
* `--all` — Output all positions, even if they are not CpG or de-novo CpG candidates or do not pass filters
* `--bed <BED>` — Output BED file with the called methylated positions
* `--bed-format <BED_FORMAT>` — Format of the output BED file

   If not specified, the format is guessed based on the file extension.

  Possible values:
  - `bed-gz`:
    BGZIP compressed file, usually `.bed.gz`
  - `bed`:
    Regular BED file, usually `.bed`


###### **Processing Options:**

* `--segment-max-length <SEGMENT_MAX_LENGTH>` — Maximum length of a segment in bases

   Used for splitting work between threads. Tweak this to adjust memory usage.

  Default value: `100000`
* `--segment-overlap <SEGMENT_OVERLAP>` — Number of bases to overlap between segments

   Helpful to avoid missing variants at the edges of segments.

  Default value: `200`
* `--error-model <ERROR_MODEL>` — The error model to use

   This should match the sequencing platform used to generate the data

  Default value: `novaseq6000`

  Possible values:
  - `miseq`:
    MiSeq <https://support.illumina.com/sequencing/sequencing_instruments/miseq.html>
  - `miniseq`:
    MiniSeq <https://support.illumina.com/sequencing/sequencing_instruments/miniseq.html>
  - `nextseq500`:
    NextSeq 500 <https://support.illumina.com/sequencing/sequencing_instruments/nextseq-500.html>
  - `nextseq550`:
    NextSeq 550 <https://support.illumina.com/sequencing/sequencing_instruments/nextseq-550.html>
  - `hiseq2500`:
    HiSeq 2500 <https://support.illumina.com/sequencing/sequencing_instruments/hiseq_2500.html>
  - `novaseq6000`:
    NovaSeq 6000 <https://support.illumina.com/sequencing/sequencing_instruments/novaseq-6000.html>
  - `hiseq-x-ten`:
    HiSeq X Ten <https://support.illumina.com/sequencing/sequencing_instruments/hiseq-x.html>

* `--skip-methylation-calling` — Calculate threshold values and filters for methylation

  Default value: `false`
* `--vcf-threads <VCF_THREADS>` — Number of threads to use for writing (and compressing) VCF files

   This is subtracted from `--threads` but never below 1. Adjust this if you think that VCF writing is a bottleneck, e.g. when the output files contain a lot of positions.

  Default value: `3`
* `-@`, `--threads <TOTAL_THREADS>` — Number of threads to use for processing the BAM file. Will use all available threads when not specified.

   Note that VCF writing might use additional threads internally for compression. This can be overwritten with `--vcf-threads`.

  Default value: `14`



## `rastair per-read`

Call methylation per-read

This will produce a bed file that list the methylation status of all CpGs in every read that overlaps a CpG, plus some other metadata

**Usage:** `rastair per-read [OPTIONS] --fasta-file <FASTA_FILE> <BAM_FILE>`

###### **Arguments:**

* `<BAM_FILE>` — Path to sorted and indexed BAM file

###### **Filter Options:**

* `-f`, `--include-flags <INCLUDE_FLAGS>` — Include reads that match all of these bit-flags

  Default value: `3`
* `-F`, `--exclude-flags <EXCLUDE_FLAGS>` — Exclude reads that match any of these bit-flags

  Default value: `3852`
* `-w`, `--max-read-length <MAX_READ_LENGTH>` — expected maximum read length. If set too short, some read positions might not get counted. Safest to set this a bit higher than the actual read length, to allow for indels in reads

  Default value: `200`
* `-q`, `--min-mapq <MIN_MAPQ>` — Minimum mapping quality per aligned read

  Default value: `1`
* `--exclude-ambiguous` — Exclude reads where the orientation cannot be unambiguously determined
* `--count-clipped` — Count clipped positions

   By default, rastair ignores the leading (soft and hard) clipped positions in the "positions in read" columns. The indices written can be seen as "position in read relative to the first base actually aligned".

   If `--count-clipped` is set, clipped positions will instead be counted. The indices written then match the sequence of the read.

###### **Input Options:**

* `-r`, `--fasta-file <FASTA_FILE>` — Path to sorted and indexed (via samtools faidx) FASTA file. Can be bgzip compressed, but requires both a gzi index and a fai index
* `-l`, `--region <REGION>` — Restrict to a specific chromosome or region of a chromosome. Format is "chr", "chr:start" or "chr:start-end", where start is 1-based and end is inclusive
* `--calls <CALLS>` — BED file Rastair wrote with methylation calls per position

###### **Output Options:**

* `-A`, `--all-reads` — Report reads with no CpGs in them
* `--bed <BED>` — Output BED file with all reads

  Default value: `-`
* `--bed-format <BED_FORMAT>` — Format of the output BED reads file

   If not specified, the format is guessed based on the file extension.

  Possible values:
  - `bed-gz`:
    BGZIP compressed file, usually `.bed.gz`
  - `bed`:
    Regular BED file, usually `.bed`


###### **Processing Options:**

* `--segment-max-length <SEGMENT_MAX_LENGTH>` — Maximum length of a segment in bases

   Used for splitting work between threads. Tweak this to adjust memory usage.

  Default value: `100000`
* `--segment-overlap <SEGMENT_OVERLAP>` — Number of bases to overlap between segments

   Helpful to avoid missing variants at the edges of segments.

  Default value: `500`
* `-@`, `--threads <TOTAL_THREADS>` — Number of threads to use for processing the BAM file. Will use all available threads when not specified.

   Note that VCF writing might use additional threads internally for compression. This can be overwritten with `--vcf-threads`.

  Default value: `14`



## `rastair bam`

Add methylation information to BAM files

This will rewrite a BAM file to add methylation information and change the methylated positions in the sequence to their original base.

**Usage:** `rastair bam [OPTIONS] --fasta-file <FASTA_FILE> <BAM_FILE> <CALLS_FILE>`

###### **Arguments:**

* `<BAM_FILE>` — Path to sorted and indexed BAM file
* `<CALLS_FILE>` — Rastair's calls to determine methylation

###### **Input Options:**

* `-r`, `--fasta-file <FASTA_FILE>` — Path to sorted and indexed (via samtools faidx) FASTA file. Can be bgzip compressed, but requires both a gzi index and a fai index
* `-l`, `--region <REGION>` — Restrict to a specific chromosome or region of a chromosome. Format is "chr", "chr:start" or "chr:start-end", where start is 1-based and end is inclusive

###### **Output Options:**

* `-o`, `--output <OUTPUT>` — Output file

  Default value: `-`

###### **Processing Options:**

* `--segment-max-length <SEGMENT_MAX_LENGTH>` — Maximum length of a segment in bases

   Used for splitting work between threads. Tweak this to adjust memory usage.

  Default value: `100000`



## `rastair convert`

Convert between different file formats

**Usage:** `rastair convert [OPTIONS]`

###### **Filter Options:**

* `--bed-ml <ML_THRESHOLD>` — Minimum ML score to consider a position as variant

   This does nothing if the input data does not contain ML scores.

  Default value: `0.80`
* `--bed-include-empty` — Include CpG positions with zero coverage

   This can be useful to get a complete list of CpG positions in the output BED file. Note that this requires the input data to contain a complete list of CpG positions, e.g. by using the `--cpgs-only` option when calling methylation.

###### **Input Options:**

* `-i`, `--input <INPUT>` — Input file

  Default value: `-`
* `-f`, `--input-format <INPUT_FORMAT>` — Input file format, guessed from file extension if not specified

  Possible values:
  - `vcf`:
    Text-based VCF format (.vcf)
  - `bcf`:
    Binary VCF format (.bcf)
  - `vcf-compressed`:
    Compressed text-based VCF format (.vcf.gz)
  - `mpk.lz4`


###### **Output Options:**

* `-o`, `--output <OUTPUT>` — Output file

  Default value: `-`
* `-F`, `--output-format <OUTPUT_FORMAT>` — Output file format, guessed from file extension if not specified

  Possible values:
  - `vcf`:
    Text-based VCF format (.vcf)
  - `bcf`:
    Binary VCF format (.bcf)
  - `vcf-compressed`:
    Compressed text-based VCF format (.vcf.gz)
  - `mpk.lz4`
  - `bed`:
    Regular BED file, usually `.bed`
  - `bed-gz`:
    BGZIP compressed file, usually `.bed.gz`




## `rastair view`

View internal format as JSON lines

**Usage:** `rastair view [OPTIONS] <INPUT>`

###### **Arguments:**

* `<INPUT>` — Message Pack file to view

###### **Output Options:**

* `-o`, `--output <OUTPUT>` — Message Pack file to view

  Default value: `-`



## `rastair mbias`

Calculate conversion per base position in read

This will produce a `mbias.html` file with information about conversion counts relative to read position.

Please note that this is currently implemented as an R script. Unless you're using the official Docker image, you need to install R and the necessary packages yourself.

**Usage:** `rastair mbias [OPTIONS] <BED_FILE>`

###### **Arguments:**

* `<BED_FILE>` — Input per-read BED file (can be gzipped)

###### **Filter Options:**

* `--region <REGION>` — Genomic region
* `--include-flag <INCLUDE_FLAG>` — Include bitflag as integer

  Default value: `3`
* `--exclude-flag <EXCLUDE_FLAG>` — Exclude bitflag as integer

  Default value: `3852`
* `--read-length <READ_LENGTH>` — Read length as integer

###### **Output Options:**

* `--output-prefix <OUTPUT_PREFIX>` — Output path prefix

  Default value: `.`

###### **Processing Options:**

* `--tabix-path <TABIX_PATH>` — Path to tabix executable

  Default value: `tabix`



<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>
