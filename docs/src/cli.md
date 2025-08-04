# Command-Line Help for `rastair2`

This document contains the help content for the `rastair2` command-line program.

**Command Overview:**

* [`rastair2`↴](#rastair2)
* [`rastair2 call`↴](#rastair2-call)
* [`rastair2 convert`↴](#rastair2-convert)
* [`rastair2 view`↴](#rastair2-view)

## `rastair2`

Rastair2

Process TAPS-sequenced BAM files for methylation calling

**Usage:** `rastair2 [OPTIONS] <COMMAND>`

###### **Subcommands:**

* `call` — Call methylated positions
* `convert` — Convert between different file formats
* `view` — View internal format as JSON lines

###### **Options:**

* `-v`, `--verbose` — Enable more logging



## `rastair2 call`

Call methylated positions

**Usage:** `rastair2 call [OPTIONS] --fasta-file <FASTA_FILE> <BAM_FILE>`

###### **Arguments:**

* `<BAM_FILE>` — Path to sorted and indexed BAM file

###### **Options:**

* `-r`, `--fasta-file <FASTA_FILE>` — Path to sorted and indexed (via samtools faidx) FASTA file. Can be bgzip compressed, but requires both a gzi index and a fai index
* `-l`, `--region <REGION>` — Restrict to a specific chromosome or region of a chromosome. Format is "chr", "chr:start" or "chr:start-end", where start is 1-based and end is inclusive
* `--segment-max-length <SEGMENT_MAX_LENGTH>` — Maximum length of a segment in bases

  Default value: `100000`
* `--segment-overlap <SEGMENT_OVERLAP>` — Number of bases to overlap between segments

  Default value: `100`
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

* `--keep-overlapping-reads` — Whether to keep overlapping reads

  Default value: `false`
* `--cpgs-only` — Only look at sites that are CpG in the reference

  Default value: `false`
* `--nOT <N_OT>` — For OT reads, exclude `[r1_start, r1_end, r2_start, r2_end]` bases from counting.

   The coordinates are relative to the read, so start is the distance from the 5' of the read, the end is the distance to the 3', irrespective of which way around the read aligns to the reference.

   Also note that the distance is relative to read length, not alignment length, so soft-clipped bases count, too!

  Default value: `0,0,0,0`
* `--nOB <N_OB>` — For OB reads, exclude `[r1_start, r1_end, r2_start, r2_end]` bases from counting.

   The coordinates are relative to the read, so start is the distance from the 5' of the read, the end is the distance to the 3', irrespective of which way around the read aligns to the reference.

   Also note that the distance is relative to read length, not alignment length, so soft-clipped bases count, too!

  Default value: `0,0,0,0`
* `-f`, `--include-flags <INCLUDE_FLAGS>` — Include reads that match all of these bit-flags
* `-F`, `--exclude-flags <EXCLUDE_FLAGS>` — Exclude reads that match any of these bit-flags
* `--cpg-novo-min-depth <CPG_NOVO_MIN_DEPTH>` — Minimum reads needed in support of de-novo CpG

  Default value: `2`
* `--cpg-novo-min-baseq <CPG_NOVO_MIN_BASEQ>` — Minimum base quality for de-novo CpGs

  Default value: `15`
* `--cpg-novo-min-mapq <CPG_NOVO_MIN_MAPQ>` — Minimum mapping quality for de-novo CpGs

  Default value: `50`
* `--cpg-novo-min-vaf <CPG_NOVO_MIN_VAF>` — Minimum variant allele frequency for de-novo CpGs

  Default value: `0.2`
* `--skip-methylation-calling` — Calculate threshold values and filters for methylation

  Default value: `false`
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
* `--ml <ML>` — Use machine learning model with this threshold value to call variants and methylation events

   When specified, a ML model will classify positions with a prediction score. Anything above this threshold is considered PASS.
* `--model-cpg <MODEL_CPG>` — Path to the model for CpG positions

   Default is the bundled model in the Rastair binary.
* `--model-denovo-cpg <MODEL_DENOVO_CPG>` — Path to the model for de novo CpG positions

   Default is the bundled model in the Rastair binary.
* `--model-others <MODEL_OTHERS>` — Path to the model for other positions

   Default is the bundled model in the Rastair binary.
* `-o`, `--vcf <VCF>` — VCF/BCF output file path (use - to write to stdout)

   Format is guessed based on the file extension: `.vcf` for VCF (uncompressed), `.vcf.gz` for VCF (compressed), `.bcf` for BCF (compressed) `.mpk.lz4` for internal format (Message Pack, LZ4-compressed)
* `--vcf-threads <VCF_THREADS>` — Number of threads to use for writing (and compressing) VCF files

   This is subtracted from `--threads` but never below 1. Adjust this if you think that VCF writing is a bottleneck, e.g. when the output files contain a lot of positions.

  Default value: `3`
* `--bed <BED>` — Output BED file with the called methylated positions
* `--bed-format <BED_FORMAT>` — Format of the output BED file

   If not specified, the format is guessed based on the file extension.

  Possible values:
  - `bed-gz`:
    BGZIP compressed file, usually `.bed.gz`
  - `bed`:
    Regular BED file, usually `.bed`

* `-@`, `--threads <TOTAL_THREADS>` — Number of threads to use for processing the BAM file. Will use all available threads when not specified.

   Note that VCF writing might use additional threads internally for compression. This can be overwritten with `--vcf-threads`.

  Default value: `14`



## `rastair2 convert`

Convert between different file formats

**Usage:** `rastair2 convert [OPTIONS]`

###### **Options:**

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




## `rastair2 view`

View internal format as JSON lines

**Usage:** `rastair2 view [OPTIONS] <INPUT>`

###### **Arguments:**

* `<INPUT>` — Message Pack file to view

###### **Options:**

* `-o`, `--output <OUTPUT>` — Message Pack file to view

  Default value: `-`



<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>
