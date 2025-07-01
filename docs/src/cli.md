# Command-Line Help for `rastair2`

This document contains the help content for the `rastair2` command-line program.

**Command Overview:**

* [`rastair2`↴](#rastair2)
* [`rastair2 call`↴](#rastair2-call)
* [`rastair2 convert`↴](#rastair2-convert)

## `rastair2`

Rastair2

Process TAPS-sequenced BAM files for methylation calling

**Usage:** `rastair2 [OPTIONS] <COMMAND>`

###### **Subcommands:**

* `call` — Call methylated positions
* `convert` — Convert between different file formats

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
* `--cpg-novo-min-depth <CPG_NOVO_MIN_DEPTH>` — Minimum reads needed in support of de-novo CpG

  Default value: `2`
* `--cpg-novo-min-baseq <CPG_NOVO_MIN_BASEQ>` — Minimum base quality for de-novo CpGs

  Default value: `15`
* `--cpg-novo-min-mapq <CPG_NOVO_MIN_MAPQ>` — Minimum mapping quality for de-novo CpGs

  Default value: `50`
* `--cpg-novo-min-vaf <CPG_NOVO_MIN_VAF>` — Minimum variant allele frequency for de-novo CpGs

  Default value: `0.2`
* `--calling <CALLING>` — The methylation calling mode

  Default value: `none`

  Possible values:
  - `none`:
    Don't perform methylation calling
  - `thresholds`:
    Call methylation events based on thresholds

* `--vaf-min <VAF_MIN>` — The minimum variant allele frequency

  Default value: `0`
* `--reads-min <READS_MIN>` — The minimum number of reads to call a variant

  Default value: `3`
* `-o`, `--vcf-output <VCF_OUTPUT>` — VCF/BCF output file path (use - to write to stdout)

   Format is guessed based on the file extension: `.vcf` for VCF (uncompressed), `.vcf.gz` for VCF (compressed), `.bcf` for BCF (compressed) `.mpk.lz4` for internal format (Message Pack, LZ4-compressed)

  Default value: `-`
* `--vcf-threads <VCF_THREADS>` — Number of threads to use for writing (and compressing) VCF files

   This is subtracted from `--threads` but never below 1

  Default value: `2`
* `-@`, `--threads <THREADS>` — Number of threads to use for processing the BAM file. Will use all available threads when not specified.

   Note that VCF writing might use additional threads internally for compression. This can be overwritten with `--vcf-threads`.

  Default value: `14`



## `rastair2 convert`

Convert between different file formats

**Usage:** `rastair2 convert [OPTIONS] --input <INPUT> --output <OUTPUT>`

###### **Options:**

* `--input <INPUT>` — Input file
* `--input-format <INPUT_FORMAT>` — Input file format, guessed from file extension if not specified

  Possible values:
  - `vcf`:
    Text-based VCF format
  - `bcf`:
    Binary VCF format (BCF)
  - `vcf-compressed`:
    Compressed text-based VCF format
  - `mpk.lz4`

* `-o`, `--output <OUTPUT>` — Output file
* `--output-format <OUTPUT_FORMAT>` — Output file format, guessed from file extension if not specified

  Possible values:
  - `vcf`:
    Text-based VCF format
  - `bcf`:
    Binary VCF format (BCF)
  - `vcf-compressed`:
    Compressed text-based VCF format
  - `mpk.lz4`
  - `bed`




<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>
