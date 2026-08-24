# Command-Line Help for `rastair`

This document contains the help content for the `rastair` command-line program.

**Version:** `2.2.0-rc.1`

**Command Overview:**

* [`rastair`↴](#rastair)
* [`rastair call`↴](#rastair-call)
* [`rastair per-read`↴](#rastair-per-read)
* [`rastair bam`↴](#rastair-bam)
* [`rastair bam standard`↴](#rastair-bam-standard)
* [`rastair bam legacy`↴](#rastair-bam-legacy)
* [`rastair convert`↴](#rastair-convert)
* [`rastair view`↴](#rastair-view)
* [`rastair mbias`↴](#rastair-mbias)
* [`rastair license`↴](#rastair-license)

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
* `license` — Show license -- rastair is licensed under a non-commercial use licence

###### **Options:**

* `-v`, `--verbose` — Enable more logging

   You can also use the `RASTAIR_LOG` environment variable to configure logging in a more precise way.

   Note that trace-level logging is disabled in production builds.



## `rastair call`

Call methylated positions

Process TAPS-sequenced BAM files and call methylated positions.

If no output file is specified, the output is written to stdout. You can use `--vcf` and `--bed` to write to files instead.

If using `-c` (`--cpgs-only`), all CpG positions in the reference as well as de-novo CpGs are written. Stdout will default to BED.

Only variants that pass all filters are written by default. Use `--all` to get a full VCF file.

**Usage:** `rastair call [OPTIONS] --fasta-file <FASTA_FILE> <BAM_FILE>`

###### **Arguments:**

* `<BAM_FILE>` — Path to sorted and indexed BAM or CRAM file

###### **Filter Options:**

* `--unpaired` — Enable unpaired mode

   In this mode, unpaired reads are accepted and strand assignment uses only the read's alignment direction (forward=OT, reverse=OB).

  Default value: `false`
* `--keep-overlapping-reads` — Whether to keep overlapping paired-end reads

   In unpaired (`--single-strand`) mode this is ignored because read-pair overlap deduplication is disabled.

  Default value: `false`
* `--v-min-depth <V_MIN_DEPTH>`

  Default value: `3`
* `--max-coverage <MAX_COVERAGE>`

  Default value: `1000`
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
* `--no-ml` — Only use hard thresholds to call variants and methylation events.

   This disables using the machine learning models. This will make rastair much faster, but at the cost of accuracy.
* `--ml <ML>` — Use machine learning model with this threshold value to call variants and methylation events

   When specified, a ML model will classify positions with a prediction score. Anything above this threshold is considered PASS.

   For consistency with `--no-ml`, this option can be also be specified as `--ml` without a value, which will use the default threshold.

  Default value: `0.50`
* `--model <MODEL>` — Path to the combined model file containing CpG, denovo, and others models

   Default is the bundled model in the Rastair binary.
* `-c`, `--cpgs-only` — Report CpGs only and default to BED output

   Only report positions that are CpGs in the reference or variants that would result in a de-novo CpG.

   If combined with `--all`, non-passing de-novo CpG positions and CpGs in the reference but without coverage in the sample will also be reported.

  Default value: `false`
* `--bed-include-empty` — Include CpG positions with zero coverage

   This can be useful to get a complete list of CpG positions in the output BED file. Note that this requires the input data to contain a complete list of CpG positions, e.g. by using the `--cpgs-only` option when calling methylation.

###### **Indel Options:**

* `--experimental-indels <EXPERIMENTAL_INDELS>` — Enable experimental indel calling

   Defaults to ML scoring.

  Possible values:
  - `ml`:
    The model decides: a low score vetoes an allele the hard chain would keep. Best precision, worst recall
  - `no-ml`:
    A fixed filter chain only, with no ML scoring of indels at call time
  - `ml-rescue`:
    The hard chain's calls all stand, and the model reinstates an allele the chain dropped as homozygous reference. Best recall

* `--min-indel-ao <MIN_INDEL_AO>` — Minimum alternate observations to call an indel

   This threshold is evaluated per indel allele after read-level indel filters are applied.

  Default value: `2`
* `--min-indel-depth <MIN_INDEL_DEPTH>` — Minimum depth to call an indel

   Depth is computed at the locus after applying indel-specific read filters and depth adjustments.

  Default value: `2`
* `--indel-max-mismatches <INDEL_MAX_MISMATCHES>` — Maximum number of non-TAPS mismatches allowed on a read supporting an indel.

   C->T mismatches on OT reads and G->A mismatches on OB reads are excluded from the count as they are expected TAPS methylation signal.

  Default value: `5`
* `--indel-end-of-read-cutoff <INDEL_END_OF_READ_CUTOFF>` — Ignore indels within this many bases of either read end.

   This stricter end-of-read filter helps suppress alignment artifacts at read starts and ends.

  Default value: `0`

###### **Input Options:**

* `-r`, `--fasta-file <FASTA_FILE>` — Path to sorted and indexed (via samtools faidx) FASTA file. Can be bgzip compressed, but requires both a gzi index and a fai index
* `-l`, `--region <REGIONS>` — Restrict processing to specific genomic regions.

   Accepts either space-separated region strings or a single BED file: - Region strings: `chr`, `chr:start`, `chr:start-end` (1-based inclusive) - Multiple regions separated by whitespace: `"chr1 chr2:100-200"` - BED files with `@` prefix: `@targets.bed`
* `--require-tags <REQUIRE_TAGS>` — Require reads to have a specific SAM tag value

   Format: TAG=VALUE, e.g. `--require-tags RG=mygroup`. Accepts one or more values (space-separated). A read is kept only if it matches all of the specified tag=value pairs.

###### **Methylation Options:**

* `--m-vaf-min <M_VAF_MIN>` — The minimum variant allele frequency

  Default value: `0.2`
* `--m-min-depth <M_MIN_DEPTH>` — The minimum number of reads to call a position as methylated

  Default value: `3`
* `--m-bq-ratio-min <M_BQ_RATIO_MIN>` — The minimum quality ratio `(ad_alt*bq_alt + 1) / (ad_ref*bq_ref + 1)`

  Default value: `0.27`
* `--m-read-position-min <M_READ_POSITION_MIN>` — The minimum relative position in read for alt allele evidence

  Default value: `0.2`
* `--m-read-position-max <M_READ_POSITION_MAX>` — The maximum relative position in read for alt allele evidence

  Default value: `0.8`
* `--m-max-coverage <M_MAX_COVERAGE>` — The maximum coverage depth for methylation calling

  Default value: `1000`

###### **Output Options:**

* `--all` — Output all positions, even if they do not pass filters.

   If combined with `--cpgs-only`, only CpG positions will be reported, including non-passing de-novo CpGs, and those without coverage.
* `-o`, `--vcf <VCF>` — VCF/BCF output file path (use - to write to stdout)

   Format is guessed based on the file extension: `.vcf` for VCF (uncompressed), `.vcf.gz` for VCF (compressed), `.bcf` for BCF (compressed) `.mpk.lz4` for internal format (Message Pack, LZ4-compressed)
* `--vcf-info-fields <VCF_INFO_FIELDS>` — Additional INFO fields to include in VCF output (comma-separated VCF field IDs)

   By default, only a minimal set is included.

  Possible values: `AD`, `BQ`, `DP`, `MQ`, `MQ0`, `NS`, `AS_SB`, `SC5`, `AF`, `ABQ`, `AMQ`, `AS_SS_BQ`, `AS_SS_MQ`, `PIR`, `ENT100`, `NAB`, `NOI`, `M5mC_Strands`, `CPG`, `CPGnovo`

* `--vcf-format-fields <VCF_FORMAT_FIELDS>` — Additional FORMAT fields to include in VCF output (comma-separated VCF field IDs)

   By default, only a minimal set is included.

  Possible values: `GT`, `GL`, `GC`, `DP`, `M5mC`, `DPM5mC`, `ADM5mC`, `ML`

* `--vcf-all-fields`

  Default value: `false`
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
* `--guess-read-orientation` — Guess OT/OB read orientation from mismatch motifs instead of SAM flags

   Scans read mismatches against the reference and counts `TG` versus `CA` motifs in a 2 bp window anchored at each mismatch (current+next and previous+current), using the htslib/reference-oriented read sequence. `TG > CA` means OT, `CA > TG` means OB, and ties / evidence-free reads are split pseudo-randomly but reproducibly per read.

   Useful for non-directional libraries, commonly from tagmentation prep.

   This option currently affects only the pileup-based `call` workflow.

  Default value: `false`
* `--error-model <ERROR_MODEL>` — The error model to use

   Accepts platform names or a custom error rate (e.g., 0.005)

  Default value: `novaseq6000`

  Possible values:
  - `miseq`:
    MiSeq <https://support.illumina.com/sequencing/sequencing_instruments/miseq.html>
  - `miniseq`:
    MiniSeq <https://support.illumina.com/sequencing/sequencing_instruments/miniseq.html>
  - `nextseq500`:
    NextSeq500 <https://support.illumina.com/sequencing/sequencing_instruments/nextseq-500.html>
  - `nextseq550`:
    NextSeq550 <https://support.illumina.com/sequencing/sequencing_instruments/nextseq-550.html>
  - `hiseq2500`:
    HiSeq2500 <https://support.illumina.com/sequencing/sequencing_instruments/hiseq_2500.html>
  - `novaseq6000`:
    NovaSeq6000 <https://support.illumina.com/sequencing/sequencing_instruments/novaseq-6000.html>
  - `hiseqxten`:
    HiSeq X Ten <https://support.illumina.com/sequencing/sequencing_instruments/hiseq-x.html>

* `--linear-dedup-threshold <LINEAR_DEDUP_THRESHOLD>` — Depth threshold below which linear name dedup is used instead of a hashmap

   At pileup positions with depth ≤ this value, read name deduplication uses a linear scan through parallel suffix/name arrays rather than an `FxHashMap`. Set to 0 to always use the hashmap.

  Default value: `30`
* `--indel-error-rate <INDEL_ERROR_RATE>` — Error rate for indel genotyping (higher than SNV due to alignment uncertainty)

  Default value: `0.05`
* `--gpu` — Use GPU-accelerated ML predictions which might speed up large datasets, but can be slower for small ones due to overhead.

   Requires a Metal/Vulkan/DX12-capable GPU.
* `--vcf-threads <VCF_THREADS>` — Number of threads to use for writing (and compressing) VCF files

   This is subtracted from `--threads` but never below 1. Adjust this if you think that VCF writing is a bottleneck, e.g. when the output files contain a lot of positions.

  Default value: `1`
* `-@`, `--threads <TOTAL_THREADS>` — Number of threads to use for processing the BAM file. Will use all available threads when not specified.

   Note that VCF writing might use additional threads internally for compression. This can be overwritten with `--vcf-threads`.

  Default value: `10`
  [env: `RASTAIR_THREADS`]



## `rastair per-read`

Call methylation per-read

This will produce a bed file that list the methylation status of all CpGs in every read that overlaps a CpG, plus some other metadata

**Usage:** `rastair per-read [OPTIONS] --fasta-file <FASTA_FILE> <BAM_FILE>`

###### **Arguments:**

* `<BAM_FILE>` — Path to sorted and indexed BAM or CRAM file

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
* `--unpaired` — Enable unpaired mode

   In this mode, unpaired reads are accepted and orientation is inferred from alignment direction only (forward=OT, reverse=OB).

  Default value: `false`
* `--count-clipped` — Count clipped positions

   By default, rastair ignores the leading (soft and hard) clipped positions in the "positions in read" columns. The indices written can be seen as "position in read relative to the first base actually aligned".

   If `--count-clipped` is set, clipped positions will instead be counted. The indices written then match the sequence of the read.

###### **Input Options:**

* `-r`, `--fasta-file <FASTA_FILE>` — Path to sorted and indexed (via samtools faidx) FASTA file. Can be bgzip compressed, but requires both a gzi index and a fai index
* `-l`, `--region <REGIONS>` — Restrict processing to specific genomic regions.

   Accepts either space-separated region strings or a single BED file: - Region strings: `chr`, `chr:start`, `chr:start-end` (1-based inclusive) - Multiple regions separated by whitespace: `"chr1 chr2:100-200"` - BED files with `@` prefix: `@targets.bed`
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
* `--guess-read-orientation` — Guess OT/OB read orientation from mismatch motifs instead of SAM flags

   Scans read mismatches against the reference and counts `TG` versus `CA` motifs in a 2 bp window anchored at each mismatch (current+next and previous+current), using the htslib/reference-oriented read sequence. `TG > CA` means OT, `CA > TG` means OB, and ties / evidence-free reads are split pseudo-randomly but reproducibly per read.

  Default value: `false`
* `-@`, `--threads <TOTAL_THREADS>` — Number of threads to use for processing the BAM file. Will use all available threads when not specified.

   Note that VCF writing might use additional threads internally for compression. This can be overwritten with `--vcf-threads`.

  Default value: `10`



## `rastair bam`

Add methylation information to BAM files

Writes a new BAM file that includes methylation tags derived from Rastair calls.

**Usage:** `rastair bam <COMMAND>`

###### **Subcommands:**

* `standard` — Write modBAM with MM/ML tags as specified by the SAM 4.5 spec This will rewrite SEQ to un-modify bases that have methylation evidence
* `legacy` — Write BAM with "legacy" XR/XG/XM tags, compatible with tools like DRAGEN and Bismark



## `rastair bam standard`

Write modBAM with MM/ML tags as specified by the SAM 4.5 spec This will rewrite SEQ to un-modify bases that have methylation evidence

**Usage:** `rastair bam standard [OPTIONS] --fasta-file <FASTA_FILE> <BAM_FILE> <CALLS_FILE>`

###### **Arguments:**

* `<BAM_FILE>` — Path to sorted and indexed BAM or CRAM file
* `<CALLS_FILE>` — Rastair's calls to determine methylation

###### **Input Options:**

* `-r`, `--fasta-file <FASTA_FILE>` — Path to sorted and indexed (via samtools faidx) FASTA file. Can be bgzip compressed, but requires both a gzi index and a fai index
* `-l`, `--region <REGIONS>` — Restrict processing to specific genomic regions.

   Accepts either space-separated region strings or a single BED file: - Region strings: `chr`, `chr:start`, `chr:start-end` (1-based inclusive) - Multiple regions separated by whitespace: `"chr1 chr2:100-200"` - BED files with `@` prefix: `@targets.bed`

###### **Output Options:**

* `-o`, `--output <OUTPUT>` — Output file

  Default value: `-`

###### **Processing Options:**

* `--segment-max-length <SEGMENT_MAX_LENGTH>` — Maximum length of a segment in bases

   Used for splitting work between threads. Tweak this to adjust memory usage.

  Default value: `100000`
* `-@`, `--threads <THREADS>` — Number of threads to use for processing the BAM file

  Default value: `10`



## `rastair bam legacy`

Write BAM with "legacy" XR/XG/XM tags, compatible with tools like DRAGEN and Bismark

**Usage:** `rastair bam legacy [OPTIONS] --fasta-file <FASTA_FILE> <BAM_FILE> <CALLS_FILE>`

###### **Arguments:**

* `<BAM_FILE>` — Path to sorted and indexed BAM or CRAM file
* `<CALLS_FILE>` — Rastair's calls to determine methylation

###### **Input Options:**

* `-r`, `--fasta-file <FASTA_FILE>` — Path to sorted and indexed (via samtools faidx) FASTA file. Can be bgzip compressed, but requires both a gzi index and a fai index
* `-l`, `--region <REGIONS>` — Restrict processing to specific genomic regions.

   Accepts either space-separated region strings or a single BED file: - Region strings: `chr`, `chr:start`, `chr:start-end` (1-based inclusive) - Multiple regions separated by whitespace: `"chr1 chr2:100-200"` - BED files with `@` prefix: `@targets.bed`

###### **Output Options:**

* `-o`, `--output <OUTPUT>` — Output file

  Default value: `-`

###### **Processing Options:**

* `--segment-max-length <SEGMENT_MAX_LENGTH>` — Maximum length of a segment in bases

   Used for splitting work between threads. Tweak this to adjust memory usage.

  Default value: `100000`
* `-@`, `--threads <THREADS>` — Number of threads to use for processing the BAM file

  Default value: `10`



## `rastair convert`

Convert between different file formats

**Usage:** `rastair convert [OPTIONS]`

###### **Filter Options:**

* `-c`, `--cpgs-only` — Report CpGs only and default to BED output

   Only report positions that are CpGs in the reference or variants that would result in a de-novo CpG.

   If combined with `--all`, non-passing de-novo CpG positions and CpGs in the reference but without coverage in the sample will also be reported.

  Default value: `false`
* `--bed-include-empty` — Include CpG positions with zero coverage

   This can be useful to get a complete list of CpG positions in the output BED file. Note that this requires the input data to contain a complete list of CpG positions, e.g. by using the `--cpgs-only` option when calling methylation.
* `--bed-ml <ML_THRESHOLD>` — Minimum ML score to consider a position as variant

   This does nothing if the input data does not contain ML scores.

  Default value: `0.50`

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

* `--all` — Output all positions, even if they do not pass filters.

   If combined with `--cpgs-only`, only CpG positions will be reported, including non-passing de-novo CpGs, and those without coverage.

###### **Processing Options:**

* `--error-model <ERROR_MODEL>`

  Default value: `novaseq6000`

  Possible values:
  - `miseq`:
    MiSeq <https://support.illumina.com/sequencing/sequencing_instruments/miseq.html>
  - `miniseq`:
    MiniSeq <https://support.illumina.com/sequencing/sequencing_instruments/miniseq.html>
  - `nextseq500`:
    NextSeq500 <https://support.illumina.com/sequencing/sequencing_instruments/nextseq-500.html>
  - `nextseq550`:
    NextSeq550 <https://support.illumina.com/sequencing/sequencing_instruments/nextseq-550.html>
  - `hiseq2500`:
    HiSeq2500 <https://support.illumina.com/sequencing/sequencing_instruments/hiseq_2500.html>
  - `novaseq6000`:
    NovaSeq6000 <https://support.illumina.com/sequencing/sequencing_instruments/novaseq-6000.html>
  - `hiseqxten`:
    HiSeq X Ten <https://support.illumina.com/sequencing/sequencing_instruments/hiseq-x.html>

* `-@`, `--threads <TOTAL_THREADS>` — Total number of threads to use (e.g. for parallel compression)

  Default value: `10`
  [env: `RASTAIR_THREADS`]



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

**Usage:** `rastair mbias [OPTIONS] <--bed <BED_FILE>|--bam <BAM_FILE>>`

###### **Filter Options:**

* `--region <REGION>` — Genomic region
* `--include-flag <INCLUDE_FLAG>` — Include bitflag as integer

  Default value: `3`
* `--exclude-flag <EXCLUDE_FLAG>` — Exclude bitflag as integer

  Default value: `3852`
* `--read-length <READ_LENGTH>` — Read length as integer

###### **Input Options:**

* `--bed <BED_FILE>` — Input per-read BED file (must be tabix indexed or indexable)
* `--bam <BAM_FILE>` — Input BAM file
* `--reference <REFERENCE>` — Reference FASTA file (required for V-bias and GC/CpG bias plots)
* `--vcf <VCF>` — VCF file with methylation calls (required for GC/CpG bias plots)

###### **Options:**

* `--r-script-dir <R_SCRIPT_DIR>` — Override directory to find R scripts

   When not set, tries to look for `$rastair_path/scripts` and `./scripts`
  [env: `R_SCRIPT_DIR`]

###### **Output Options:**

* `--output-prefix <OUTPUT_PREFIX>` — Output path prefix

  Default value: `.`

###### **Processing Options:**

* `--no-vbias` — Do not generate V-bias plots (faster)
* `--no-gc` — Do not generate GC/CpG bias plots
* `--tabix-path <TABIX_PATH>` — Path to tabix executable

  Default value: `tabix`
* `--bcftools-path <BCFTOOLS_PATH>` — Path to bcftools executable

  Default value: `bcftools`
* `--rastair-path <RASTAIR_PATH>` — Path to rastair executable

  Default value: `rastair`
* `--threads <THREADS>` — Number of threads to use

  Default value: `1`
* `--wgbs` — Treat the input as inverted, i.e. mod=unmod and unmod=mod



## `rastair license`

Show license -- rastair is licensed under a non-commercial use licence

**Usage:** `rastair license`



<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>
