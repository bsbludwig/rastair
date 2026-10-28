# Changelog

This is the changelog for Rastair 2.

## Unreleased

Calling:

- Pileups are read by [seqair](https://github.com/Softleif/seqair) instead of htslib.
  Along with some other changes, this should be another great performance improvement for Rastair.
- `call` calls indels in every run, with a TAPS-aware pair-HMM instead of the ML model; `--experimental-indels` is gone, and the bundled model no longer contains indel models.
  Indels at most 20 bp apart are called together, against haplotypes carrying every combination of them, with positions showing a nearby SNV left open.
  An indel allele needs at least 5% of the column's depth, besides `--min-indel-ao`, to be scored.
- An SNV and an indel at one position can sit on the same haplotype, as a complex allele often is spelled; reads carrying both were charged against every haplotype.
- An indel record's `GT` is the likeliest genotype of its own `PL` among those carrying a called allele, so `GT`, `GQ` and `PL` agree.
  Its `AD` counts the reads that favour the record's reference and each of its alleles, and `DP` those plus the reads of another allele at the same position; a read that fits two alleles about equally counts for neither.
- Removed `--min-indel-depth`, `--indel-error-rate`, `--indel-max-mismatches` and `--indel-end-of-read-cutoff`.
- The `indel_strand` FILTER is replaced by `low_indel_qual`.
- The pair-HMM also re-scores the model's SNV candidates: alts the model scores at 0.005 or more, the two it scores highest at a position, are scored against the haplotypes, and a second set of forests in the bundled model scores them again with that evidence. On by default; `--snv-engine=false` scores with the first forests alone.
- SNV genotypes come from the same pair-HMM; the read-count genotyper is gone. An SNV's `GT` is the likeliest genotype of its own `PL` among those carrying a called alt, so `GT`, `GQ` and `PL` agree.
- Each SNV candidate is scored by the pair-HMM together with its nearest neighbours within 20 bp, on either side. Runs of more than three candidates were cut left to right, so an SNV's `GQ` and `PL` depended on where the cut fell: one scored beside a linked neighbour has `0/0` far less likely than one scored alone.
- An SNV alternate is called only when the model passes it and the pair-HMM's genotype holds a copy: when the model passes two alternates but the genotype holds one, when it passes a third, or when the pair-HMM could not score it (an `N` reference, no usable read), it is not in the passing record; under `--all` it is written separately with the new `no_genotype_copy` FILTER.
- Two called SNVs next to each other that form a CpG (`TA>CG`) are a de-novo CpG, for genotyping, methylation and `CPGnovo` alike, when the pair-HMM puts the two on one chromosome. A called SNV makes no CpG with a reference neighbour whose own called SNV the pair-HMM puts on the same chromosome (`TG>CA`, `CA>TG`).
- A called insertion that puts a `C` before a reference `G`, or a `G` after its reference `C` anchor, makes a de-novo CpG as a called SNV does: the reference base gets `CPGnovo` and an `M5mC` from the reads whose neighbour is the inserted base, and the model scores a `G>A` or `C>T` there as methylation, not as an SNV. Such a site's methylation was called a PASS `G>A` or `C>T`. A deletion joining a `C` and a `G` makes no CpG yet.
- `call` links heterozygous variants that share read pairs into phase blocks, written as phased genotypes (`0|1`) with new `PS` and `PQ` FORMAT fields; each block has its own `PS`. Each read's vote comes from the pair-HMM's likelihoods, so a converted read at a CpG casts none. `--phase-min-switch-quality` sets where a block is cut.
- With `--gpu`, the pair-HMM scores on the GPU as well; the calls are the same as on the CPU.
- After a GPU failure (an error, a crash of the inference thread, or fewer scores than it was sent), the region and the rest of the run are scored on the CPU. An alt whose features cannot be computed is left unscored (`low_ml_score`).
- Removed `--no-ml`. Without the model the pair-HMM had to decide every SNV, which was slower than the model (2.6x the wall time at 26x) and less accurate (−0.04 SNV F1 at 26x, −0.07 at 10x).
- The model alone decides which alleles pass. Removed the hard thresholds `--v-min-depth`, `--m-vaf-min`, `--m-min-depth`, `--m-bq-ratio-min`, `--m-read-position-min`, `--m-read-position-max` and the four `--cpg-novo-min-*` options, and their FILTER codes (`lowDp`, `m_vaf`, `m_bq_ratio`, `m_pos`, `dnCpG_*`). Rejected records under `--all` carry `pre_ml` or `low_ml_score`, and one with too little evidence to score no longer carries both.
- `--m-max-coverage` stays, and a methylation candidate (a `C>T` at a CpG's `C`, a `G>A` at its `G`) that more reads show is rejected as `m_highDp` whatever the model scores it. It was a FILTER code next to the model's verdict, so such an alt could pass.
- The default `--ml` threshold is 0.4 instead of 0.5, which recovers SNVs at low depth (+0.0024 SNV F1 at 10x) and is neutral at 26x and 73x.
- `--ml` below 0.005 is an error: an SNV is called only with the pair-HMM's genotype, and the pair-HMM scores only alts the model scores that high.
- `ml train` and `call` score the same rows: training skips the alts `call` filters before the model and rows with an infinite feature, and both take a neighbouring column only when it is the adjacent position.
- Three of the model's features at a CpG's neighbour read the wrong strand or column: at a de-novo `X>G` CpG both β features read the strand that shows no methylation there, at a reference CpG's `G` the β was 1−β, and the neighbour's alt score weighed its reference reads by this position's base quality.
- The model's β feature at a de-novo CpG counts the converted reads also where they read as the reference: it was 0 at every `T>C` and `A>G` de-novo CpG.
- The bundled model is retrained on the fixes above and the pair-HMM's evidence, in a flat file format: a model trained with an earlier version cannot be loaded and has to be retrained with `rastair ml train`.
- `ml train` no longer accepts `--ml-features`; the standard features are the only set.
- `ml train` takes one or more BAMs before the truth VCF, pooled into one training collection, so one library down-sampled to several depths trains one model.
- A mapping quality of 255 ("not available" in the SAM spec) counts as 60 in `MQ`, `AMQ`, the model's features and the pair-HMM; it counted as a near-certain mapping. The per-read BED still writes 255.
- `--max-coverage` counts the reads that pass the flag, mapping-quality and tag filters, so a pileup of MAPQ 0 reads cannot crowd out the usable ones; `0` means no limit.
- `--require-tags` compares the value as the tag's SAM type: integers and floats by value (`NM=0` matches `NM:i:0`), strings exactly. It compared every tag as a string, so a typed tag such as `HP:i:1` matched no read.
- With `--guess-read-orientation`, an aligned base outside the reference window is skipped, as `per-read` does; it discarded the read's whole guess, so reads longer than the segment overlap fell back to their flag strand.
- New `--rescue-soft-clip-cpg`: a single soft-clipped base that is the missing partner of a reference CpG is counted as a normal observation.
- New `--segment-max-bytes`: a segment whose index-estimated size exceeds this budget is read in smaller pieces.
- `--segment-max-length` accepts 1,000 to 1,000,000; 0 made `per-read` and `bam` allocate without bound.
- Removed `--linear-dedup-threshold`: it only tuned the htslib pileup path.

Output:

- The `GL` and `GC` FORMAT fields are replaced by the standard `GQ` (Integer, capped at 99) and `PL` (`Number=G`), which vcfeval, hap.py, GLnexus and imputation tools read. Both come from one vector, the record's genotype likelihoods over its own alleles: `GQ` is the chance the written genotype is wrong among all of them, `0/0` included, so a genotype the reads doubt says so; a two-alt record has one `GQ` over its six genotypes. A reference-only record has no `PL`, and its `GQ` is the chance it is het or hom-alt, from its reads. The BED `gt_p_score` column repeats `gt_conf_score`, which takes the `GQ`. `--vcf-format-fields GL` and `GC` are errors.
- A position with no reads to weigh (an uncovered CpG under `--all`) is written `./.` with no `GQ` or `QUAL`, instead of a confident `0/0`; `call --bed` and `convert` write its `GQ` as 0.
- The `0/0` `GQ` at the partner of a de-novo CpG candidate weighs the CpG's reads only once that candidate is called; at the partner of a rejected candidate every read weighs, as anywhere else.
- `QUAL` is the confidence of whatever decided the record: the model's score for an SNV record, the pair-HMM's margin for an indel record (a likelihood ratio, not a probability), and for a reference-only record its `GQ` instead of a constant 99.
- `ML` is `.` for an alt the model did not score, instead of 0.
- Multi-allelic records list their `ALT` alleles by read depth, deepest first.
- `M5mC` (and `DPM5mC`, `ADM5mC`) are written exactly where `CPG` or `CPGnovo` is set ([#12](https://github.com/bsbludwig/rastair/issues/12)).
- Reference-only records (`ALT=.`) are only written at CpG and de-novo CpG positions.
- `M5mC`, `DPM5mC` and `ADM5mC` hold one value per CpG a position takes part in, the reference CpG first and the de-novo one second. A reference CpG without reads dropped its slot, so the de-novo CpG's β was read as the reference one's.
- A CpG that no read shows has no β: `M5mC` is `.` and the BED `beta_est` is empty, instead of 0, so averaging no longer counts it unmethylated. `DPM5mC` and `ADM5mC` stay 0. `bam` and `per-read` read such a BED row as a CpG without a β.
- The unmodified and modified counts of `M5mC_Strands` and the BED `unmod` and `mod` columns are the reads β is computed from, those whose neighbour in the read is the CpG partner. They counted every read with that base on the strand, so β could not be recomputed from them.
- `M5mC_Strands` is `0,0,0,0` at a position that is no CpG, next to a de-novo candidate that was not called, and at a de-novo CpG it counts the side the called alt creates.
- `--all` reports every reference CpG, also those without coverage and also without `--cpgs-only`; with `--cpgs-only` it reports the de-novo CpGs that were called, not the rejected candidates. `call --bed-include-empty` is gone: `--all` writes the uncovered CpGs to the BED as to the VCF.
- CpG-only output (`-c`, or `--bed` alone) reports exactly the CpGs of a full run: it does the work of a full run and drops the other records only when writing. It writes no indel records.
- `convert` no longer accepts `--error-model`, `--all`, `-c`/`--cpgs-only` and `--bed-ml`, which it never used.
- `convert` takes `M5mC` as the VCF writes it; it set it to 0 at a CpG with a heterozygous C>T or G>A SNV. It writes one BED row per position, also from an `--all` VCF, which gave duplicates.
- Removed the MessagePack output format and the `view` subcommand.
- `call` writes a CSI index next to a `.vcf.gz` or `.bcf` output.
- `bam` replaces existing `XM`, `XR` and `XG` tags instead of failing ([#15](https://github.com/bsbludwig/rastair/issues/15)).
- `per-read` no longer writes duplicate rows ([#18](https://github.com/bsbludwig/rastair/issues/18)).
- `bam` annotates a CpG whose BED β is 0. It took every such row for an SNP, so a fully unmethylated CpG had no `XM` (`z`), and a CpG whose few methylated reads round β to `0.00` lost their MM/ML entries. A BED row is now an SNP only where its genotype holds no copy of the CpG's base.
- `per-read`'s `mod_denovos` and `unmod_denovos` go by each read's base, like `mod_cpgs`, `unmod_cpgs` and `bam`'s tags: `T` or `C` on an OT read, `A` or `G` on an OB read. They listed every read covering a de-novo CpG by the site's β (methylated above 0.5), whatever its base or strand.
- `per-read`'s `num_cpg` counts the de-novo CpGs a read shows as well, and `num_mod` the methylated ones, so a read whose only CpGs are de-novo is written without `--all-reads`. Both counted reference CpGs only.

Regions:

- A `--region`'s first base is processed (it was skipped, as was the first base of a whole contig), and a G on a region's first base sees the base before it, so a CpG there is found.
- Overlapping or adjacent `--region`s are merged; they were processed, and written, twice.
- A region BED line with fewer than three fields is an error with its line number, like a line with a bad number; it was skipped.
- A BAM contig missing from the FASTA (a decoy or alt contig) is skipped instead of failing the run; one a `--region` names is still an error.

Running:

- Progress is shown in the terminal tab/taskbar (OSC 9;4) in terminals that support it, e.g. Windows Terminal, Ghostty, WezTerm, iTerm2 and Konsole.
- `call`, `per-read` and `bam` write their outputs (and indices) as `<name>.partial` and rename them only when the run succeeded, so a file under its final name is always complete.
- A crash (panic) in any thread, or an error writing the output, stops the run instead of processing the remaining segments; output files are still closed properly, also when the writer itself crashed.
- `call` and `per-read` skipped a segment that failed to process and exited successfully, leaving a gap in the output. They still process all other segments, but exit with an error that says how many segments failed.
- `call`, `per-read` and `bam` no longer accumulate finished segments without limit behind a slow one (deep coverage, a repeat); at most four segments per thread wait to be written.
- Crashes and internal errors end with a link to open a pre-filled GitHub issue, with the error, Rastair's version and the region being processed. With `--verbose`, only crashes and internal errors show a backtrace.

## Version 2.2.0 (2026-08-24)

Highlights:

- The reported beta values might change when updating to rastair 2.2. Methylation beta is now calculated by taking into account both positions of a CpG, meaning only reads containing `TG`/`CA` are counted as methylated.
- Reporting insertion and deletion calls.
  This is not enabled by default, while we're refining the model. Enable with `--experimental-indels`.
- Support guessing the read orientation using `--guess-read-orientation`

Further changes:

- Filtering reads by multiple tags now means a read needs to have _all_ of the specified tags.
- Support CRAM input for `rastair bam legacy` rewrites.

## Version 2.1.1 (2026-04-15)

Fixes for mbias plots.

## Version 2.1 (2026-03-19)

Highlights:

- Running Rastair's calling model on GPU.
  Using `--gpu` gives a significant speedup and works cross-platform (tested with Vulkan on Linux and Metal on macOS).
  (Other optimizations also improve CPU-only performance.)
- A new subcommand, `rastair bam` to add methylation annotations to exiting BAM files.

Further changes:

- Support single-stranded reads
- Support filtering reads by tags
- Fix BED output sometimes reporting misleading beta values

## Version 2.0 (2026-02-05)

This is a complete rewrite of Rastair.
While supporting the same methylation calling output, the main new addition is **variant calling** (outputting VCF/BCF).
Rastair now uses a bundles ML model to produce accurate calls while still being very performant.
