# Changelog

This is the changelog for Rastair 2.

## Unreleased

- The default `--ml` threshold is 0.4 instead of 0.5, which recovers SNVs at low depth (+0.0024 SNV F1 at 10x) and is neutral at 26x and 73x.
- Pileups are read by [seqair](https://github.com/Softleif/seqair) instead of htslib.
  Along with some other changes, this should be another great performance improvement for Rastair.
- Removed `--linear-dedup-threshold`: it only tuned the htslib pileup path.
- `M5mC` (and `DPM5mC`, `ADM5mC`) are now written exactly where `CPG` or `CPGnovo` is set ([#12](https://github.com/bsbludwig/rastair/issues/12)).
- Reference-only records (`ALT=.`) are only written at CpG and de-novo CpG positions.
- Progress is shown in the terminal tab/taskbar (OSC 9;4) in terminals that support it, e.g. Windows Terminal, Ghostty, WezTerm, iTerm2 and Konsole.
- `call`, `per-read` and `bam` write their outputs (and indices) as `<name>.partial` and rename them only when the run succeeded, so a file under its final name is always complete.
- A crash (panic) in any thread, or an error writing the output, now stops the run instead of processing the remaining segments; output files are still closed properly, also when the writer itself crashed.
- `call` and `per-read` used to skip a segment that failed to process and exit successfully, leaving a gap in the output. They still process all other segments, but now exit with an error that says how many segments failed.
- `call`, `per-read` and `bam` no longer accumulate finished segments without limit behind a slow one (deep coverage, a repeat); at most four segments per thread wait to be written.
- Crashes and internal errors end with a link to open a pre-filled GitHub issue, with the error, Rastair's version and the region being processed. With `--verbose`, only crashes and internal errors show a backtrace.
- Fixed a hang when the GPU inference thread crashes: queued regions now fall back to CPU scoring, like after any other GPU failure.
- `call` links heterozygous variants that share read pairs into phase blocks, written as phased genotypes (`0|1`) with a new `PS` FORMAT field. This is always on, except for CpG-only output. Each read's vote comes from the pair-HMM's likelihoods, so a converted read at a CpG casts none.
- `call` calls indels with a TAPS-aware pair-HMM instead of the ML model. This is always on, except for CpG-only output.
- SNV genotypes (`GT`, `GC`) come from the same pair-HMM, which scores every alt the model can pass; the read-count genotyper and `--snv-genotype` are gone.
- With `--gpu`, the pair-HMM scores on the GPU as well; the calls are the same as on the CPU.
- An indel allele needs at least 5% of the column's depth, besides `--min-indel-ao`, to be scored.
- Indels at most 20 bp apart are called together, against haplotypes carrying every combination of them, with positions showing a nearby SNV left open; two indels on one chromosome are now both called, and an SNV next to an indel no longer turns into an indel call.
- Removed `--min-indel-depth`, `--indel-error-rate`, `--indel-max-mismatches` and `--indel-end-of-read-cutoff`.
- The `indel_strand` FILTER is replaced by `low_indel_qual`.
- The bundled model no longer contains indel models.
- Removed `--guess-read-orientation` from `call` and `per-read`: strands always come from the SAM flags.
- Removed `--keep-overlapping-reads`: overlapping mates always count once per fragment.
- Removed `--rescue-soft-clip-cpg`.
- The `GL` FORMAT field is gone: on variant calls it repeated `GC`, on hom-ref records it was a second number from the same binomial test. `GC` is declared `Number=1`, matching the one value it carries. The BED `gt_p_score` column now repeats `gt_conf_score`. `--vcf-format-fields GL` is an error.
- `QUAL` is the confidence of whatever decided the record: a reference-only record carries its `GC` instead of a constant 99, and with `--no-ml` an SNV record carries the pair-HMM's margin instead of a number derived from the error rate and depth.
- Removed `--m-max-coverage` and its `m_highDp` FILTER: `--max-coverage` caps the reads first.

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
