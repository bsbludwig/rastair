# Changelog

This is the changelog for Rastair 2.

## Unreleased

- `M5mC` (and `DPM5mC`, `ADM5mC`) are now written exactly where `CPG` or `CPGnovo` is set ([#12](https://github.com/bsbludwig/rastair/issues/12)).
- Reference-only records (`ALT=.`) are only written at CpG and de-novo CpG positions.

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
