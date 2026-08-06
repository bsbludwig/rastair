# Indel Calling

As of version 2.2.0, rastair supports calling @indel:pl variants in the `call` workflow. This is still an experimental feature that is enabled with the `--experimental-indels` command-line flag.

This mode extends normal variant/methylation processing with indel extraction from pileups, per-read filtering, and per-allele genotyping.

## How Indel Candidates Are Built

At each pileup position, Rastair collects insertion and deletion observations from aligned reads.
Each observed indel is represented as an allele and aggregated with supporting read-level metrics.

Read-level filters include:

- End-of-read cutoff (`--indel-end-of-read-cutoff`)
- Maximum non-TAPS mismatches (`--indel-max-mismatches`)
- Existing base-quality, mapping-quality, flag, and overlap filters

The mismatch filter is TAPS-aware:

- `C→T` mismatches on @OT reads are ignored
- `G→A` mismatches on @OB reads are ignored

These conversions are expected methylation signal and should not penalize indel-supporting reads.

## Calling Thresholds

Per-allele indel calls require:

- Minimum alternate observations (`--min-indel-ao`)
- Minimum filtered depth (`--min-indel-depth`)

When @ML is disabled (`--no-ml`), indels are hard-filtered by binomial genotyping.
When ML is enabled, alleles passing basic AO/depth thresholds are forwarded to ML scoring, while binomial genotype is still emitted for informational genotyping fields.

## Scope and Status

Indel calling is marked experimental and enabled only with `--experimental-indels`.
The default behavior of `rastair call` remains SNP + methylation calling without indels.
