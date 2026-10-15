# Indel Calling

As of version 2.2.0, rastair supports calling @indel:pl variants in the `call` workflow.

```admonish warning
This is still an experimental feature!
The `--experimental-indels` command-line flag must be used to enable it.
```

## How Indel Candidates Are Built

At each pileup position, Rastair collects insertion and deletion observations from aligned reads.
Each observed indel is represented as an allele and aggregated with supporting read-level metrics.

Read-level filters include:

- End-of-read cutoff (`--indel-end-of-read-cutoff`)
- Maximum non-TAPS mismatches (`--indel-max-mismatches`)
- Existing base-quality, mapping-quality, flag, and overlap filters

See the [CLI reference][cli-indels] for a full list.

[cli-indels]: ../cli.md#indel-options

The mismatch filter is TAPS-aware:

- `C→T` mismatches on @OT reads are ignored
- `G→A` mismatches on @OB reads are ignored

These conversions are expected methylation signal and should not penalize indel-supporting reads.

## Modes

A optional "mode" can be chosen when specifying the flag, `--experimental-indels=[ml|no-ml|ml-rescue]`.
By default, a @ML model is used to call indels, same as for @SNP:pl.

While this feature is experimental,
the calling behaviour can be fine-tuned by selecting a "mode".

| Mode           | Implementation                                               | Profile              |
| -------------- | ------------------------------------------------------------ | -------------------- |
| `no-ml`        | A fixed filter chain                                         | More false positives |
| `ml` (default) | Reject an allele based on model alone                        | More false negatives |
| `ml-rescue`    | Run both, reinstate an allele the filters dropped as hom-ref | Best recall          |

Which is best depends on depth.
Measured on GIAB HG001 across 20 libraries,
`ml-rescue` won on every library above ~32x and `ml` on every library below ~25x,
with the crossover confirmed within a single subsampled library.
Below ~10x the arms are not distinguishable from replicate noise.

```admonish note
These settings do not affect calling non-indel positions.
Running `rastair call --no-ml` will however mean that `ml` and `ml-rescue` fall back to `no-ml` behaviour.
```

## Molecular Genotyping and Canonical Positions

A model with `Haplotypes` feature metadata enables full PairHMM scoring and
molecular genotype assignment automatically. Simple indels are left-normalized
before feature extraction and calling. Reference context, repeat lengths,
per-molecule counts and boundary metrics use that canonical anchor; read
sequences retain their original CIGAR geometry for likelihood calculation.
A missing or low-quality normalized boundary is represented as uncertainty,
rather than dropping an otherwise accepted indel observation.

Reference support is local to that anchor. A noisy mate that does not cover the
position cannot remove its partner's usable reference observation. A rejected
indel is uninformative at its raw and normalized anchors; it cannot be counted
as clean reference support. Overlapping mates still count as one molecule.

The molecular model evaluates reference/alternate copy number jointly. Multiple
indel alternatives at an anchor can produce a `1/2` genotype, with their actual
sequences distinguished in likelihood scoring, priors and pruning. Unsafe
complex events and search limits use an explicit count-based fallback.

`GTM=pairhmm` identifies molecular genotypes and `GTM=counts` the count path.
`GQ` is genotype confidence; on ML-filtered calls, QUAL measures variant presence.
The forest does not supply additional independent genotype evidence.
In `ml-rescue` mode, a high forest score cannot force a molecular `0/0` decision
to become `0/1`. See [genotyping](genotyping.md) for confidence limitations and
fallback conditions.
