# Indel Calling

As of version 2.2.0, rastair supports calling @indel:pl variants in the `call` workflow. This is still an experimental feature.

There are two pathways, selected by mutually exclusive flags:

- `--experimental-indels` runs a fixed chain of hard filters with no @ML scoring.
- `--experimental-indels-ml` scores candidates with the indel @ML models instead.
  Combined with `--no-ml` it degrades to the hard-filter chain.

The default behaviour of `rastair call` remains SNP + methylation calling without indels.

> **The shipped indel models need a retrain.** `--experimental-indels-ml` warns about
> this at startup. Indel counting moved from alignment level to fragment level and the
> terminal-repeat feature was redefined (see below); the models learned their splits
> against the old values. The hard-filter pathway is unaffected.

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

## Fragment-Level Counting

Indel support is counted **per fragment**, not per alignment, matching the granularity of
the depth it is divided by. Overlapping mates of the same read pair therefore contribute a
single vote. Without this, an indel inside the mate-overlap window would have its VAF
roughly doubled relative to the depth denominator and be systematically over-called.

`--keep-overlapping-reads` disables the deduplication, and then everything — alt counts and
depth alike — stays at the alignment level, so the ratio remains internally consistent.

### Known gap: mate discordance

Deduplication picks the surviving mate from the anchor base and the `second` flag, blind to
indels. When the two mates' CIGARs disagree about an indel, the fragment's vote is decided
arbitrarily and the losing mate's evidence is dropped — measured at roughly 1,562 lost votes
per 5 Mb, 97% of it genuine mate-level alignment disagreement.

## Noise Exclusion

The hard-filter pathway excludes *noisy* fragments from genotyping: those with a soft-clip,
or with a tandem repeat at either read terminus (a homopolymer of ≥4 bp or a dinucleotide
repeat of ≥3 units). These are the alignments prone to slippage, and a read that slipped is
poor evidence in either direction.

Two properties matter:

- **The exclusion is symmetric.** Noisy fragments leave both the alternate count and the
  depth. Dropping them from the denominator alone raises VAF by roughly the noise rate; the
  binomial genotype flips from hom-ref to het at VAF ≈ 0.218, so a one-sided haircut can walk
  a low-VAF slippage artifact across that boundary on its own.
- **It does not look at the observed base**, so — unlike the mismatch filter — it cannot vary
  with local methylation.

The strand-bias test below is the one place all supporting fragments are used, noise
included: its null is the locus' full strand mix, so both sides have to be on the same footing.

The @ML pathway keeps its own, one-sided `depth_offset` so that feature distributions do not
move under the models.

## Calling Thresholds

Per-allele indel calls require:

- Minimum alternate observations (`--min-indel-ao`, default 2), counted after noise exclusion
- Minimum filtered depth (`--min-indel-depth`, default 2)

The binomial genotype (`--indel-error-rate`) then classifies the allele as hom-ref, het or
hom-alt. On the hard-filter pathway a hom-ref verdict is a rejection (`indel_hom_ref`); on
the ML pathway the genotype is informational and the @ML score decides.

## Strand Bias

The default strand gate is **both-strand concordance**: an allele supported on only one of
@OT/@OB is rejected (`indel_strand`). This is a *presence* rule — genuine TAPS indels are
frequently strand-asymmetric (alt-allele reference bias, not methylation), so the *degree* of
skew does not separate real calls from artifacts, whereas single-strandedness is strongly
enriched in them.

The split is on OT/OB, never the alignment reverse flag: both mates of a fragment share an
OT/OB assignment but have opposite reverse flags, so OT/OB is the only strand notion that
survives fragment deduplication.

### Cost at low AO

Concordance is deterministic, which makes it aggressive where support is thin. At
`--min-indel-ao 2` a genuine 2/0 split is a coin flip yet still fails, and an allele supported
entirely on one strand fails at any AO. Measured on chr12:1-6 Mb of an NA12878 TAPS BAM the
`ot > 0 && ob > 0` rule rejects 26.5% of all indel candidates, of which only ~15% are excess
over the chance rate; it also rejects homozygous indels at loci with no coverage on one strand.
To recover those, raise `--min-indel-ao` — so a one-sided split is no longer a coin flip —
rather than loosening the gate.

### Optional significance test

A two-sided exact binomial test of the OT/OB split against the strand mix of the *rest* of the
locus remains available as an additional gate via `--indel-strand-bias-alpha` (default 0 =
off). Taking the null from the non-supporting fragments means a locus whose coverage is itself
strand-skewed does not make every allele on it look biased. Because it keys on the *degree* of
skew — which on TAPS data also flags genuine strand-asymmetric indels — it is opt-in rather
than on by default.

## Output

Indels obey the same emission contract as SNVs:

- By default, only PASS indels are written.
- `--all` additionally writes failing alleles, each carrying the FILTER that rejected it
  (`indel_strand`, `indel_hom_ref`, `indel_no_ml`).
- `--cpgs-only` writes no indels at all.

`indel_no_ml` marks an allele on the ML pathway that the model could not score — a missing
model or a failed feature extraction. It is distinct from PASS on purpose: the `ML` format
field is empty in both cases, so without the FILTER "we could not judge this" would be
indistinguishable from "we judged it and it passed".

Verdicts travel on the call through `.mpk`, so `rastair convert` renders the same FILTERs as
a direct VCF run.
