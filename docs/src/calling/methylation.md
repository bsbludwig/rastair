# Methylation calling in Rastair

A core feature of Rastair is the ability to call @methylation from @TAPS sequencing data.
In @TAPS (and @5Base) sequences, a methylated position is a `CG` sequence
where the @OT `C` is read as `T` and the @OB `G` is read as `A`.

For example, given this reference and @pileup (`.` means reference base):

```text
Reference:     A T C G C C T  Strand
Reads:         . . . . . . .    +
               . . . A . .      -
               . . . A . . .    -
               . . T . . . .    +
               . . . . . . .    -
               . . T . . . .    +
               . . T . . . .    +
```

We can see that in the `CG` context, the `C` is read as `T` in 3 OT reads, and the `G` is read as `A` in 2 OB reads.
This gives us a good indication that the `C` is methylated.

## Output (VCF)

When methylation calling is enabled, Rastair will include all `CpG` sites in the output VCF file.
It will set the beta value(s) in the `M5mC` format field for each site.

In most cases, a single beta value is reported. However, when a position is both:

1. Part of an original CpG site in the reference genome, AND
2. Affected by a variant that creates a de-novo CpG site

Then **two beta values** will be reported in the `M5mC` field:

- The first value represents the methylation level of the original CpG context
- The second value represents the methylation level of the de-novo CpG context

This allows for accurate methylation quantification in complex scenarios where multiple CpG contexts overlap at a single position.

If a certain threshold of confidence is met,
the alt allele will be set to `.` since it is not considered a variant in the traditional sense.

## Criteria for methylation calling

- Only @CpG sites (including [de-novo CpGs](denovo.md)) are considered for methylation calling.
- @OT reads showing `TG` on a ref `CG` position
- @OB reads showing `CA` on the ref `CG` position
- Not a homozygous variant away from C/T
- Sufficient read depth (controlled with `--m-min-depth`)
- High base quality (controlled with `--min-baseq`)
- High mapping quality (controlled with `--min-mapq`)

### Both positions matter

Rastair only counts those positions as methylated
where both the position itself and the other side of the CpG
match the expected pattern.
That means, that an @OT read showing `TA` is not considered methylation evidence,
only `TG` is.
This also applies to de-novo CpGs.

This method of "pairwise" counting shows different beta values
compared to a trivial count looking at just one position.
But it can be seen as biologically more correct,
and the beta at heterozygous positions and the neighbour are more similar,
with $\beta_{snp} - \beta_{neighbor}$ approximating $0$.

## Beta value correction

When a CpG also overlaps a heterozygous variant, then only the non-variant C/G bases are available for methylation. Where the variant allele is a T/A at a C/G position, respectively, this thus requires that we do not count the "genetic T/A alleles", but only the chemically converted T/A. Unfortunately, we do not know which ones are which, so rastair currently uses a simple heuristic: we assume a diploid genome, and thus half of all reads are expected to derive from the T/A allele, while the other half represents the C/G. We therefore only consider the _excess_ of T/A reads beyond half the reads observed as bona-fide "modified":

$$
M_{excess} = \max \{ M_{total} - \frac{M + U}{2}, 0 \} \\
   \beta = \frac{M_{excess}}{M_{excess} + U}
$$

Where $M$ are the number of modified reads and $U$ are the number of unmodified reads.
