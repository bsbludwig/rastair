# Methylation calling in Rastair

A core feature of Rastair is the ability to call @methylation from @TAPS sequencing data.
In @TAPS sequences, a mehtylated position is a `CG` sequence
where the @OT `C` is read as `T` and the @OB `G` is read as `A`.

For example, given this reference and @pileup (`.` means reference base):

```text
Ref     A T C G C C T  Strand
Reads   . . . . . . .    +
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

- Only @CpG sites are considered for methylation calling.
- @OT reads with `T` on the ref `C` position
- @OB reads with `A` on the ref `G` position
- Low number of other read bases at same positions
- Good read depth
- High base quality
- High mapping quality
- Coverage far away from start/end of reads
- Coverage not close to indels

## Filters

In @VCF, filters are used to indicate whether we have reliable evidence supporting a variant/methylation call.
We currently don't use filters in this implementation,
but the criteria above will become filters in the future.
