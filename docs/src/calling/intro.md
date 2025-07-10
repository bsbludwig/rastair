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
It will set the beta value in the `M5mC` format field for each site.

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
