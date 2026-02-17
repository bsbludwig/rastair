# modBAM

```admonish note
This is a preview feature.
```

The modBAM format is a @BAM/@SAM file
that contains modified base information.
Writing modified BAM files is useful for storing per-read modification information alongside the alignment data.
Right now, Rastair only allows storing @CpG @methylation information,
which is represented as 5-methylcytosine (5mC) on the C and G on the opposite strand.

Rastair supports two annotation formats:
1. the "standard" format using `MM` tags, as specified in [The SAM tags reference][samtags] version 4.5.
1. the "legacy" format using `XM`/`XR`/`XG` tags as used by @dragen and @bismark, as described in the [Illumina docs][dragen-docs].

## Standard mode

Writes the `MM` and `ML` tags as specified by the [Optional Fields Specification][samtags].
This is the format expected by tools like [modkit].

### Change of read sequence

The read sequence in the modBAM file will differ from the @read sequence in the input BAM file.

Since Rastair deals with reads from @TAPS,
methylated `C`s are represented as `T`s in the read sequence.
However, both the `MM` tag specification and other tools like [modkit]
expect the fundamental base to be in the read sequence.
That means, for a methylated CpG,
the read sequence should contain a `C` on the forward strand and a `G` on the reverse strand,
with only the `MM` tag indicating the modification
instead of the presence of a `T` or `A` in the sequence.

To be compatible with this,
Rastair will rewrite the read sequence at methylated positions
when writing the standard modBAM output.

## Legacy mode

Writes `XR`/`XG`/`XM` tags as used by @dragen and @bismark.
The `XM` tag marks methylated positions as `Z`, unmethylated target-base
positions as `z`, and everything else as `.`.
(The read sequence is not rewritten.)

Use this mode when downstream tools expect the legacy Bismark-style tags.

[samtags]: https://samtools.github.io/hts-specs/SAMtags.pdf
[modkit]: https://nanoporetech.github.io/modkit/
[dragen-docs]: https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/MPipelineMethBAM_fDG.htm
