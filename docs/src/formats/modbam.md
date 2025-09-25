# modBAM

The modBAM format is a @BAM/@SAM file
that contains modified base information in the `MM` and `ML` tags
as specified by the [Optional Fields Specification][samtags].
This format is useful for storing per-read modification information alongside the alignment data.
In Rastair, we're mainly interested in storing @CpG @methylation information, which is represented
as 5-methylcytosine (5mC) on the C and G on the opposite strand.

## Change of read sequence

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
when writing the modBAM output.

[samtags]: https://samtools.github.io/hts-specs/SAMtags.pdf
[modkit]: https://nanoporetech.github.io/modkit/
