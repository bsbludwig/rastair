# modBAM

> [!NOTE]
> This is a preview feature.

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

## Haplotype tags

`call` links heterozygous variants into phase blocks, written as phased genotypes with a `PS` phase set in the @VCF.
When the calls given to `rastair bam` are a VCF or BCF file instead of a @BED file,
both modes also tag each @read with the haplotype it comes from,
with the tags [WhatsHap] writes and tools like IGV and [modkit] read:

| Tag  | Type | Description                                                             |
| ---- | ---- | ----------------------------------------------------------------------- |
| `HP` | `i`  | Haplotype, 1 or 2: 1 carries the allele left of the `|` in the VCF `GT` |
| `PS` | `i`  | Phase set the haplotype belongs to, the `PS` of the VCF records         |

```bash
rastair call -r reference.fa.gz test.bam --vcf test.vcf.gz
rastair bam standard -r reference.fa.gz test.bam test.vcf.gz -o test_mm.bam
```

The VCF has to be compressed and indexed (`.vcf.gz` or `.bcf`); `call` writes the index next to it.
A BED file holds no phase, so with one the reads keep whatever `HP` and `PS` they have.

A read pair is tagged as a whole, from the phased heterozygous variants either of its reads shows,
so a read that covers none of them takes its mate's tag (if the mate starts within 1,000 bases).
Each variant counts with the quality of the base that shows it, and the phase set the read pair agrees with most clearly gives the tag.
Read pairs that cover no phased variant, or support both haplotypes equally, get no tags, and lose any they came with.

A converted base is not taken for a variant:
at a `C>T` variant, a `T` on the strand where a methylated `C` reads as `T` says which allele the read has
only if the read's next base is not a `G` (and likewise an `A` at a `G>A` variant on the other strand, after a `C`).
At an @indel, the read's sequence after the variant has to match one haplotype exactly and not the other,
so a read that ends inside a repeat the indel sits in does not count.

[samtags]: https://samtools.github.io/hts-specs/SAMtags.pdf
[WhatsHap]: https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag
[modkit]: https://nanoporetech.github.io/modkit/
[dragen-docs]: https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/MPipelineMethBAM_fDG.htm
