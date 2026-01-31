# Format converter

Rastair includes a `convert` subcommand
that allows converting between many of the supported file formats.

## Streaming

Rastair's `convert` command supports streaming input and output,
so that it can be used in a pipeline with other commands.
Since it only guesses the formats based on file names,
they need to be specified explicitly when streaming.

For example,
to convert the first entries from a VCF file to BED format
you can use the following command:

```bash
head -n1000 test.vcf | rastair convert -f vcf -F bed | less
```

## Possible conversions

Note that not all formats contain the same information.
Rastair can only convert to formats that contain the same or less information than the source.

That means, for a given input, it can convert to a format in the same group or a lower group in this list:

1. Internal format (`.mpk.lz4`)

   Rastair's internal data structures, serialized.
   Contains all metrics in their original structure.
   Not guaranteed to be stable across versions.
2. @VCF/@BCF (`.vcf`/`.vcf.gz`, `.bcf`)

   Contains basically all information, serialized in a standard way.
   Can be processed using @bcftools.
3. @BED (`.bed`) and per-read BED (`.reads.bed`)

   Contains the most important calling information.
   Can be processed using @bedtools (or anything else that can read tab-separated-values files).
