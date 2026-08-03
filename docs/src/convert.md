# Format converter

Rastair includes a `convert` subcommand
that allows converting between many of the supported file formats.

## Streaming

Rastair's `convert` command supports streaming input and output,
so that it can be used in a pipeline with other commands.
Since it only guesses the formats based on file names,
they need to be specified explicitly when streaming.

For example,
to convert the first entries from a @VCF file to @BED format
you can use the following command:

```bash
head -n1000 test.vcf | rastair convert -f vcf -F bed | less
```

## Possible conversions

Note that not all formats contain the same information.
Rastair can only convert to formats that contain the same or less information than the source.
Concretely, these are:

- Convert between @VCF (incl. `.vcf.gz`) and @BCF
- Convert @VCF (incl. `.vcf.gz`) and @BCF to @BED
