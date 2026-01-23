<img src="img/logo_white.png" alt="Rastair logo" id="logo">
<!-- <img src="img/logo_black.png" alt="Rastair logo" class="logo"> -->

## Background
Rastair is a command-line tool that allows the simultaneous detection of genetic variants and methylated positions from short-read sequencing data that was generated using a "mod-C&rarr;T" method, such as @TAPS or Illumina's @5Base technology.

Traditional @bisulfite sequencing (BS-seq) converts *all* non-modified cytosine (C) to thymine (T). This results in reads that differ substantially from the reference and are thus harder to align. Coverting most C to T also reduces the available information for variant identification. While several tools have been developed to overcome this problem, genetic variant calls from BS-seq remain substantially worse than those derived from whole-genome sequencing data.

In contrast, mod-C&rarr;T methods only affect around 60M positions in the human genome, equivalent to only approx. 2% of all nucleotides. This leads to greatly improved sequencing quality, higher mapping rates, and better yield from low-input DNA. It also makes it possible to identify genetic variation - in addition to epigenetic changes - with much higher accuracy. **Rastair implements a fast and accurate algorithm to simultaneously provide such high-quality variant *and* methylation calls.**

## License

Rastair is free for academic and other non-commercial use, and the [code is available on bitbucket](https://www.bitbucket.org/{{BITBUCKET_REPO_FULL_NAME}}/). You can read the details of the license [here](https://www.bitbucket.org/{{BITBUCKET_REPO_FULL_NAME}}/src/main/LICENSE.txt).

```admonish info
For commercial entities that would like to use rastair, please contact [enquiries@innovation.ox.ac.uk](mailto:enquiries@innovation.ox.ac.uk?cc=benjamin.schuster-boeckler%40ludwig.ox.ac.uk&subject=Rastair%20%2F%20reference%2024811) quoting reference 24811.
```

## Getting started

Briefly, we provide pre-built binaries for [Linux (x86)](https://s3.{{S3_REGION}}.amazonaws.com/{{S3_BUCKET}}/build/release-{{VERSION}}/rastair-{{VERSION}}-x86_64-unknown-linux-gnu.tar.gz), [Mac (Apple Silicon)](https://s3.{{S3_REGION}}.amazonaws.com/{{S3_BUCKET}}/build/release-{{VERSION}}/rastair-{{VERSION}}-aarch64-apple-darwin.zip) and [Mac (Intel)](https://s3.{{S3_REGION}}.amazonaws.com/{{S3_BUCKET}}/build/release-{{VERSION}}/rastair-{{VERSION}}-x86_64-apple-darwin.zip). We also provide a [docker image](https://hub.docker.com/r/sbludwig/rastair). Conda integration is still work in progress, but will happen soon. For build instructions and more details, see the [installation page](installation.md).

## Usage

Call methylation at all CpG positions (including CpGs formed by SNPs) from a bam file and output as a [tabix-indexed](https://www.htslib.org/doc/tabix.html) [bed file](formats/bed):

```bash
rastair call --bed output.bed.gz -r reference.fasta.gz input.bam
```
```admonish tip
By default, rastair will use all available CPU cores. You can restrict this with `-@ 1`.
```

Rastair can also produce variant and methylation calls in [VCF format](formats/vcf.md):

```bash
rastair call --vcf output.vcf.gz -r reference.fasta.gz input.bam
```

For a more in-depth look at different use-cases of rastair with practical examples, see the [examples](examples.md) section. For an explanation of the output file formats, see [BED](formats/bed.md) and [VCF](formats/vcf.md) sections.

## Citing rastair
A publication for rastair is in progress. We will update this page with a reference to the biorxiv preprint shortly.