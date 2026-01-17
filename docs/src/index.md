# Rastair

## Background
Rastair is a command-line tool that allows the simultaneous detection of genetic variants and methylated positions from short-read sequencing data that was generated using a "mod-C -> T" method, such as @TAPS or Illumina's @5Base technology.

Traditional @bisulfite sequencing (BS-seq) converts *all* non-modified cytosine (C) to thymine (T). This results in reads that differ substantially from the reference and are thus harder to align. Coverting most C to T also reduces the available information for variant identification. While several tools have been developed to overcome this problem, genetic variant calls from BS-seq remain substantially worse than those derived from whole-genome sequencing data.

In contrast, mod-C -> T methods only affect around 60M positions in the human genome, equivalent to only approx. 2% of all nucleotides. This leads to greatly improved sequencing quality, higher mapping rates, and better yield from low-input DNA. It also makes it possible to identify genetic variation - in addition to epigenetic changes - with much higher accuracy. **Rastair implements a fast and accurae algorithm to simultaneously provide such high-quality variant *and* methylation calls.**

## License

Rastair is free for academic and other non-commercial use, and the [code is available on bitbucket](https://www.bitbucket.org/bsblabludwig/rastair/). For commercial entities that would like to use rastair beyond initial evaluation, please contact [enquiries@innovation.ox.ac.uk](mailto:enquiries@innovation.ox.ac.uk?cc=benjamin.schuster-boeckler%40ludwig.ox.ac.uk&subject=Rastair%20%2F%20reference%2024811) quoting reference 24811. You can read the details of the license [here](https://bitbucket.org/bsblabludwig/rastair/src/main/LICENSE.txt).

## Getting started

Detailed installation instructions can be found [here](installation.md). Briefly, we provide pre-built binaries for [Linux (x86)]() and [Mac (Apple Silicon)](). You can also build rastair from source by cloning the [code repository](https://www.bitbucket.org/bsblabludwig/rastair/) and running

```bash
cargo build -r
```

We also provide a [docker image](https://hub.docker.com/r/sbludwig/rastair). Conda integration is still work in progress, but will happen soon.

## Usage

For a brief introduction to the main use-cases of rastair with practical examples, see the [examples](examples.md) section. For an explanation of the output file formats, see [BED](formats/bed.md) and [VCF](formats/vcf.md) sections.
