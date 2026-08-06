# Rastair 2

Rastair is a CLI application that allows
the simultaneous detection of genetic variants and methylated positions
from short-read sequencing data created using
[TET-Assisted Pyridine-Borane Sequencing (TAPS)](https://www.watchmakergenomics.com/taps.html), [Illumina 5-Base sequencing](https://www.illumina.com/science/genomics-research/articles/5-base-solution.html) or
any other method that converts 5mC to T.

## Installation

### From Source

To build rastair, you'll need:

- Rust (version 1.88 or later)
- `libclang-dev` and `cmake`

Then you can build it with:

```bash
cargo xtask release
```

### Binary releases

We provide binary releases for Linux and Mac. Visit [www.rastair.com](https://www.rastair.com/installation.html) for details.

### Docker

You can install the pre-built docker image in the usual way. For example, to install version 2.2.0, you could do:

```bash
docker pull sbludwig/rastair:version-2.2.0
```

## Usage examples

### 1. Call genomic variants and methylated positions at all genomic positions

```bash
rastair call -r reference.fa.gz -o test.vcf.gz test.bam
```

### 2. Only call methylation at CpGs and variants that generate a CpG

```bash
rastair call -r reference.fa.gz --bed test.bed.gz test.bam
```

This is significantly faster than calling variants across all genomic loci, but still runs the ML-based variant caller.

### 2. Call methylation at CpGs without ML filtering (ultra-fast)

```bash
rastair call --no-ml -r reference.fa.gz --bed test.bed.gz test.bam
```

This will be very similar to the output from other methylation callers.

## Contributing

If you want to contribute, please read the [CONTRIBUTING.md](CONTRIBUTING.md) file.
To start coding, have a look at [the `xtask` documentation](tools/xtask/README.md).

## License

This software is made available under a "non-commercial-use" license. Please refer to LICENSE.txt for details.
