# Installation

## Homebrew

For installing rastair locally on macOS and Linux,
you can use [homebrew](https://brew.sh/):

```sh
brew tap bsblabludwig/rastair git@bitbucket.org:bsblabludwig/homebrew-rastair.git
brew install rastair
```

## Docker

You can use the pre-built docker image in the usual way:

```bash
docker pull sbludwig/rastair:version-{{DOCKER_VERSION}}
```

### Building using Docker

You can also build Rastair using Docker.
Ensure you have Docker installed and running on your system,
then, clone the repo and build the image:

```bash
docker build -t rastair .
```

This image is based on the R base image and includes all necessary dependencies
to also run the bundled R scripts.

## Conda

```bash
conda create -n rastair_env -c conda-forge -c bioconda rastair={{VERSION}}
conda activate rastair_env
```

## Pre-built binaries

You can also download pre-built binaries for a number of platforms.

### Linux

We provide a pre-built binary for Linux (x86):

- [rastair-{{VERSION}}-x86_64-unknown-linux-gnu.tar.gz](https://s3.{{S3_REGION}}.amazonaws.com/{{S3_BUCKET}}/build/release-{{VERSION}}/rastair-{{VERSION}}-x86_64-unknown-linux-gnu.tar.gz) (64-bit, x86)

This binary was built on Ubuntu 20.04 but should work on most recent distributions.

```admonish info
While rastair itself is hard-linked and therefore independent of system libraries, this is unfortunately not yet the case for htslib. If your system uses a GLIBC older than 2.30, then you will have to [compile from source](#building-from-source).
```

### Mac OSX

For Apple users, we provide both Apple Silicon and Intel binaries:

- [rastair-{{VERSION}}-aarch64-apple-darwin.zip](https://s3.{{S3_REGION}}.amazonaws.com/{{S3_BUCKET}}/build/release-{{VERSION}}/rastair-{{VERSION}}-aarch64-apple-darwin.zip) (ARM "Apple Silicon", for all "M" line processors)

- [rastair-{{VERSION}}-x86_64-apple-darwin.zip](https://s3.{{S3_REGION}}.amazonaws.com/{{S3_BUCKET}}/build/release-{{VERSION}}/rastair-{{VERSION}}-x86_64-apple-darwin.zip) (64-bit, older x86 Intel Macs)

## Building from source

### Pre-requisites
To compile from source, you need a working [Rust installation](https://www.rust-lang.org/tools/install) (version 1.88 or later). Rastair depends on [rust-htslib](https://github.com/rust-bio/rust-htslib), which currently requires a working [clang library](https://clang.llvm.org/get_started.html) as well as [cmake](https://cmake.org/download/). On most systems, these are either already available or can be installed using a standard package manager:

First, clone the repo:

```bash
git clone git@bitbucket.org:{{BITBUCKET_REPO_FULL_NAME}}.git
```

Then follow the platform-specific instructions below.

#### Ubuntu
```bash
sudo apt install cmake
```

#### Fedora
```bash
sudo dnf install -y clang cmake
```

#### Mac OSX (Homebrew)
We assume that you have Xcode developer tools installed. In that case, you only need

```bash
brew install cmake
```

### Compile
Finally, build the project using:

```bash
cargo xtask release
```

The binary will be located in `target/release/rastair`.

```admonish tip
On some systems, you might get performance improvements by allowing the compiler to use platform-specific optimisations:
> ```bash
> RUSTFLAGS="-C target-cpu=native" cargo xtask release
> ```
```
