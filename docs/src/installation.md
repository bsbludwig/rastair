# Installation

## Homebrew

For installing rastair locally on macOS and Linux,
you can use [homebrew](https://brew.sh/):

```sh
brew tap bsblabludwig/rastair https://bitbucket.org/bsblabludwig/homebrew-rastair.git
brew install rastair
```

## Docker

You can use the pre-built docker image in the usual way:

```bash
docker pull sbludwig/rastair:version-{{VERSION}}
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

- [rastair-v{{VERSION}}-x86_64-unknown-linux-gnu.tar.gz](https://s3.eu-west-2.amazonaws.com/com.rastair.releases/build/release-v{{VERSION}}/rastair-v{{VERSION}}-x86_64-unknown-linux-gnu.tar.gz) (64-bit, x86)

This binary is built on RHEL 8 and works on distributions with glibc 2.28 or newer.

> [!NOTE]
> The binary brings the libraries it needs with it; the system's C library is the only one it takes from your machine. If your system uses a GLIBC older than 2.28, then you will have to [compile from source](#building-from-source).

### macOS

For Apple users, we provide both Apple Silicon and Intel binaries:

- [rastair-v{{VERSION}}-aarch64-apple-darwin.zip](https://s3.eu-west-2.amazonaws.com/com.rastair.releases/build/release-v{{VERSION}}/rastair-v{{VERSION}}-aarch64-apple-darwin.zip) (ARM "Apple Silicon", for all "M" line processors)

- [rastair-v{{VERSION}}-x86_64-apple-darwin.zip](https://s3.eu-west-2.amazonaws.com/com.rastair.releases/build/release-v{{VERSION}}/rastair-v{{VERSION}}-x86_64-apple-darwin.zip) (64-bit, older x86 Intel Macs)

## Building from source

### Pre-requisites
To compile from source, you need a working [Rust installation](https://www.rust-lang.org/tools/install). You do not need to pick a version: the repository contains a `rust-toolchain.toml`, so `rustup` downloads and uses the exact toolchain rastair is built and tested with. Building also needs a C compiler, since a few dependencies compile C code they bring along. Running the test suite additionally needs [cmake](https://cmake.org/download/): the tests read rastair's output back with [htslib](https://github.com/rust-bio/rust-htslib) as an independent check, and building htslib takes it. Rastair itself no longer uses htslib. On most systems, both are either already available or can be installed using a standard package manager:

First, clone the repo:

```bash
git clone git@github.com:bsbludwig/rastair.git
```

Then follow the platform-specific instructions below.

#### Ubuntu
```bash
sudo apt install clang cmake
```

#### Fedora
```bash
sudo dnf install -y clang cmake
```

#### macOS (Homebrew)
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

> [!TIP]
> On some systems, you might get performance improvements by allowing the compiler to use platform-specific optimisations:
> ```bash
> RUSTFLAGS="-C target-cpu=native" cargo xtask release
> ```
