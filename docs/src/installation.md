# Installation

Right now, Rastair2 is not distributed as a pre-built binary.
You can build it from source from this repository.

## Building from source

Install the following dependencies:
- [Rust](https://www.rust-lang.org/tools/install) (version 1.88 or later)
- System dependencies for htslib and zlib: `libclang-dev cmake`

Clone the repository and build the project:

```bash
cargo build --release
```

The binary will be located in `target/release/rastair2`.
