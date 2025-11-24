# Rastair 2

Rastair is a CLI application that allows
the simultaneous detection of genetic variants and methylated positions
from short-read sequencing data created using the TET-Assisted Pyridine-Borane Sequencing (TAPS) method.

This repository contains the code for the **in-progress** development of Rastair version 2.

## Status

The code is in the early stages of development.

## Installation

There are no official releases yet.

### From Source

To build rastair, you'll need:

- Rust (version 1.88 or later)
- `libclang-dev` and `cmake`

Then you can build it with:

```bash
cargo xtask release
```

## Contributing

If you want to contribute, please read the [CONTRIBUTING.md](CONTRIBUTING.md) file.
To start coding, have a look at [the `xtask` documentation](tools/xtask/README.md).

## License

This software is made available under a "non-commercial-use" license. Please refer to LICENSE.txt for details.
