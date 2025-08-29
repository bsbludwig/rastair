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

This software is made available under the terms of the [GNU Affero General Public License v3](https://www.gnu.org/licenses/agpl-3.0.html). If you require a more restrictive license for commercial purposes, please contact the authors to discuss alternative arrangements.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
