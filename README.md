_`rastair` provides a set of tools to process [TAPS sequencing data](https://www.nature.com/articles/s41587-019-0041-2)._

# Installation

`rastair` is written in Rust. If the `cargo` command is available on your system, you can install it like this:

```bash
cargo install --git https://bitbucket.org/bsblabludwig/rastair.git
```

We are planning to release statically compiled binaries for different systems in the future.

# Available commands

`rastair` offers a number of sub-commands, each with their own set of options. You can get a list of available commands by typing

```bash
rastair --help
```

Currently, the following sub-commands are available:

```text
Fast and flexible extraction of methylation information from bam files

Usage: rastair [OPTIONS] [COMMAND]

Commands:
  call      Call methylation at CpG positions from a bam file
  per-read  Call methylation per read. This will produce a bed file that list the methylation status of all CpGs in every read that overlaps a CpG, plus some other metadata
  mbias     Calculate conversion per base position in read
  map-cpgs  print a map of all CpGs in a fasta file and their indices as a bed file
  bed2pat   Utility function to convert per-read bed files to PAT files compatible with wgbstools and UXMtools
  help      Print this message or the help of the given subcommand(s)

Options:
  -v, --verbosity...
  -h, --help          Print help
  -V, --version       Print version
```

In general, you can use `--help` on all `rastair` sub-commands to get detailed instructions on the available options.

# Notes
## Soft-clipping option syntax
The correction of "m-bias", ie the loss of conversion around read ends, is an important aspect of several rastair sub-commands (e.g. `call`, `bed2pat`). The command-line argument for this was inspired by [MethylDackel](https://github.com/dpryan79/MethylDackel).
However, we decided to simplify the argument and make it more consistent (in our perspective, at least). This means we only have one set of parameters `--nOT` and `--nOB`. Each of them takes a comma-separated list of 4 integers, denoting the number of bases from the start/end of read 1/2 that should be _ignored_. However, unlike MethylDackel, rastair accounts for read orientation, so the command like arguments always refer to the physical start/end position of the read, not the position in the alignment. To give an example: imagine the following read pair

```text
    000000000111111111122
    123456789012345678901
R1: CG--------TG--------->
R2:                         <AC-------------GC-
                             876543210987654321
														 111111111000000000
```

This "F1R2" read represents the OT (ie R1 is the OT, and R2 is the reverse complement of the OT). A parameter of `--nOT 0,5,0,5` will exclude the `A` at position 18 in R2, because it ocurrs within 5 bases from the end of R2 _in read coordinates_, not in reference coordinates.

## Performance considerations
Most commands in `rastair` can be run in multi-threaded mode using the `-@ <ncores>` parameter. However, the performance increase has diminishing returns, as the threads have to eventually synchronise for writing to the output. We have found that increasing the number of BAM read threads to 2 (`--read-threads 2`) in combination with 4 to 6 processing threads (`-@ 4`) seems to perform best on a high-io-throughput system.

# License

This software is made available under the terms of the [GNU Affero General Public License v3](https://www.gnu.org/licenses/agpl-3.0.html). If you require a more restrictive license for commercial purposes, please contact the authors to discuss alternative arrangements.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

