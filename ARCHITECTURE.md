# Architecture and Code Style

This document describes the architecture and code style of the `rastair` project.

## Structure

This project is a Rust CLI application.
It is not meant to be released or used as a library,
even though it is structured as a library for internal organization.

### Crates

Aside from the main `rastair` crate, we factored out some functionality into separate crates for better organization and reusability:

- `crates/rastair-types`: Library crate with common types.
- `crates/rastair-vcf`: Library crate for handling VCF files.
- `tools/*`: Internal tooling for developing rastair.

## Testing

NOTE: `cargo clippy` is used to check for common mistakes and code style issues beyond default compiler warnings.

Run the tests with `cargo xtask test`.
- Unit tests are added as test modules in the same file as the code they test.
- Integration tests for the CLI are in the `tests` directory.
- We use `insta` for snapshot testing in various places.
  You can run the tests with `cargo insta test` to check the snapshots,
  and `cargo insta review` to review changes to the snapshots.

### Tools

- [`proptest`](https://github.com/proptest-rs/proptest) is used for property-based testing.
- [`cargo insta`](https://insta.rs/) is used for snapshot testing in various places.

### Test Coverage

The code is tested using [`cargo llvm-cov`](https://github.com/taiki-e/cargo-llvm-cov):

```bash
cargo xtask test --coverage
```

This is also run in CI, making the tests fail with insufficient coverage (currently less than 70% of lines covered).

You can use `#[cfg_attr(coverage_nightly, coverage(off))]` to ignore certain functions or modules from coverage analysis, e.g. `Debug` or `Display` implementations that are only used in logs.

## Performance

Good reference for general topics:
[The Rust Performance Book](https://nnethercote.github.io/perf-book/build-configuration.html)

To measure performance, a representative dataset is needed, not just the small test files in this repository.
For the following, we assume you have a good "call" command as `$call`, e.g. `call test.bam -r test.fa.gz -o tmp/test.bcf`

- Use [samply](https://github.com/mstange/samply/) to quickly get profiling data: `cargo build --profile profiling && samply record $CARGO_TARGET_DIR/profiling/rastair $call`
- You can use [cargo-pgo](https://github.com/Kobzol/cargo-pgo) for building with profile-guided optimizations (PGO):
  `cargo xtask release --pgo -- $call`

### Reducing memory pressure

Rastair uses [mimalloc](https://github.com/microsoft/mimalloc) instead of the system allocator.

To reduce allocations, we use
- `SmallVec` for lists where we know the maximum number of elements is often small
- `SmolStr` for short strings or those that are often reused (note that these strings are immutable)

## Code Style

### Formatting

Usage of `rustfmt` is required.
Best set up formatting code on save.
Alternatively, run `cargo fmt` before committing.

### Strictly-typed data

We use `struct`s and "newtypes" to represent data
with as much information in the type system as possible.
E.g., we use a `RootMeanSquare` type instead of `f64`.

### Multiple `impl` blocks

Instead of putting all logic in one big `impl` block for type,
we add an `impl` block where the methods are most likely to be used.
This might seem confusing in some cases, but often leads to better reading locality.

### CLI argument composition

We use `clap` for CLI argument parsing.
Instead of defining a massive struct with all possible options,
we define smaller structs in the places where they are needed,
e.g. in the segmenter, the variant caller, and the methylation caller.
These structs are then added as fields with `#[command(flatten)]` to the struct for the subcommand.
This allows us to keep the code modular and maintainable,
without the need to convert types and arguments.
