# Architecture and Code Style

## Structure

This project is a Rust CLI application.
It is structured as a library with a binary crate that serves as the entry point.

## Testing

NOTE: `cargo clippy` is used to check for common mistakes and code style issues beyond default compiler warnings.

Run the tests with `cargo test`.
- Unit tests are added as test modules in the same file as the code they test.
- Integration tests for the CLI are in the `tests` directory.

### Tools

- [`proptest`](https://github.com/proptest-rs/proptest) is used for property-based testing.
- [`cargo insta`](https://insta.rs/) is used for snapshot testing in various places.

### Test Coverage

The code is tested using `cargo tarpaulin`:

```bash
cargo tarpaulin -o html --output-dir tmp/coverage
```

## Code Style

Usage of `rustfmt` is required.
Best set up formatting code on save.
Alternatively, run `cargo fmt` before committing.
