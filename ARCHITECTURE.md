# Architecture and Code Style

## Structure

This project is a Rust CLI application.
It is structured as a library with a binary crate that serves as the entry point.

## Testing

The code is tested using `cargo test`.
Additionally, `cargo clippy` is used to check for common mistakes and code style issues.

### Test Coverage

The code is tested using `cargo tarpaulin`:

```bash
cargo tarpaulin -o html --output-dir tmp/coverage
```

## Code Style

Usage of `rustfmt` is required.
Best set up formatting code on save.
Alternatively, run `cargo fmt` before committing.
