# Cargo tasks for building and testing Rastair

This crate is follows the [xtask pattern](https://github.com/matklad/cargo-xtask)
to add custom tasks to the `cargo` command.

## Usage

```
cargo xtask <task>
```

## Future Tasks

- [ ] Make sure local tools are installed, e.g. using [cargo-run-bin](https://github.com/dustinblackman/cargo-run-bin)
  - [ ] mdbook
  - [ ] cargo llvm-cov
  - [ ] nextest
- [ ] Add watch tasks for development
