# Cargo tasks for building and testing Rastair

This crate is follows the [xtask pattern](https://github.com/matklad/cargo-xtask)
to add custom tasks to the `cargo` command.

## Installation

To install all tools, run `cargo bin --install`.

## Usage

```
cargo xtask <task>
```

where `<task>` is one of the tasks below.
Run with `--help` or `<task> --help` to see all options.

## Tasks

### `pre-commit`

Run quick checks.

On first run, it also adds itself as a hook that is run by `git` before writing commits.
This prevents commits with failing tests that would lead to `CI` failures anyway.
(If you need an escape hatch, use `git commit --no-verify -m WIP`.)

### `test`

Run all tests.

If there are changes in snapshots, this will fail.

### `insta`

Run tests and also accept new snapshots.
Make sure to review the snapshots before commiting them.

### `docs`

Generate documentation with `mdbook`.

Run with `cargo xtask docs --serve` to have it start a webserver
that automatically updates when the documentation files are changed.
You can access it under [localhost:3000](http://localhost:3000/).

Make sure to run `cargo bin --install` to (locally) install `mdbook` and its extensions.
