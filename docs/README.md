# User-facing documentation

This folder contains user-facing documentation for the project.
It is written in Markdown and built with [mdbook](https://rust-lang.github.io/mdBook/).

## CLI Documentation

The CLI documentation is generated using `cargo run -- generate-cli-docs docs/src/cli.md`.

## Version

`{{VERSION}}` in any chapter becomes `[package].version` from the root `Cargo.toml`
(the release tag without its `v`); the `mdbook-glossary` preprocessor substitutes it.
