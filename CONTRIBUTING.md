# Contributing to Rastair 2

Code contributions are welcome, and can be submitted for review as pull requests. Please see below for some tips how to get started.

Note that as of version 2.0.0, rastair is distributed under a non-commerical use license (see LICENSE.txt). By submitting a pull request, you accept that your edits are re-released under the same terms, and you grant the University of Oxford (as owners of rastair) a non-exclusive, perpetual, irrevocable, worldwide and royalty-free license, with rights to sublicense through multiple levels of sublicensees, to use, copy, distribute, modify and create derivative works of your contribution.

See [ARCHITECURE.md](ARCHITECTURE.md) for the architecture of the project.

## Tools

Install tools with

```bash
cargo bin -s
cargo bin -i
```

## Code checks

Run `cargo xtask test` to run the tests.

Run `cargo xtask pre-commit` to set up a hook that runs formatting and tests before committing.
