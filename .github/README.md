# CI/CD on GitHub Actions

| Workflow | Trigger | What it does |
| --- | --- | --- |
| [`ci.yml`](workflows/ci.yml) | push to `main`, any pull request | `toolchain-pins`, `fmt`, `clippy`, `test`, `external-tools` (all parallel) |
| [`release.yml`](workflows/release.yml) | push of a `v*` tag | runs CI, builds+signs Linux/macOS binaries, uploads to S3 + GitHub Releases, rebuilds and deploys the docs site |

`release.yml` can also be started manually from the Actions tab (**Run workflow**)
with a tag name. In that mode it builds and notarizes everything but **skips**
every upload, which is the safe way to test changes to the release pipeline.

## CI jobs

* **fmt** — `cargo fmt --all --check`. On failure the diff goes to the log and
  each offending file gets a file-level annotation (rustfmt has no line-level
  machine-readable check output).
* **clippy** — `cargo clippy --locked --workspace --all-targets`. The lint set
  lives in `[lints.clippy]` in `Cargo.toml`. Deliberately no `-D warnings`: the
  `deny`-level lints there are clean and do fail the job, but the `warn`-level
  ones still fire ~170 times, so they surface as annotations to burn down.
  Add `-- -D warnings` once they are gone. Diagnostics become inline
  annotations via [`problem-matchers/rust.json`](problem-matchers/rust.json), a
  plain regex config registered with `::add-matcher::` — no third-party action.
* **test** — `cargo-nextest nextest run --profile ci` plus `cargo test --doc`
  (nextest does not run doctests). The `ci` profile in `.config/nextest.toml`
  writes `target/nextest/ci/junit.xml`, which is turned into per-test
  annotations. nextest is invoked as `cargo-nextest nextest …`, **not**
  `cargo nextest …`, because `.cargo/config.toml` aliases the latter to
  cargo-run-bin, which would recompile nextest from source on every run.
  Without `--features external-tool-tests`, the two test targets that shell
  out to third-party tools are not built here.
* **external-tools** — the same runner (no Docker), running only
  `tests/bam_external_tools.rs` and `tests/mbias_report.rs`, which cross-check
  rastair against samtools, bismark, modkit, modbedtools and the R QC report.
  Kept separate so third-party CLI drift or a flaky download does not redden
  the main test signal; it runs in parallel with `test`.
  * It reuses the `test` job's `shared-key: ci` cache with `save-if: false`,
    so the dependency tree comes back warm without a second cache entry.
  * `--test <name>` rather than a nextest filter expression, because it also
    stops cargo compiling the test binaries the `test` job already covers.
  * `.github/scripts/install-external-tools.sh` (runnable locally) fetches
    bismark, modkit and uv; neither bismark nor modkit is packaged for
    Debian/Ubuntu. ~19 MB, ~4 s — not worth caching. R and its `r-cran-*`
    packages come from apt (~59 s), which is what lets `mbias_report.rs`
    render the report instead of self-skipping; drop them if that minute is
    not worth it.
  * Two tool-version traps are pinned in that script: the **Perl** Bismark
    (`v0.25.x`), because the `bismark-rust-v3.x` rewrite cannot do paired-end
    extraction yet and the test BAM is paired; and `--sampling-frac 1` for
    `modkit summary`, which dropped `--no-sampling`.
  * `Dockerfile.ci` still exists for running these tests on macOS, which has no
    modkit or Bismark build, but **CI does not build it** — this job installs
    onto the runner instead so it can share the cargo cache. That image has no
    R (so `mbias_report` self-skips there) and takes its tools from bioconda,
    so it is a fallback rather than a mirror; this job is the source of truth.

## Conventions

**Toolchain.** [`rust-toolchain.toml`](../rust-toolchain.toml) is the single
source of truth for the Rust version. No workflow and no Dockerfile pins one:
they all install rustup with `--default-toolchain none` and let that file decide.
[`scripts/check-toolchain-pins.sh`](scripts/check-toolchain-pins.sh) enforces
this (run it locally; the `toolchain-pins` job runs it in CI) and also checks
that `Cargo.toml`'s `rust-version` matches. That is why there is a local
[`actions/setup-rust`](actions/setup-rust/action.yml) instead of
`dtolnay/rust-toolchain`: that action takes the version from its own `@rev` or a
`toolchain:` input, so the version would have to be duplicated in every workflow,
and it assumes rustup already exists on the runner — which is false inside the
`redhat/ubi8:8.8` release container. The composite action covers both cases.

**Action pinning.** Every third-party action is pinned to a full commit SHA with
the human-readable tag in a trailing comment. Tags are mutable; SHAs are not, and
`softprops/action-gh-release` in particular runs with `contents: write`. When
bumping, resolve the new SHA with
`gh api repos/<owner>/<repo>/git/ref/tags/<tag>`.

## Not automated yet

**Docker image publishing.** `docs/src/installation.md` tells users to
`docker pull sbludwig/rastair:version-X.Y.Z`, but nothing in CI builds or pushes
that image — it is still a manual `docker build -f Dockerfile . && docker push`.
Adding it to `release.yml` needs `DOCKERHUB_USERNAME` + `DOCKERHUB_TOKEN`
secrets and, ideally, `docker/build-push-action` with
`platforms: linux/amd64,linux/arm64`.

**Build provenance attestation.** `actions/attest-build-provenance` would make
GitHub sign a statement — "this exact binary, identified by its SHA-256, was
produced by this workflow file at this commit, on a GitHub-hosted runner" — and
publish it to a public transparency log. Anyone could then run
`gh attestation verify rastair-v2.2.0-x86_64-unknown-linux-gnu.tar.gz --repo <owner>/rastair`
and get a yes/no answer on whether that download really came from this pipeline.
It needs `id-token: write` + `attestations: write` and about five lines in the
publish job. It is complementary to the Apple notarization, which only covers
the macOS binaries and only proves the Developer ID signed them, not which
pipeline built them.

## Required configuration

None of this is needed for `ci.yml` — only for tagged releases.

### Repository secrets

Settings → Secrets and variables → Actions → **Secrets**

| Name | Used by | Notes |
| --- | --- | --- |
| `MACOS_CERT_B64` | macOS build | Base64 of the "Developer ID Application" `.p12` (`base64 -i cert.p12 \| pbcopy`) |
| `CERT_PASSWORD` | macOS build | Password of that `.p12`; also used as the temporary keychain password |
| `APPLE_ID` | macOS build | Apple ID of the notarization account |
| `APPLE_TEAM_ID` | macOS build | 10-character Apple Developer Team ID |
| `APPLE_NOTARIZATION_PASSWORD` | macOS build | App-specific password, **not** the Apple ID password |
| `MACOS_SIGN_IDENTITY` | macOS build | e.g. `Developer ID Application: … (TEAMID)` |
| `AWS_ACCESS_KEY_ID` | release upload | Key with write access to the release bucket |
| `AWS_SECRET_ACCESS_KEY` | release upload | |
| `AWS_ACCESS_KEY_ID_DOC` | docs deploy | Key with write access to the docs bucket + CloudFront |
| `AWS_SECRET_KEY_DOC` | docs deploy | Note the name: `..._SECRET_KEY_DOC`, not `..._SECRET_ACCESS_KEY_DOC` |

The names are identical to the ones used in the old Bitbucket Pipelines and
CircleCI setups, so the values can be copied across verbatim.

### Repository variables

Settings → Secrets and variables → Actions → **Variables**

| Name | Used by | Example |
| --- | --- | --- |
| `S3_BUCKET` | release upload | release artifact bucket |
| `AWS_REGION` | release upload | e.g. `eu-west-2` |
| `S3_BUCKET_REL` | docs | bucket the docs link to for downloads (usually same as `S3_BUCKET`) |
| `AWS_REGION_REL` | docs | region for those download URLs |
| `S3_BUCKET_DOC` | docs deploy | bucket serving the docs site |
| `AWS_REGION_DOC` | docs deploy | region of that bucket |
| `CLOUDFRONT_DISTRIBUTION_ID` | docs deploy | `E12QS42FTTJKFT` (was hardcoded in `bitbucket-pipelines.yml`) |
| `MSTDN_INSTANCE` | docs | Mastodon instance shown on the landing page |
| `MSTDN_ACCOUNT` | docs | Mastodon account shown on the landing page |
| `ALL_APPLE_CERTS` | macOS build | whitespace-separated URLs of extra Apple root/intermediate `.cer` files |

## Cutting a release

1. Bump the version everywhere listed under "Release version bump checklist" in
   [`AGENTS.md`](../AGENTS.md).
2. Push the commit, wait for CI to go green.
3. `git tag v2.2.1 && git push origin v2.2.1`.

The tag must match `[package].version` in `Cargo.toml` exactly (`v` + version) or
the workflow stops before building anything.
