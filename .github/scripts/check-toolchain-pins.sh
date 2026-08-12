#!/usr/bin/env bash
# Guards the "one Rust version" rule.
#
# rust-toolchain.toml is the single source of truth. Everything that builds this
# project must either derive the version from it or be checked against it here.
# Run locally with: .github/scripts/check-toolchain-pins.sh
set -euo pipefail

cd "$(dirname "$0")/../.."

fail() {
  echo "::error::$*"
  echo "FAIL: $*" >&2
  status=1
}
status=0

channel="$(sed -n 's/^[[:space:]]*channel[[:space:]]*=[[:space:]]*"\(.*\)".*/\1/p' rust-toolchain.toml)"
if [ -z "$channel" ]; then
  echo "Could not read [toolchain].channel from rust-toolchain.toml" >&2
  exit 1
fi
echo "rust-toolchain.toml channel: $channel"

# 1. Cargo.toml's declared MSRV must equal the pinned toolchain. Anything lower
#    is unverifiable, because rust-toolchain.toml forces this exact version on
#    everyone who builds from source with rustup.
msrv="$(sed -n 's/^rust-version[[:space:]]*=[[:space:]]*"\(.*\)".*/\1/p' Cargo.toml)"
if [ "$msrv" != "$channel" ]; then
  fail "Cargo.toml rust-version is '$msrv' but rust-toolchain.toml pins '$channel'"
else
  echo "OK: Cargo.toml rust-version = $msrv"
fi

# 2. No Dockerfile may hardcode a Rust version. They all install rustup with
#    --default-toolchain none and let rust-toolchain.toml decide.
for dockerfile in Dockerfile Dockerfile.ci Dockerfile.release; do
  [ -f "$dockerfile" ] || continue
  if hits="$(grep -nE 'FROM[[:space:]]+rust:[0-9]|--default-toolchain[[:space:]]+[0-9]' "$dockerfile")"; then
    fail "$dockerfile hardcodes a Rust version, use --default-toolchain none: $hits"
  else
    echo "OK: $dockerfile pins no Rust version"
  fi
done

# 3. Same for the workflows.
for workflow in .github/workflows/*.yml; do
  if hits="$(grep -nE 'toolchain:[[:space:]]*[0-9]|rust-toolchain@[0-9]|image:[[:space:]]*rust:[0-9]' "$workflow")"; then
    fail "$workflow hardcodes a Rust version: $hits"
  fi
done
echo "OK: no workflow hardcodes a Rust version"

exit "$status"
