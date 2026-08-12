#!/usr/bin/env bash
# Installs the third-party tools that `tests/bam_external_tools.rs` cross-checks
# rastair against, on Linux x86_64.
#
# samtools, tabix/bgzip and R come from apt; bismark and modkit are not packaged
# for Debian/Ubuntu at all, but both publish prebuilt release tarballs, which is
# far cheaper than the conda image in Dockerfile.ci.
#
# Prints the directory to prepend to PATH on stdout. Usage:
#   echo "$(.github/scripts/install-external-tools.sh)" >> "$GITHUB_PATH"
set -euo pipefail

MODKIT_VERSION="0.6.4"
# The Perl Bismark, deliberately not the `bismark-rust-v3.x` rewrite: that one
# fails with "not yet implemented in this build: paired-end extraction (input
# has PAIRED flag set); PE arrives in Phase C", and rastair's test BAM is paired.
# Revisit once the rewrite ships PE support.
BISMARK_VERSION="0.25.1"

prefix="${EXTERNAL_TOOLS_PREFIX:-${RUNNER_TEMP:-/tmp}/rastair-external-tools}"
mkdir -p "$prefix/bin"

log() { echo "$*" >&2; }

# Bismark is a set of Perl scripts, so this is the source tarball; it needs
# `perl` and `samtools` on PATH at runtime.
if [ ! -x "$prefix/bismark/bismark_methylation_extractor" ]; then
  log "Installing bismark ${BISMARK_VERSION}"
  mkdir -p "$prefix/bismark"
  curl -fsSL \
    "https://github.com/FelixKrueger/Bismark/archive/refs/tags/v${BISMARK_VERSION}.tar.gz" \
    | tar -xz --strip-components=1 -C "$prefix/bismark"
fi

if [ ! -x "$prefix/modkit/modkit" ]; then
  log "Installing modkit ${MODKIT_VERSION}"
  mkdir -p "$prefix/modkit"
  curl -fsSL \
    "https://github.com/nanoporetech/modkit/releases/download/v${MODKIT_VERSION}/modkit_v${MODKIT_VERSION}_u16_x86_64.tar.gz" \
    | tar -xz --strip-components=1 -C "$prefix/modkit"
fi

for tool in "$prefix"/bismark/*; do
  if [ -f "$tool" ] && [ -x "$tool" ]; then ln -sf "$tool" "$prefix/bin/$(basename "$tool")"; fi
done
ln -sf "$prefix/modkit/modkit" "$prefix/bin/modkit"

# `uvx modbedtools` is used by one test; uv has no apt package either.
if ! command -v uvx >/dev/null 2>&1 && [ ! -x "$prefix/bin/uvx" ]; then
  log "Installing uv"
  UV_INSTALL_DIR="$prefix/bin" INSTALLER_NO_MODIFY_PATH=1 \
    sh -c "curl -LsSf https://astral.sh/uv/install.sh | sh" >&2
fi

export PATH="$prefix/bin:$PATH"
log "--- installed ---"
log "$(modkit --version 2>&1 || echo 'modkit MISSING')"
log "$(bismark_methylation_extractor --version 2>&1 | head -1 || echo 'bismark MISSING')"
log "$(uvx --version 2>&1 || echo 'uvx MISSING')"

echo "$prefix/bin"
