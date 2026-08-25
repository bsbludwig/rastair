################################################################################
# BUILD STAGE
################################################################################
FROM debian:bookworm-slim AS builder

WORKDIR /app/rastair

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl git build-essential pkg-config unzip libclang-dev cmake && apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy source code
COPY . /app/rastair

# `--default-toolchain none` makes rustup use rust-toolchain.toml
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- --default-toolchain none --profile minimal --no-modify-path -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Build for release
RUN cargo xtask release

################################################################################
# RELEASE STAGE
################################################################################
FROM rocker/r-ver:4.3.3 AS release

RUN apt-get update && apt-get install -y --no-install-recommends \
        bash-completion bcftools bedtools gzip pandoc procps samtools tabix vcftools \
        libbz2-dev libcurl4-openssl-dev liblzma-dev libssl-dev libuv1-dev libxml2-dev zlib1g-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# `Ncpus` only affects build speed: every R package below is compiled from
# source, since the pinned CRAN versions postdate the base image's snapshot.
ENV R_INSTALL_OPTS="options(Ncpus = parallel::detectCores())"

# Install BiocManager for Bioconductor package management
RUN Rscript -e "$R_INSTALL_OPTS" -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

RUN Rscript -e "$R_INSTALL_OPTS" -e "BiocManager::install(c('Rsamtools', 'Biostrings', 'GenomicRanges'), version = '3.18', ask = FALSE)"

RUN Rscript -e "$R_INSTALL_OPTS" \
    -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('argparser', version = '0.7.2', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('knitr', version = '1.50', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('rmarkdown', version = '2.29', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('R.utils', version = '2.13.0', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('data.table', version = '1.17.8', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('ggplot2', version = '4.0.2', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('gtable', version = '0.3.6', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('ggside', version = '0.4.1', repos = 'https://cloud.r-project.org')"

RUN Rscript -e "invisible(lapply(c('argparser', 'Biostrings', 'data.table', 'GenomicRanges', 'ggplot2', 'ggside', 'gtable', 'knitr', 'R.utils', 'rmarkdown', 'Rsamtools'), function(p) library(p, character.only = TRUE)))"

# Copy the compiled binary from the build stage
# (up until here both stages can run in parallel)
COPY --from=builder /app/rastair/target/release/rastair /usr/local/bin/rastair

# Generate bash completion script in case people want to use it
RUN rastair internal shell-completions bash > /usr/share/bash-completion/completions/rastair
RUN echo "if ! shopt -oq posix; then if [ -f /usr/share/bash-completion/bash_completion ]; then . /usr/share/bash-completion/bash_completion; elif [ -f /etc/bash_completion ]; then . /etc/bash_completion; fi; fi" >> /root/.bashrc

# Copy R scripts
RUN mkdir -p /usr/local/share/rastair/scripts
COPY ./scripts /usr/local/share/rastair/scripts

# Set working directory
WORKDIR /app

# Set runtime env
ENV R_SCRIPT_DIR=/usr/local/share/rastair/scripts

# Default command is rastair
CMD ["/bin/bash"]
