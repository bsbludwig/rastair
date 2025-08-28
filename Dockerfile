################################################################################
# BUILD STAGE
################################################################################
FROM rust:1.89.0 AS builder

WORKDIR /app/rastair

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends curl unzip libclang-dev cmake && apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy source code
COPY . /app/rastair
# Build for release
RUN cargo xtask release

################################################################################
# RELEASE STAGE
################################################################################
FROM r-base:4.3.3 AS release
# Copy the compiled binary from the build stage
COPY --from=builder /app/rastair/target/release/rastair2 /usr/local/bin/rastair

# Install useful dependencies
RUN apt update && apt-get -y upgrade && apt-get -y --no-install-recommends install procps bash-completion && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN rastair internal shell-completions bash > /usr/share/bash-completion/completions/rastair

# Install BiocManager for Bioconductor package management
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

# Install required R and Bioconductor packages
RUN Rscript -e "BiocManager::install('Rsamtools', version = '3.18', ask = FALSE)" \
    -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('ggplot2', version = '3.5.1', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('gtable', version = '0.3.6', repos = 'https://cloud.r-project.org')"

# Set working directory
WORKDIR /app

# Copy R scripts
COPY ./Analysis /app/scripts
