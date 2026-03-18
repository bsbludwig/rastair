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

# Install useful dependencies
RUN apt update && apt-get -y upgrade && apt-get -y --no-install-recommends install procps bash-completion samtools bedtools bcftools vcftools tabix pandoc gzip libcurl4-openssl-dev libssl-dev && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install BiocManager for Bioconductor package management
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

# Install required R and Bioconductor packages
RUN Rscript -e "BiocManager::install(c('Rsamtools', 'Biostrings', 'GenomicRanges'), version = '3.18', ask = FALSE)"
RUN Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('argparser', version = '0.7.2', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('knitr', version = '1.50', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('rmarkdown', version = '2.29', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('R.utils', version = '2.13.0', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('data.table', version = '1.17.8', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('ggplot2', version = '4.0.2', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('gtable', version = '0.3.6', repos = 'https://cloud.r-project.org')" \
    -e "remotes::install_version('ggside', version = '0.4.1', repos = 'https://cloud.r-project.org')"

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
