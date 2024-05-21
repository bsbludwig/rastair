################################################################################
# BUILD STAGE
################################################################################
FROM rust:1.78-buster AS builder

RUN cargo new --bin /app/rastair
WORKDIR /app/rastair

# Cache dependencies, only recompile when Cargo.toml or Cargo.lock changes
COPY ./Cargo.lock ./Cargo.lock
COPY ./Cargo.toml ./Cargo.toml
RUN cargo build --release
RUN rm ./src/*.rs
RUN rm ./target/release/deps/rastair*

# Copy source code
COPY ./src ./src
# Build for release
RUN cargo build --release

################################################################################
# RELEASE STAGE
################################################################################
FROM rust:1.78-slim-buster as release
# Copy the compiled binary from the build stage
COPY --from=builder /app/rastair/target/release/rastair /usr/local/bin/rastair

# Install `ps` command for process monitoring
RUN apt update
RUN apt-get -y upgrade
RUN apt-get -y install procps
