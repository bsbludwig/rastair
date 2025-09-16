use clap::{Parser, value_parser};
use clio::ClioPath;
use color_eyre::{
    Result,
    eyre::{Context as _, ensure},
};
use noodles::tabix;
use rastair::utils::logging::setup_logging;
use std::io::{Read as _, Write};
use tracing::{info, warn};

#[derive(Debug, Parser)]
struct Cli {
    #[clap(subcommand)]
    command: Subcommand,

    /// Enable more logging
    ///
    /// You can also use the `RASTAIR_LOG` environment variable to configure
    /// logging in a more precise way. See the documentation of the
    /// `tracing-subscriber` library to learn more.
    #[arg(short, long, global = true)]
    verbose: bool,
}

#[derive(Debug, clap::Subcommand)]
#[command(version)]
enum Subcommand {
    Query {
        /// Path to a bgzip-compressed and tabix-indexed file
        #[arg(value_parser=value_parser!(ClioPath).exists().is_file())]
        source: ClioPath,

        /// Region to query, e.g. `chr1:1000-2000`
        region: String,
    },
    Index {
        /// Path to a bgzip-compressed BED file
        #[arg(value_parser=value_parser!(ClioPath).exists().is_file())]
        source: ClioPath,
    },
}

fn main() -> Result<()> {
    let args = Cli::parse();
    setup_logging(args.verbose);

    match args.command {
        Subcommand::Query { source, region } => {
            let region = region.parse().wrap_err("Failed to parse region")?;

            let mut reader = tabix::io::indexed_reader::Builder::default()
                .build_from_path(source.path())
                .wrap_err("Failed to open tabix file")?;
            let query = reader.query(&region).wrap_err("Failed to query tabix file")?;

            let stdout = std::io::stdout();
            let mut stdout = stdout.lock();
            for result in query {
                match result {
                    Ok(record) => {
                        writeln!(stdout, "{}", record.as_ref())?;
                    }
                    Err(error) => {
                        warn!(%error, "Failed to read record from tabix file");
                    }
                }
            }
        }
        Subcommand::Index { source } => {
            check_file_is_bgz(&source).wrap_err_with(|| {
                format!("Source file `{}` is not bgzipped", source.path().display())
            })?;

            let tabix_index = source.path().with_extension("bed.gz.tbi");
            if tabix_index.exists() {
                warn!(path = ?tabix_index, "Tabix index already exists, not overwriting");
                return Ok(());
            }
            let status = std::process::Command::new("tabix")
                .arg("-p")
                .arg("bed")
                .arg(source.path())
                .status()
                .wrap_err("Failed to execute `tabix` to create index")?;
            ensure!(status.success(), "Tabix exited with status: {status}");
            ensure!(
                tabix_index.exists(),
                "Tabix index file not found after creating it: {tabix_index:?}"
            );
            info!(path = ?tabix_index, "Created tabix index file");
        }
    }

    Ok(())
}

fn check_file_is_bgz(source: &ClioPath) -> Result<()> {
    let mut file = std::fs::File::open(source.path()).wrap_err("Failed to open source file")?;
    let mut buffer = [0; 2];
    file.read_exact(&mut buffer).wrap_err("Failed to read from source file")?;
    ensure!(buffer == [0x1f, 0x8b], "Source file is not bgzipped (missing gzip magic number)");
    Ok(())
}
