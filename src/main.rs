use std::io::Write as _;

use clap::{CommandFactory as _, Parser as _};
use clio::ClioPath;
use color_eyre::{
    Section,
    eyre::{Context, Result},
};
use rastair2::{
    call::{CallParams, call},
    convert::ConvertParams,
};
use tracing::{debug, info, warn};
use tracing_subscriber::{EnvFilter, layer::SubscriberExt as _};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[derive(Debug, clap::Parser)]
struct Cli {
    #[clap(subcommand)]
    command: Subcommand,

    /// Enable more logging
    #[arg(short, long, global = true)]
    verbose: bool,
}

/// Rastair2
///
/// Process TAPS-sequenced BAM files for methylation calling
#[derive(Debug, clap::Subcommand)]
#[command(version, about)]
#[allow(clippy::large_enum_variant)]
enum Subcommand {
    /// Call methylated positions
    Call(CallParams),
    /// Convert between different file formats
    Convert(ConvertParams),
    /// Write shell completions
    #[command(hide = true)]
    GenerateShellCompletions {
        /// The shell to generate the completions for
        #[arg(value_enum)]
        shell: clap_complete_command::Shell,
    },
    /// Write CLI help as markdown file
    #[command(hide = true)]
    GenerateCliDocs {
        /// The output file to write the markdown to
        #[arg()]
        output: ClioPath,
    },
}

fn main() -> Result<()> {
    color_eyre::install()
        .wrap_err("Failed to set up panic handler")
        .note("Seeing this error message is somewhat ironic, we know")?;

    let args = Cli::parse();

    let subscriber = {
        let default_log_settings =
            if args.verbose { "info,rastair2=debug" } else { "warn,rastair2=info" };
        let mut env_filter = EnvFilter::new(default_log_settings);
        if let Ok(env) = std::env::var("RASTAIR_LOG") {
            for directive in env.split(',') {
                if directive.is_empty() {
                    continue;
                }
                match directive.parse() {
                    Ok(parsed_directive) => {
                        env_filter = env_filter.add_directive(parsed_directive);
                    }
                    Err(error) => {
                        warn!(%error, "Warning: Invalid log directive `{directive}`");
                    }
                }
            }
        }

        tracing_subscriber::Registry::default()
            .with(tracing_error::ErrorLayer::default())
            .with(env_filter)
            .with(
                tracing_subscriber::fmt::Layer::default()
                    .with_target(true)
                    .with_thread_names(args.verbose)
                    // .with_span_events(FmtSpan::CLOSE) // maybe enable with flag
                    .with_writer(std::io::stderr),
            )
    };
    tracing::subscriber::set_global_default(subscriber)?;

    match args.command {
        Subcommand::Call(params) => {
            // track execution time
            let start = std::time::Instant::now();
            debug!(?params, "Running call command");
            call(&params)?;
            let duration = start.elapsed();
            info!(?duration, "Call finished");
        }
        Subcommand::Convert(params) => {
            // track execution time
            let start = std::time::Instant::now();
            debug!(?params, "Running convert command");
            rastair2::convert::convert(&params)?;
            let duration = start.elapsed();
            info!(?duration, "Convert finished");
        }
        Subcommand::GenerateShellCompletions { shell } => {
            shell.generate(&mut Cli::command(), &mut std::io::stdout());
        }
        Subcommand::GenerateCliDocs { output } => {
            let mut file = output.clone().create().wrap_err("Failed to create output")?;
            let markdown = clap_markdown::help_markdown::<Cli>();
            file.write_all(markdown.as_bytes())
                .wrap_err_with(|| format!("Failed to write CLI help to {output}"))?;
        }
    }

    Ok(())
}
