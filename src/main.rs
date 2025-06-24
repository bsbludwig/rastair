use clap::{CommandFactory as _, Parser as _};
use color_eyre::{
    Section,
    eyre::{Context, Result},
};
use rastair2::call::{CallParams, call};
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
enum Subcommand {
    /// Call methylated positions
    Call(CallParams),
    /// Write shell completions
    #[command(hide = true)]
    GenerateShellCompletions {
        /// The shell to generate the completions for
        #[arg(value_enum)]
        shell: clap_complete_command::Shell,
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
                        warn!(?error, "Warning: Invalid log directive `{directive}`");
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
        Subcommand::GenerateShellCompletions { shell } => {
            shell.generate(&mut Cli::command(), &mut std::io::stdout());
        }
    }

    Ok(())
}
