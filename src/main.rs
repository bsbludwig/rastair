use clap::{CommandFactory as _, Parser as _};
use color_eyre::eyre::Result;
use rastair2::call::{CallParams, call};
use tracing::{debug, info};
use tracing_subscriber::{fmt::format::FmtSpan, layer::SubscriberExt as _};

#[derive(Debug, clap::Parser)]
struct Cli {
    #[clap(subcommand)]
    command: Subcommand,
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
    color_eyre::install()?;
    let subscriber = tracing_subscriber::Registry::default()
        .with(tracing_error::ErrorLayer::default())
        .with(tracing_subscriber::EnvFilter::from_default_env())
        .with(
            tracing_subscriber::fmt::Layer::default()
                .with_target(true)
                .with_span_events(FmtSpan::CLOSE)
                .with_writer(std::io::stderr),
        )
        .with(tracing_subscriber::EnvFilter::from_default_env());
    tracing::subscriber::set_global_default(subscriber)?;

    let args = Cli::parse();

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
