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
    io::mpk::viewer::MpkViewParams,
};
use tracing::{debug, info, warn};

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
    /// View internal format as JSON lines
    View(MpkViewParams),
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
    /// Write VCF fields as markdown file
    #[command(hide = true)]
    GenerateVcfDocs {
        /// The output file to write the markdown to
        #[arg()]
        output: ClioPath,
    },
}

fn main() -> Result<()> {
    color_eyre::install()
        .wrap_err("Failed to set up panic handler")
        .note("Seeing this error message is somewhat ironic, we know")?;
    reset_sigpipe();

    let args = Cli::parse();
    rastair2::utils::setup_tracing(args.verbose);

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
        Subcommand::View(params) => {
            warn!("This format is for internal use only and may change without notice.");
            rastair2::io::mpk::viewer::view(&params)?;
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
        Subcommand::GenerateVcfDocs { output } => {
            let mut file = output.clone().create().wrap_err("Failed to create output")?;
            rastair2::vcf::Record::description()
                .to_markdown(&mut file)
                .wrap_err("Failed to generate VCF docs")?;
        }
    }

    Ok(())
}

/*
 * some super-hacky sh*t to make this behave like a normal unix program and quit when the pipe ends
 */
fn reset_sigpipe() {
    #[cfg(unix)]
    // SAFETY: Calls libc
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
}
