use std::io::Write as _;

use clap::{CommandFactory as _, Parser as _};
use clio::ClioPath;
use color_eyre::eyre::{Context, Result};
use rastair::{
    BamRewriteArgs, CallParams, ConvertParams, MBiasParams, MpkViewParams, PerReadParams, call,
    call_reads, setup_logging,
};
use tracing::{debug, info, warn};

/// Use mimalloc as the global allocator, which proves to be faster than the
/// default system allocator in our benchmarks.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[derive(Debug, clap::Parser)]
struct Cli {
    #[command(subcommand)]
    command: Subcommand,

    /// Enable more logging
    ///
    /// You can also use the `RASTAIR_LOG` environment variable to configure
    /// logging in a more precise way. See the documentation of the
    /// `tracing-subscriber` library to learn more.
    #[arg(short, long, global = true)]
    verbose: bool,
}

/// Rastair -- detect genetic variants and methylated positions from short-read
/// sequencing data created using TET-Assisted Pyridine-Borane Sequencing.
///
/// See <https://docs.rastair.com/> for more information.
#[derive(Debug, clap::Subcommand)]
#[command(version)]
#[allow(clippy::large_enum_variant)]
enum Subcommand {
    /// Call methylated positions
    ///
    /// Process TAPS-sequenced BAM files and call methylated positions.
    ///
    /// If no output file is specified, the output is written to stdout. You can
    /// use `--vcf` and `--bed` to write to files instead.
    ///
    /// If using `-c` (`--cpgs-only`), all CpG positions in the reference as
    /// well as de-novo CpGs are written. Stdout will default to BED.
    ///
    /// Only variants that pass all filters are written by default. Use `--all`
    /// to get a full VCF file.
    Call(CallParams),
    /// Call methylation per-read
    ///
    /// This will produce a bed file that list the methylation status of all
    /// CpGs in every read that overlaps a CpG, plus some other metadata
    PerRead(PerReadParams),
    /// Add methylation information to BAM files
    ///
    /// This will rewrite a BAM file to add methylation information and change
    /// the methylated positions in the sequence to their original base.
    Bam(BamRewriteArgs),
    /// Convert between different file formats
    Convert(ConvertParams),
    /// View internal format as JSON lines
    View(MpkViewParams),
    /// Calculate conversion per base position in read
    ///
    /// This will produce a `mbias.html` file with information about conversion
    /// counts relative to read position.
    ///
    /// Please note that this is currently implemented as an R script. Unless
    /// you're using the official Docker image, you need to install R and the
    /// necessary packages yourself.
    Mbias(MBiasParams),
    #[command(hide = true)]
    Internal {
        /// Generate documentation files
        #[command(subcommand)]
        command: Generate,
    },
}

#[derive(Debug, clap::Subcommand)]
enum Generate {
    /// Write shell completions
    ShellCompletions {
        /// The shell to generate the completions for
        #[arg(value_enum)]
        shell: clap_complete_command::Shell,
    },
    /// Write CLI help as markdown file
    #[command(hide = true)]
    CliDocs {
        /// The output file to write the markdown to
        output: ClioPath,
    },
    /// Write VCF fields as markdown file
    #[command(hide = true)]
    VcfDocs {
        /// The output file to write the markdown to
        output: ClioPath,
    },
}

fn main() -> Result<()> {
    let args = Cli::parse();
    setup_logging(args.verbose);

    #[cfg(unix)]
    // make sure we quit when the pipe closes
    // SAFETY: Calls libc, but at the start of the program
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }

    match args.command {
        Subcommand::Call(params) => {
            // track execution time
            let start = std::time::Instant::now();
            debug!(?params, "Running call command");
            call(params)?;
            let duration = start.elapsed();
            info!(?duration, "Call finished");
        }
        Subcommand::PerRead(params) => {
            // track execution time
            let start = std::time::Instant::now();
            debug!(?params, "Running call reads command");
            call_reads(&params)?;
            let duration = start.elapsed();
            info!(?duration, "Calling reads finished");
        }
        Subcommand::Bam(params) => {
            // track execution time
            let start = std::time::Instant::now();
            debug!(?params, "Running bam command");
            rastair::rewrite_bam(&params)?;
            let duration = start.elapsed();
            info!(?duration, "Bam rewrite finished");
        }
        Subcommand::Convert(params) => {
            // track execution time
            let start = std::time::Instant::now();
            debug!(?params, "Running convert command");
            rastair::convert(&params)?;
            let duration = start.elapsed();
            info!(?duration, "Convert finished");
        }
        Subcommand::View(params) => {
            warn!("This format is for internal use only and may change without notice.");
            rastair::mpk_view(&params)?;
        }
        Subcommand::Mbias(params) => {
            debug!(?params, "Running mbias command");
            rastair::mbias(&params)?;
        }
        Subcommand::Internal { command } => match command {
            Generate::ShellCompletions { shell } => {
                shell.generate(&mut Cli::command(), &mut std::io::stdout());
            }
            Generate::CliDocs { output } => {
                let mut file = output.clone().create().wrap_err("Failed to create output")?;
                let markdown = clap_markdown::help_markdown::<Cli>();
                file.write_all(markdown.as_bytes())
                    .wrap_err_with(|| format!("Failed to write CLI help to {output}"))?;
            }
            Generate::VcfDocs { output } => {
                let mut file = output.clone().create().wrap_err("Failed to create output")?;
                rastair::VcfRecord::description()
                    .to_markdown(&mut file)
                    .wrap_err("Failed to generate VCF docs")?;
            }
        },
    }

    Ok(())
}
