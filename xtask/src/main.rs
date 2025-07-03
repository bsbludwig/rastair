use cargo_pgo::{
    get_cargo_ctx,
    pgo::{instrument::PgoInstrumentShortcutArgs, optimize::PgoOptimizeArgs},
};
use clap::Parser;
use color_eyre::{
    Result,
    eyre::{Context, ContextCompat, ensure, eyre},
};
use std::{
    env::set_current_dir,
    path::{Path, PathBuf},
    process::Command as StdCommand,
};
use tracing::info;

#[derive(Debug, clap::Parser)]
struct Cli {
    #[clap(subcommand)]
    command: Command,
}

#[derive(Debug, clap::Subcommand)]
enum Command {
    /// Run tests
    Test {
        #[clap(long)]
        coverage: bool,
    },
    /// Generate documentation with mdbook
    Docs {
        #[clap(long)]
        serve: bool,
    },
    /// Build a release version using cargo-pgo
    Release {
        /// Enable PGO (Profile-Guided Optimization) instrumentation
        #[clap(long, requires("args"))]
        pgo: bool,

        /// Additional arguments to pass to run the binary for profiling
        #[clap(last = true)]
        args: Vec<String>,
    },
}

fn main() -> Result<()> {
    color_eyre::install()?;
    tracing_subscriber::fmt().init();
    let workspace_root = workspace_dir().wrap_err("Could not get workspace root path")?;
    set_current_dir(workspace_root)
        .wrap_err("Could not change current directory to workspace root")?;

    let cli = Cli::parse();
    match cli.command {
        Command::Test { coverage } => {
            info!("Running tests...");
            run_tests(coverage)?;
        }
        Command::Docs { serve: open } => {
            info!("Generating documentation...");
            generate_docs(open)?;
        }
        Command::Release { pgo, args } => {
            if pgo {
                info!("Building release version with PGO...");
                build_pgo_release(&args)?;
            } else {
                info!("Building release version without PGO...");
                build_release()?;
            }
        }
    }

    Ok(())
}

fn workspace_dir() -> Result<PathBuf> {
    let output = std::process::Command::new(env!("CARGO"))
        .arg("locate-project")
        .arg("--workspace")
        .arg("--message-format=plain")
        .output()
        .wrap_err("Could not call `cargo locate-project`")?
        .stdout;
    let cargo_path = Path::new(
        std::str::from_utf8(&output).wrap_err("Could not convert workspace path to string")?.trim(),
    );
    Ok(cargo_path.parent().wrap_err("could not get parent of Cargo.toml path")?.to_path_buf())
}

fn build_release() -> Result<(), color_eyre::eyre::Error> {
    let status = StdCommand::new(env!("CARGO"))
        .arg("build")
        .arg("--release")
        .status()
        .wrap_err("Failed to call cargo")?;
    ensure!(status.success(), "Cargo build failed");
    Ok(())
}

fn run_tests(with_coverage: bool) -> Result<()> {
    let ci = std::env::var("CI").is_ok();
    if with_coverage {
        ensure!(
            StdCommand::new("cargo-llvm-cov").arg("llvm-cov").arg("--version").status()?.success(),
            "cargo-llvm-cov is not installed. Please install it with: cargo install cargo-llvm-cov"
        );
        info!("Running tests with coverage...");
        StdCommand::new("cargo-llvm-cov")
            .arg("test")
            .arg("--workspace")
            .arg("--doctests")
            .env("RUST_BACKTRACE", "1")
            .env("RUSTC_BOOTSTRAP", "1") // for doctests, don't worry about it
            .env("INSTA_UPDATE", if ci { "auto" } else { "always" })
            .status()
            .wrap_err("Failed to run tests with coverage")?;
        return Ok(());
    } else {
        StdCommand::new("cargo")
            .arg("test")
            .arg("--all")
            .env("INSTA_UPDATE", if ci { "auto" } else { "always" })
            .env("RUST_BACKTRACE", "1")
            .status()
            .wrap_err("Failed to run tests")?;
    }

    Ok(())
}

fn generate_docs(serve: bool) -> Result<()> {
    ensure!(
        StdCommand::new("mdbook").arg("--version").status()?.success(),
        "mdbook is not installed. Please install it with: cargo install mdbook"
    );

    let (sender, receiver) = std::sync::mpsc::channel();
    ctrlc::set_handler(move || {
        info!("Received Ctrl+C, shutting down...");
        let _ = sender.send(());
    })?;

    let mut cli_docs = StdCommand::new(env!("CARGO"))
        .arg("run")
        .arg("--")
        .arg("generate-cli-docs")
        .arg("docs/src/cli.md")
        .spawn()
        .wrap_err("Failed to generate CLI docs")?;

    let mut child = StdCommand::new("mdbook")
        .arg(if serve { "serve" } else { "build" })
        .current_dir("docs")
        .spawn()
        .wrap_err("Failed to start mdbook serve")?;

    let _pls_exit = receiver.recv();
    let _ = child.kill();
    let _ = cli_docs.kill();

    Ok(())
}

fn build_pgo_release(args: &[String]) -> Result<()> {
    let args = PgoInstrumentShortcutArgs::try_parse_from(
        ["rastair2"].into_iter().chain(args.iter().map(|x| x.as_str())),
    )?;
    let ctx = get_cargo_ctx(&[]).map_err(|e| eyre!("{e}"))?;
    cargo_pgo::pgo::instrument::pgo_instrument(
        ctx,
        args.into_full_args(cargo_pgo::build::CargoCommand::Run),
    )
    .map_err(|e| eyre!("{e}"))?;
    let ctx = get_cargo_ctx(&[]).map_err(|e| eyre!("{e}"))?;
    cargo_pgo::pgo::optimize::pgo_optimize(ctx, PgoOptimizeArgs::parse_from(["build"]))
        .map_err(|e| eyre!("{e}"))?;
    Ok(())
}
