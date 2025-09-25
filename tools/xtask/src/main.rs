use cargo_pgo::{
    get_cargo_ctx,
    pgo::{instrument::PgoInstrumentShortcutArgs, optimize::PgoOptimizeArgs},
};
use clap::Parser;
use color_eyre::{
    Result, Section,
    eyre::{Context, ContextCompat, bail, ensure, eyre},
};
use std::{
    env, fs,
    path::{Path, PathBuf},
    process::{Child, Command as StdCommand},
    thread,
};
use tracing::{debug, info};

#[derive(Debug, clap::Parser)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, clap::Subcommand)]
enum Command {
    /// Run quick checks, also adds itself as pre-commit hook
    PreCommit,
    /// Run tests
    Test {
        #[arg(long)]
        coverage: bool,
    },
    /// Accept new snapshots and run doctests
    Insta,
    /// Generate documentation with mdbook
    Docs {
        #[arg(long)]
        serve: bool,
    },
    /// Build a release version using cargo-pgo
    Release {
        /// Enable PGO (Profile-Guided Optimization) instrumentation
        #[arg(long, requires("args"))]
        pgo: bool,

        /// Additional arguments to pass to run the binary for profiling
        #[arg(last = true)]
        args: Vec<String>,
    },
}

fn main() -> Result<()> {
    color_eyre::install()?;
    tracing_subscriber::fmt().init();
    let workspace_root = workspace_dir().wrap_err("Could not get workspace root path")?;
    env::set_current_dir(&workspace_root)
        .wrap_err("Could not change current directory to workspace root")?;

    // add `cargo-bin` shims to PATH
    let shims_path = workspace_root.join(".bin/.shims");
    let path = std::env::var("PATH").unwrap_or_else(|_| String::new());
    // add custom shims
    let path = format!("{}:{path}", shims_path.display());
    // SAFETY: This is safe because we are at the beginning of the program
    unsafe {
        env::set_var("PATH", path);
    }

    let cli = Cli::parse();
    match cli.command {
        Command::PreCommit => {
            info!("Running pre-commit checks...");
            pre_commit().wrap_err("Pre-commit checks failed")?;
        }
        Command::Test { coverage } => {
            info!("Running tests...");
            let status = run_tests(coverage)?.wait().wrap_err("Failed to wait for test process")?;
            ensure!(status.success(), "Tests failed, please fix the issues before committing.");
        }
        Command::Insta => {
            info!("Running insta tests and doctests...");
            let status = StdCommand::new("cargo")
                .arg("insta")
                .arg("test")
                .arg("--workspace")
                .arg("--test-runner=nextest")
                .arg("--accept")
                .status()
                .wrap_err("Failed to run insta tests")?;
            ensure!(
                status.success(),
                "Insta tests failed, please fix the issues before committing."
            );
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

fn pre_commit() -> Result<()> {
    fn install_pre_commit_hook() -> Result<()> {
        let path = Path::new(".git/hooks/pre-commit");
        if path.exists() {
            debug!("Pre-commit hook already exists, skipping installation.");
            return Ok(());
        }
        let hook_content = "#! /usr/bin/env bash\nset -euo pipefail\ncargo xtask pre-commit\n";
        fs::write(path, hook_content).wrap_err("Failed to write pre-commit hook")?;

        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let mut perms = fs::metadata(path)?.permissions();
            perms.set_mode(0o755);
            fs::set_permissions(path, perms)
                .wrap_err("Failed to make pre-commit hook executable")?;
        }

        info!("Pre-commit hook installed successfully.");
        Ok(())
    }

    install_pre_commit_hook().wrap_err("Failed to install pre-commit hook")?;

    let mut test_handle = StdCommand::new("cargo")
        .arg("nextest")
        .arg("run")
        .arg("--all")
        .spawn()
        .wrap_err("Failed to run tests")?;

    let fmt_handle = thread::spawn(|| {
        let fmt = StdCommand::new("cargo")
            .arg("fmt")
            .arg("--all")
            .arg("--check")
            .status()
            .wrap_err("Failed to run cargo fmt --check")?;
        ensure!(fmt.success(), "Cargo fmt check failed");
        Ok(())
    });

    if let Err(error) =
        fmt_handle.join().map_err(|e| eyre!("Failed to join fmt thread: {e:?}")).and_then(|res| res)
    {
        test_handle.kill().wrap_err("Failed to kill test process")?;
        return Err(error).suggestion("Run `cargo fmt` to fix formatting issues.");
    }

    let status = test_handle.wait().wrap_err("Failed to wait for test process")?;
    ensure!(status.success(), "Tests failed, please fix the issues before committing.");

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

fn build_release() -> Result<()> {
    let status = StdCommand::new(env!("CARGO"))
        .arg("build")
        .arg("--release")
        .status()
        .wrap_err("Failed to call cargo")?;
    ensure!(status.success(), "Cargo build failed");
    Ok(())
}

fn run_tests(with_coverage: bool) -> Result<Child> {
    let ci = std::env::var("CI").is_ok();
    if with_coverage {
        ensure!(
            StdCommand::new("cargo").arg("llvm-cov").arg("--version").status()?.success(),
            "cargo-llvm-cov is not installed. Please install it with: cargo install cargo-llvm-cov"
        );
        info!("Running tests with coverage...");
        StdCommand::new("cargo")
            .arg("llvm-cov")
            .arg("test")
            .arg("--workspace")
            .arg("--doctests")
            .env("RUSTC_BOOTSTRAP", "1") // for doctests, don't worry about it
            .env("INSTA_UPDATE", if ci { "auto" } else { "always" })
            .spawn()
            .wrap_err("Failed to run tests with coverage")
    } else {
        StdCommand::new("cargo")
            .arg("test")
            .arg("--all")
            .env("INSTA_UPDATE", if ci { "auto" } else { "always" })
            .spawn()
            .wrap_err("Failed to run tests")
    }
}

fn generate_docs(serve: bool) -> Result<()> {
    StdCommand::new("mdbook")
        .arg("--version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .is_success()
        .wrap_err("Failed to run mdbook")
        .suggestion("You might need to install it with: cargo bin --install")?;

    let (sender, receiver) = std::sync::mpsc::channel();
    ctrlc::set_handler(move || {
        info!("Received Ctrl+C, shutting down...");
        let _ = sender.send(());
    })?;

    info!("Generating CLI docs...");
    StdCommand::new(env!("CARGO"))
        .arg("run")
        .arg("--quiet")
        .arg("--")
        .arg("internal")
        .arg("cli-docs")
        .arg("docs/src/cli.md")
        .is_success()
        .wrap_err("Failed to generate CLI args docs")?;

    info!("Generating format docs...");
    StdCommand::new(env!("CARGO"))
        .arg("run")
        .arg("--quiet")
        .arg("--")
        .arg("internal")
        .arg("vcf-docs")
        .arg("docs/src/formats/vcf-fields.md")
        .is_success()
        .wrap_err("Failed to generate VCF format args docs")?;

    let mut mdbook = StdCommand::new("mdbook");
    mdbook.current_dir("docs");
    mdbook.env("CARGO_TERM_QUIET", "true");

    info!("Running mdbook...");
    if serve {
        let mut child = mdbook.arg("serve").spawn().wrap_err("Failed to start mdbook serve")?;
        let _pls_exit = receiver.recv();
        let _ = child.kill();
    } else {
        mdbook.arg("build").is_success().wrap_err("Failed to build mdbook")?;
    }

    Ok(())
}

fn build_pgo_release(args: &[String]) -> Result<()> {
    let args = PgoInstrumentShortcutArgs::try_parse_from(
        ["rastair"].into_iter().chain(args.iter().map(|x| x.as_str())),
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

pub trait ExitStatusResultExt {
    fn is_success(&mut self) -> Result<()>;
}

impl ExitStatusResultExt for StdCommand {
    #[track_caller]
    fn is_success(&mut self) -> Result<()> {
        let status = self.status().wrap_err("Failed to run command")?;
        if !status.success() {
            bail!("Command failed with status: {}", status)
        }
        Ok(())
    }
}

impl ExitStatusResultExt for std::process::ExitStatus {
    #[track_caller]
    fn is_success(&mut self) -> Result<()> {
        if !self.success() {
            bail!("Command failed with status: {}", self)
        }
        Ok(())
    }
}
