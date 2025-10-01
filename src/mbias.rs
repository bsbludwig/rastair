use crate::utils::cli;
use clap::value_parser;
use clio::ClioPath;
use color_eyre::{
    Result, Section as _,
    eyre::{Context, ContextCompat, ensure, eyre},
};
use rastair_types::RegionString;
use std::{env, path::PathBuf, process::Command};
use tracing::{debug, info, warn};

#[derive(Debug, Clone, clap::Args)]
pub struct MBiasParams {
    /// Input per-read BED file (can be gzipped)
    #[arg(value_name="BED_FILE", value_parser=value_parser!(ClioPath).exists().is_file())]
    #[arg(help_heading = cli::sections::INPUT)]
    pub bed_file: ClioPath,

    /// Genomic region
    #[arg(long)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub region: Option<RegionString>,

    /// Include bitflag as integer
    #[arg(long = "include-flag", default_value_t = 3)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub include_flag: i32,

    /// Exclude bitflag as integer
    #[arg(long = "exclude-flag", default_value_t = 3852)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub exclude_flag: i32,

    /// Read length as integer
    #[arg(long = "read-length")]
    #[arg(help_heading = cli::sections::FILTER)]
    pub read_length: Option<i32>,

    /// Path to tabix executable
    #[arg(long = "tabix-path", default_value = "tabix")]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub tabix_path: String,

    /// Output path prefix
    #[arg(long = "output-prefix", default_value = ".")]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub output_prefix: String,

    /// Override directory to find R scripts
    ///
    /// When not set, tries to look for `$rastair_path/scripts` and `./scripts`
    #[arg(long = "r-script-dir", env = "R_SCRIPT_DIR")]
    pub r_script_dir: Option<ClioPath>,
}

impl MBiasParams {
    fn as_flags(&self) -> Vec<String> {
        let mut params = vec![
            "--include-flag".into(),
            self.include_flag.to_string(),
            "--exclude-flag".into(),
            self.exclude_flag.to_string(),
            "--tabix-path".into(),
            self.tabix_path.clone(),
            "--output-prefix".into(),
            self.output_prefix.clone(),
        ];
        if let Some(region) = &self.region {
            params.push("--region".into());
            params.push(region.to_string());
        }
        if let Some(read_length) = self.read_length {
            params.push("--read-length".into());
            params.push(read_length.to_string());
        }
        params
    }

    fn find_scripts_dir(&self) -> Result<PathBuf> {
        // User provided directory
        if let Some(dir) = &self.r_script_dir {
            let dir = dir.to_path_buf();
            ensure!(dir.exists(), "Given R script directory not found: {}", dir.display());
            return Ok(dir);
        }

        let paths_to_try = [
            // Check scripts next to this executable
            env::current_exe()
                .wrap_err("Could not find path of current executable")
                .and_then(|exe_path| {
                    exe_path
                        .parent()
                        .wrap_err("Could not determine parent directory of executable")
                        .map(|p| p.to_path_buf())
                })
                .map(|p| p.join("scripts"))
                .wrap_err("Failed to find scripts directory next to executable"),
            // Check current working directory ./scripts
            env::current_dir()
                .wrap_err("Could not get current working directory")
                .map(|cwd| cwd.join("scripts"))
                .wrap_err("Failed to find scripts directory in current working directory"),
        ];

        for path in &paths_to_try {
            match path {
                Err(error) => {
                    // Log the error and continue
                    let error = format!("{error:#}");
                    warn!(error, "Error checking scripts directory");
                    continue;
                }
                Ok(path) => {
                    if path.exists() {
                        return Ok(path.clone());
                    } else {
                        debug!(?path, "Scripts directory not found");
                    }
                }
            }
        }

        Err(eyre!(
            "Could not find R scripts directory. Tried the following paths: {:?}",
            paths_to_try.iter().flatten().collect::<Vec<_>>()
        )).suggestion("You can set the R script directory using the --r-script-dir option or the R_SCRIPT_DIR environment variable.")
    }
}

// Call R script with right parameters
pub fn mbias(params: &MBiasParams) -> Result<()> {
    let scripts_dir = params.find_scripts_dir().wrap_err("Failed to find R scripts directory")?;
    let r_script = scripts_dir.join("mbias.R");
    ensure!(r_script.exists(), "mbias script not found in {scripts_dir:?}");

    ensure_tabix_index_exists(params).wrap_err("Failed to create tabix index file")?;

    let status = Command::new(r_script)
        .arg(params.bed_file.path())
        .args(params.as_flags())
        .status()
        .wrap_err("Failed to execute mbias script")?;

    ensure!(status.success(), "R script exited with {status}");
    Ok(())
}

fn ensure_tabix_index_exists(params: &MBiasParams) -> Result<()> {
    let tabix_index = params.bed_file.with_file_name({
        let mut name = params
            .bed_file
            .file_name()
            .wrap_err("Failed to get file name from BED file")
            .note("Make sure the BED file exists and is accessible")?
            .to_owned();
        name.push(".tbi");
        name
    });
    if !tabix_index.exists() {
        // try to create it
        let status = Command::new(&params.tabix_path)
            .arg("-p")
            .arg("bed")
            .arg(params.bed_file.path())
            .status()
            .wrap_err("Failed to execute `tabix` to create index")
            .with_note(|| {
                if params.tabix_path != "tabix" {
                    format!("Using executable at `{}`", params.tabix_path)
                } else {
                    String::new()
                }
            })
            .suggestion("Make sure tabix is installed and in your PATH")?;
        ensure!(status.success(), "Tabix exited with status: {status}");
        ensure!(
            tabix_index.exists(),
            "Tabix index file not found after creating it: {tabix_index:?}"
        );
        info!(path = ?tabix_index, "Created tabix index file");
    }
    Ok(())
}
