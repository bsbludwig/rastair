use crate::utils::cli;
use clap::{ArgGroup, value_parser};
use clio::ClioPath;
use color_eyre::{
    Result, Section as _,
    eyre::{Context, ContextCompat, ensure, eyre},
};
use rastair_types::RegionString;
use std::{env, path::PathBuf, process::Command};
use tracing::{debug, info, warn};

#[derive(Debug, Clone, clap::Args)]
#[command(group(
    ArgGroup::new("input")
        .args(["bed_file", "bam_file"])
        .required(true)
))]
pub struct MBiasParams {
    /// Input per-read BED file (must be tabix indexed or indexable)
    #[arg(long = "bed", value_name="BED_FILE", value_parser=value_parser!(ClioPath).exists().is_file(), value_hint=clap::ValueHint::FilePath)]
    #[arg(help_heading = cli::sections::INPUT)]
    pub bed_file: Option<ClioPath>,

    /// Input BAM file
    #[arg(long = "bam", value_name="BAM_FILE", value_parser=value_parser!(ClioPath).exists().is_file(), value_hint=clap::ValueHint::FilePath)]
    #[arg(help_heading = cli::sections::INPUT)]
    pub bam_file: Option<ClioPath>,

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

    /// Reference FASTA file (required for V-bias and GC/CpG bias plots)
    #[arg(long, value_parser=value_parser!(ClioPath).exists().is_file(), value_hint=clap::ValueHint::FilePath)]
    #[arg(help_heading = cli::sections::INPUT)]
    pub reference: Option<ClioPath>,

    /// VCF file with methylation calls (required for GC/CpG bias plots)
    #[arg(long, value_parser=value_parser!(ClioPath).exists().is_file(), value_hint=clap::ValueHint::FilePath)]
    #[arg(help_heading = cli::sections::INPUT)]
    pub vcf: Option<ClioPath>,

    /// Do not generate V-bias plots (faster)
    #[arg(long = "no-vbias")]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub no_vbias: bool,

    /// Do not generate GC/CpG bias plots
    #[arg(long = "no-gc")]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub no_gc: bool,

    /// Path to tabix executable
    #[arg(long = "tabix-path", default_value = "tabix")]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub tabix_path: String,

    /// Path to bcftools executable
    #[arg(long = "bcftools-path", default_value = "bcftools")]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub bcftools_path: String,

    /// Path to rastair executable
    #[arg(long = "rastair-path", default_value = "rastair")]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub rastair_path: String,

    /// Number of threads to use
    #[arg(long, default_value_t = 1)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub threads: i32,

    /// Treat the input as inverted, i.e. mod=unmod and unmod=mod
    #[arg(long)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub wgbs: bool,

    /// Output path prefix
    #[arg(long = "output-prefix", default_value = ".")]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub output_prefix: String,

    /// Override directory to find R scripts
    ///
    /// When not set, tries to look for `$rastair_path/scripts` and `./scripts`
    #[arg(long = "r-script-dir", env = "R_SCRIPT_DIR", value_hint=clap::ValueHint::DirPath)]
    pub r_script_dir: Option<ClioPath>,
}

impl MBiasParams {
    fn as_flags(&self) -> Vec<String> {
        let rastair_path = self.current_rastair_executable();
        let mut params = vec![
            "--include-flag".into(),
            self.include_flag.to_string(),
            "--exclude-flag".into(),
            self.exclude_flag.to_string(),
            "--tabix-path".into(),
            self.tabix_path.clone(),
            "--bcftools-path".into(),
            self.bcftools_path.clone(),
            "--rastair-path".into(),
            rastair_path,
            "--output-prefix".into(),
            self.output_prefix.clone(),
            "--threads".into(),
            self.threads.to_string(),
        ];
        if let Some(region) = &self.region {
            params.push("--region".into());
            params.push(region.to_string());
        }
        if let Some(read_length) = self.read_length {
            params.push("--read-length".into());
            params.push(read_length.to_string());
        }
        if let Some(reference) = &self.reference {
            params.push("--reference".into());
            params.push(reference.path().to_string_lossy().into_owned());
        }
        if let Some(vcf) = &self.vcf {
            params.push("--vcf".into());
            params.push(vcf.path().to_string_lossy().into_owned());
        }
        if self.wgbs {
            params.push("--wgbs".into());
        }
        // auto-disable plots that require inputs the user hasn't provided
        if self.no_vbias || self.reference.is_none() {
            params.push("--no-vbias".into());
        }
        if self.no_gc || self.reference.is_none() || self.vcf.is_none() {
            params.push("--no-gc".into());
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

    fn current_rastair_executable(&self) -> String {
        if self.rastair_path != "rastair" {
            return self.rastair_path.clone();
        }

        match env::current_exe() {
            Ok(path) => path.to_string_lossy().into_owned(),
            Err(error) => {
                warn!(
                    ?error,
                    "Failed to determine current rastair executable path, falling back to default"
                );
                self.rastair_path.clone()
            }
        }
    }
}

// Call R script with right parameters
pub fn mbias(params: &MBiasParams) -> Result<()> {
    let scripts_dir = params.find_scripts_dir().wrap_err("Failed to find R scripts directory")?;
    let r_script = scripts_dir.join("mbias.R");
    ensure!(r_script.exists(), "mbias script not found in {scripts_dir:?}");

    ensure_tabix_index_exists(params).wrap_err("Failed to create tabix index file")?;

    let mut command = Command::new(r_script);
    if let Some(bed_file) = &params.bed_file {
        command.arg("--bed").arg(bed_file.path());
    }
    if let Some(bam_file) = &params.bam_file {
        command.arg("--bam").arg(bam_file.path());
    }

    let status =
        command.args(params.as_flags()).status().wrap_err("Failed to execute mbias script")?;

    ensure!(status.success(), "R script exited with {status}");
    Ok(())
}

fn ensure_tabix_index_exists(params: &MBiasParams) -> Result<()> {
    let Some(bed_file) = &params.bed_file else {
        return Ok(());
    };
    let tabix_index = bed_file.with_file_name({
        let mut name = bed_file
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
            .arg(bed_file.path())
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

#[cfg(test)]
mod tests {
    use super::MBiasParams;

    #[test]
    fn as_flags_includes_extended_script_options() {
        let params = MBiasParams {
            bed_file: None,
            bam_file: None,
            region: None,
            include_flag: 3,
            exclude_flag: 3852,
            read_length: None,
            reference: None,
            vcf: None,
            no_vbias: false,
            no_gc: false,
            tabix_path: "tabix".into(),
            bcftools_path: "bcftools".into(),
            rastair_path: "rastair".into(),
            threads: 4,
            wgbs: true,
            output_prefix: ".".into(),
            r_script_dir: None,
        };

        let flags = params.as_flags();

        assert!(flags.windows(2).any(|window| window == ["--bcftools-path", "bcftools"]));
        assert!(flags.windows(2).any(|window| window == ["--threads", "4"]));
        assert!(flags.iter().any(|flag| flag == "--wgbs"));
        assert!(flags.iter().any(|flag| flag == "--no-vbias"));
        assert!(flags.iter().any(|flag| flag == "--no-gc"));
        assert!(flags.windows(2).any(|window| window[0] == "--rastair-path"));
    }

    #[test]
    fn as_flags_preserves_explicit_rastair_path_override() {
        let params = MBiasParams {
            bed_file: None,
            bam_file: None,
            region: None,
            include_flag: 3,
            exclude_flag: 3852,
            read_length: None,
            reference: None,
            vcf: None,
            no_vbias: false,
            no_gc: false,
            tabix_path: "tabix".into(),
            bcftools_path: "bcftools".into(),
            rastair_path: "/tmp/custom-rastair".into(),
            threads: 1,
            wgbs: false,
            output_prefix: ".".into(),
            r_script_dir: None,
        };

        let flags = params.as_flags();

        assert!(flags.windows(2).any(|window| window == ["--rastair-path", "/tmp/custom-rastair"]));
    }
}
