//! `rastair mbias`: the QC report.
//!
//! This module has no analysis logic: it shells out to `scripts/mbias.R`,
//! which renders `scripts/QC_report.Rmd`. Cutoffs, plots and the per-contig
//! `{chrom}_cutoffs.txt` files are computed there, so fixes to them go in the
//! `.Rmd`. Things the R code relies on:
//!
//! - OT/OB in the M-bias table matches the Rust pair orientation: OT if
//!   `bitwAnd(flag, 96) == 96` (F1R2) or `bitwAnd(flag, 144) == 144` (R2F1),
//!   OB otherwise (`80/160` swaps them). Read 5'→3' positions are flipped for
//!   OT `Second` and OB `First` only.
//! - Cutoffs are per `(chr, read_pair, strand)` group, which needs at least
//!   `MIN_MBIAS_OBS` (3) covered positions. A genome-wide run skips a sparse
//!   contig with a `warning()`; an explicitly requested one is a `stop()`.
//! - argparser flags default to `FALSE`, not `NA`, so test them with `isTRUE`.
//!   A `--bed`-only render skips the V-bias/GC/off-target chunks, and this
//!   module adds `--no-vbias`/`--no-gc`/`--no-off-target` when there is no
//!   `--reference`, and `--no-off-target` when there is no `--bam`.
//! - The GC/CpG-bias chunks compare reads with all read-long `N`-free windows
//!   of the reference, by CpG count and GC share taken from the reference
//!   (conversion changes the read's Cs), among those with at least one CpG,
//!   as the per-read data has no others. The GC-adjusted CpG line weights
//!   each window by its GC bin's observed/expected ratio. Reads and windows
//!   are also binned by their presumed number of methylated CpGs, the sum of
//!   their CpGs' levels; levels come from a random half of the reads and only
//!   the other half is counted, so a read never sets its own bin. A read's
//!   own calls are not used for this: a reference window has none, and a
//!   per-molecule count against level sums would be skewed by methylation
//!   being concordant along a molecule. Per-read offsets land on the CpG's C
//!   (OT) or G (OB).
//! - Off-target conversion: one `call --all` over `--region` (padded by the
//!   mask's flank) or the whole input, streamed through a `pipe()` and
//!   tallied a contig at a time, since a genome-wide `--all` VCF holds every
//!   reference CpG; reads are fetched in 2 Mb tiles. No `bcftools`, no
//!   temp files. Mismatch counts per strand come from `AS_SB_OT`/`AS_SB_OB`,
//!   split per alt allele in R. Positions near a `CPG` or `CPGnovo` record
//!   are masked: `--all` writes every reference CpG, and `CPGnovo` marks both
//!   bases of a CpG that a called variant creates. Coverage comes from reads
//!   via `scanBam`: mates are sorted by name and an overlapping pair merged
//!   into one fragment, since `call` counts fragments. The converted strand
//!   is OT at a reference C and OB at a G. Skipped with `--wgbs`, where
//!   `call`'s TAPS model does not apply.
//! - M-bias positions are parsed in 250,000-row chunks (`position_counts()`):
//!   pasting a whole `(chr, read_pair, strand)` group's position strings into
//!   one `scan()` scales with the largest chromosome and ran out of memory on
//!   deep data. V-bias reuses the parsed table, with the deduplicated per-read
//!   row index as `alignment_id`; a full-table `.GRP` over read names cost
//!   seconds and a chromosome-sized workspace.
//! - A chromosome-only report needs `--region`: a whole-genome BED is otherwise
//!   read whole (`bgzip -cd`), not through `tabix`.
//!
//! `tests/mbias_report.rs` renders a synthetic per-read BED (behind
//! `external-tool-tests`, skipped without R and its packages, `tabix` or
//! `bgzip`). On macOS the `.Rmd`'s `zcat` only reads `.Z`, so put a `zcat`
//! shim running `gzip -dc` first on `PATH`. To render without `mbias.R`:
//! `Rscript -e "rmarkdown::render('scripts/QC_report.Rmd', params=list(input_bgz=...,
//! output_dir=..., region=NA, plot_vbias=FALSE, plot_gc=FALSE))"` (`region=NA`,
//! not NULL).

use crate::utils::cli;
use clap::{ArgGroup, value_parser};
use clio::ClioPath;
use color_eyre::{
    Result, Section as _,
    eyre::{Context, ContextCompat, ensure, eyre},
};
use seqair_types::RegionString;
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

    /// Reference FASTA file (required for V-bias, GC/CpG bias and off-target plots)
    #[arg(long, value_parser=value_parser!(ClioPath).exists().is_file(), value_hint=clap::ValueHint::FilePath)]
    #[arg(help_heading = cli::sections::INPUT)]
    pub reference: Option<ClioPath>,

    /// Do not generate V-bias plots (faster)
    #[arg(long = "no-vbias")]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub no_vbias: bool,

    /// Do not generate GC/CpG bias plots
    #[arg(long = "no-gc")]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub no_gc: bool,

    /// Do not estimate off-target conversion
    ///
    /// The estimate needs `--bam` and `--reference`, and is skipped without them.
    #[arg(long = "no-off-target")]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub no_off_target: bool,

    /// Path to tabix executable
    #[arg(long = "tabix-path", default_value = "tabix")]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub tabix_path: String,

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
        if self.wgbs {
            params.push("--wgbs".into());
        }
        // auto-disable plots that require inputs the user hasn't provided
        if self.no_vbias || self.reference.is_none() {
            params.push("--no-vbias".into());
        }
        if self.no_gc || self.reference.is_none() {
            params.push("--no-gc".into());
        }
        if self.no_off_target || self.reference.is_none() || self.bam_file.is_none() {
            params.push("--no-off-target".into());
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
    use clio::ClioPath;
    use std::path::Path;

    fn params() -> MBiasParams {
        MBiasParams {
            bed_file: None,
            bam_file: None,
            region: None,
            include_flag: 3,
            exclude_flag: 3852,
            read_length: None,
            reference: None,
            no_vbias: false,
            no_gc: false,
            no_off_target: false,
            tabix_path: "tabix".into(),
            rastair_path: "rastair".into(),
            threads: 1,
            wgbs: false,
            output_prefix: ".".into(),
            r_script_dir: None,
        }
    }

    fn test_data(name: &str) -> ClioPath {
        ClioPath::new(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name))
            .expect("test data path should be valid")
    }

    #[test]
    fn as_flags_includes_extended_script_options() {
        let params = MBiasParams { threads: 4, wgbs: true, ..params() };

        let flags = params.as_flags();

        assert!(flags.windows(2).any(|window| window == ["--tabix-path", "tabix"]));
        assert!(flags.windows(2).any(|window| window == ["--threads", "4"]));
        assert!(flags.iter().any(|flag| flag == "--wgbs"));
        assert!(flags.iter().any(|flag| flag == "--no-vbias"));
        assert!(flags.iter().any(|flag| flag == "--no-gc"));
        assert!(flags.iter().any(|flag| flag == "--no-off-target"));
        assert!(flags.windows(2).any(|window| window[0] == "--rastair-path"));
    }

    #[test]
    fn as_flags_disables_off_target_only_without_bam_or_reference() {
        let reference = Some(test_data("test.fasta.gz"));
        let bam_file = Some(test_data("test.bam"));
        let has_flag =
            |params: MBiasParams| params.as_flags().iter().any(|flag| flag == "--no-off-target");

        assert!(!has_flag(MBiasParams {
            reference: reference.clone(),
            bam_file: bam_file.clone(),
            ..params()
        }));
        assert!(has_flag(MBiasParams { reference: reference.clone(), ..params() }));
        assert!(has_flag(MBiasParams { bam_file: bam_file.clone(), ..params() }));
        assert!(has_flag(MBiasParams { reference, bam_file, no_off_target: true, ..params() }));
    }

    #[test]
    fn as_flags_preserves_explicit_rastair_path_override() {
        let params = MBiasParams { rastair_path: "/tmp/custom-rastair".into(), ..params() };

        let flags = params.as_flags();

        assert!(flags.windows(2).any(|window| window == ["--rastair-path", "/tmp/custom-rastair"]));
    }
}
