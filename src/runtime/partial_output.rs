//! Output files that only appear at their path once the run has succeeded.
//!
//! A run that fails closes its outputs properly, so they are valid files, just
//! incomplete. At the final path, nothing would tell them apart from complete
//! ones: not htslib (the BGZF EOF block is there), and not a workflow that
//! reruns samples whose output is missing. So outputs are written to
//! `<path>.partial` (and their indices next to that), and moved into place by
//! [`PartialOutput::commit`] only when the run succeeded. This also covers a
//! run that is killed (by a signal, the OOM killer, a power cut), which leaves
//! its `.partial` files cut off wherever it was.

use clio::ClioPath;
use color_eyre::eyre::{Context as _, Result};
use seqair::vcf::CoordinateIndex;
use std::path::{Path, PathBuf};
use tracing::{info, warn};

const PARTIAL_SUFFIX: &str = "partial";

/// Indices that writers create next to their output, as `<output>.<suffix>`.
const INDEX_SUFFIXES: [&str; 2] = [CoordinateIndex::SUFFIX, "tbi"];

/// An output file written under a temporary name until [`Self::commit`].
///
/// Dropping it without committing leaves the partial file where it is, and
/// logs where that is.
#[derive(Debug)]
#[must_use = "the output is only moved into place by `commit`"]
pub struct PartialOutput {
    /// What the output is, for logging, e.g. "VCF output"
    what: &'static str,
    target: ClioPath,
    /// `None` when the target is not a regular file (stdout, a pipe), which is
    /// written directly.
    partial: Option<Partial>,
}

#[derive(Debug)]
struct Partial {
    path: PathBuf,
    /// The target, with symlinks resolved: the output is written through a
    /// link, not over it.
    destination: PathBuf,
}

impl PartialOutput {
    pub fn new(what: &'static str, target: &ClioPath) -> Self {
        let partial = is_regular_file_target(target).then(|| {
            let destination =
                std::fs::canonicalize(target.path()).unwrap_or_else(|_| target.path().to_owned());
            let path = with_suffix(&destination, PARTIAL_SUFFIX);
            // An index left by an earlier run would be committed with this
            // output if this run fails to write one
            for suffix in INDEX_SUFFIXES {
                remove_if_exists(&with_suffix(&path, suffix));
            }
            Partial { path, destination }
        });
        Self { what, target: target.clone(), partial }
    }

    /// Where to write the output to.
    pub fn write_path(&self) -> Result<ClioPath> {
        match &self.partial {
            Some(Partial { path, .. }) => ClioPath::new(path)
                .wrap_err_with(|| format!("Invalid output path `{}`", path.display())),
            None => Ok(self.target.clone()),
        }
    }

    /// Move the output, and any index written next to it, into place.
    ///
    /// Indices are moved after the output, so they are never older than it
    /// (tabix warns about that). The target's old indices are removed first:
    /// if moving the new ones fails, a missing index is safe, a wrong one isn't.
    pub fn commit(mut self) -> Result<()> {
        if let Some(Partial { path, destination }) = self.partial.take() {
            for suffix in INDEX_SUFFIXES {
                let index = with_suffix(&destination, suffix);
                std::fs::remove_file(&index).or_else(ignore_not_found).wrap_err_with(|| {
                    format!("Failed to remove the old index `{}`", index.display())
                })?;
            }
            rename(&path, &destination)?;
            for suffix in INDEX_SUFFIXES {
                let partial_index = with_suffix(&path, suffix);
                if partial_index.exists() {
                    rename(&partial_index, &with_suffix(&destination, suffix))?;
                }
            }
        }
        info!(file = %self.target, "Wrote {}", self.what);
        Ok(())
    }
}

impl Drop for PartialOutput {
    fn drop(&mut self) {
        if let Some(Partial { path, .. }) = &self.partial
            && path.exists()
        {
            warn!(file = %path.display(), "{} left behind under its partial name", self.what);
        }
    }
}

/// Stdout, pipes and device files like `/dev/null` can't be renamed into.
fn is_regular_file_target(target: &ClioPath) -> bool {
    if target.is_std() || !target.is_local() {
        return false;
    }
    match std::fs::metadata(target.path()) {
        Ok(metadata) => metadata.is_file(),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => true,
        Err(error) => {
            warn!(?error, path = %target, "Can't tell what kind of file the output is, writing to it directly");
            false
        }
    }
}

fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut path = path.as_os_str().to_owned();
    path.push(".");
    path.push(suffix);
    PathBuf::from(path)
}

fn ignore_not_found(error: std::io::Error) -> std::io::Result<()> {
    match error.kind() {
        std::io::ErrorKind::NotFound => Ok(()),
        _ => Err(error),
    }
}

fn remove_if_exists(path: &Path) {
    if let Err(error) = std::fs::remove_file(path).or_else(ignore_not_found) {
        warn!(?error, file = %path.display(), "Failed to remove a stale file");
    }
}

fn rename(from: &Path, to: &Path) -> Result<()> {
    std::fs::rename(from, to)
        .wrap_err_with(|| format!("Failed to move `{}` to `{}`", from.display(), to.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn clio_path(path: &Path) -> ClioPath {
        ClioPath::new(path).expect("valid path")
    }

    #[test]
    fn commit_moves_output_and_index_into_place() -> Result<()> {
        let dir = tempfile::TempDir::new()?;
        let target = dir.path().join("out.vcf.gz");
        let output = PartialOutput::new("VCF output", &clio_path(&target));

        let write_path = output.write_path()?;
        assert_eq!(write_path.path(), dir.path().join("out.vcf.gz.partial"));
        std::fs::write(write_path.path(), "records")?;
        std::fs::write(dir.path().join("out.vcf.gz.partial.csi"), "index")?;

        output.commit()?;
        assert_eq!(std::fs::read_to_string(&target)?, "records");
        assert_eq!(std::fs::read_to_string(dir.path().join("out.vcf.gz.csi"))?, "index");
        assert!(!dir.path().join("out.vcf.gz.partial").exists());
        assert!(!dir.path().join("out.vcf.gz.partial.csi").exists());
        assert!(!dir.path().join("out.vcf.gz.tbi").exists());
        Ok(())
    }

    #[test]
    fn stale_indices_are_never_committed() -> Result<()> {
        let dir = tempfile::TempDir::new()?;
        let target = dir.path().join("out.bed.gz");
        std::fs::write(dir.path().join("out.bed.gz.tbi"), "old index")?;
        std::fs::write(dir.path().join("out.bed.gz.partial.tbi"), "stopped run's index")?;

        let output = PartialOutput::new("BED output", &clio_path(&target));
        std::fs::write(output.write_path()?.path(), "records")?;
        output.commit()?;
        assert_eq!(std::fs::read_to_string(&target)?, "records");
        assert!(!dir.path().join("out.bed.gz.tbi").exists());
        Ok(())
    }

    #[cfg(unix)]
    #[test]
    fn writes_through_a_symlink() -> Result<()> {
        let dir = tempfile::TempDir::new()?;
        let (file, link) = (dir.path().join("file.bed"), dir.path().join("link.bed"));
        std::fs::write(&file, "old")?;
        std::os::unix::fs::symlink(&file, &link)?;

        let output = PartialOutput::new("BED output", &clio_path(&link));
        std::fs::write(output.write_path()?.path(), "records")?;
        output.commit()?;
        assert!(std::fs::symlink_metadata(&link)?.file_type().is_symlink());
        assert_eq!(std::fs::read_to_string(&file)?, "records");
        Ok(())
    }

    #[test]
    fn uncommitted_output_stays_partial() -> Result<()> {
        let dir = tempfile::TempDir::new()?;
        let target = dir.path().join("out.bed");
        let output = PartialOutput::new("BED output", &clio_path(&target));
        std::fs::write(output.write_path()?.path(), "records")?;
        drop(output);
        assert!(!target.exists());
        assert!(dir.path().join("out.bed.partial").exists());
        Ok(())
    }

    #[test]
    fn stdout_and_devices_are_written_directly() -> Result<()> {
        let stdout = PartialOutput::new("VCF output", &ClioPath::new("-")?);
        assert!(stdout.write_path()?.is_std());
        stdout.commit()?;

        #[cfg(unix)]
        {
            let dev_null = clio_path(Path::new("/dev/null"));
            let output = PartialOutput::new("VCF output", &dev_null);
            assert_eq!(output.write_path()?.path(), Path::new("/dev/null"));
            output.commit()?;
        }
        Ok(())
    }
}
