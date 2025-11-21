#![allow(unused_imports, dead_code, reason = "test code")]

use color_eyre::eyre::ensure;
pub use color_eyre::eyre::{bail, eyre};
pub use color_eyre::{Result, eyre::Context as _};
pub use insta::{assert_debug_snapshot, assert_snapshot};
pub use insta_cmd::assert_cmd_snapshot;
use std::path::Path;
pub use std::{collections::BTreeSet, process::Command};
pub use tempfile::TempDir;

pub fn rastair() -> Command {
    let mut cmd = Command::new(insta_cmd::get_cargo_bin("rastair"));
    cmd.env("NO_COLOR", "1");
    cmd.env("RASTAIR_THREADS", "2");
    cmd
}

#[macro_export]
macro_rules! apply_common_filters {
    {} => {
        let mut settings = insta::Settings::clone_current();
        settings.add_filter(r"\w{4}-[0-9T\-:.]+Z\s", "[TIME]");
        settings.add_filter(r"duration=[\w.]+", "[DURATION]");
        settings.add_filter(r": close time.*", " [CLOSE]");
        settings.add_filter(r#"file="/.*/test.vcf"#, "file=[PATH]");
        settings.add_filter(r#"file="/.*/test.vcf.gz"#, "file=[PATH]");
        settings.add_filter(r#"file="/.*/test.bcf"#, "file=[PATH]");
        settings.add_filter(r#"file="/.*/test.mpk.lz4"#, "file=[PATH]");
        settings.add_filter(r#"file="/.*/test.bed"#, "file=[PATH]");
        settings.add_filter(r#"/var/.*/test.bam"#, "[PATH]");
        settings.add_filter(r#"/tmp/.*/test.bam"#, "[PATH]");
        settings.add_filter(r#"/var/.*/calls.bed.gz"#, "[PATH]");
        settings.add_filter(r#"/tmp/.*/calls.bed.gz"#, "[PATH]");
        let _bound = settings.bind_to_scope();
    }
}

pub trait ExitStatusResultExt {
    fn succeeds(&mut self) -> Result<()>;
}

impl ExitStatusResultExt for std::process::Command {
    #[track_caller]
    fn succeeds(&mut self) -> Result<()> {
        let mut status = self.status().wrap_err("Failed to run command")?;
        status.succeeds()
    }
}

impl ExitStatusResultExt for std::process::ExitStatus {
    #[track_caller]
    fn succeeds(&mut self) -> Result<()> {
        if !self.success() {
            bail!("Command failed with status: {}", self)
        }
        Ok(())
    }
}

impl ExitStatusResultExt for std::process::Output {
    #[track_caller]
    fn succeeds(&mut self) -> Result<()> {
        self.status.succeeds()?;
        Ok(())
    }
}

pub trait StringOutputExt {
    fn stdout(&self) -> String;
    fn stderr(&self) -> String;
}

impl StringOutputExt for std::process::Output {
    fn stdout(&self) -> String {
        String::from_utf8_lossy(&self.stdout).to_string()
    }

    fn stderr(&self) -> String {
        String::from_utf8_lossy(&self.stderr).to_string()
    }
}

pub trait StrInteratorToStringExt {
    fn collect_string(&mut self) -> String;
}

impl<'a, I: Iterator<Item = &'a str>> StrInteratorToStringExt for I {
    fn collect_string(&mut self) -> String {
        self.map(|s| s.to_string() + "\n").collect()
    }
}

pub fn vcf_content_lines(vcf_text: &str) -> impl Iterator<Item = &str> {
    vcf_text.lines().filter(|line| !line.starts_with("#"))
}

/// Read a BCF file and ensure it has contigs and at least one record.
pub fn read_bcf(path: &Path) -> Result<rust_htslib::bcf::Reader> {
    use rust_htslib::bcf::{Read, Reader};
    let mut bcf = Reader::from_path(path).wrap_err("open bcf file")?;

    ensure!(bcf.header().contig_count() > 0, "bcf file has no contigs");

    bcf.records()
        .next()
        .ok_or_else(|| eyre!("no records in bcf file"))?
        .wrap_err("failed to read first record")?;

    Ok(bcf)
}
