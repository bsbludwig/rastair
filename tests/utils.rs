#![allow(unused_imports, dead_code)]

pub use color_eyre::eyre::{bail, eyre};
pub use color_eyre::{Result, eyre::Context as _};
pub use insta::{assert_debug_snapshot, assert_snapshot};
pub use insta_cmd::assert_cmd_snapshot;
pub use std::{collections::BTreeSet, process::Command};
pub use tempfile::TempDir;

pub fn rastair() -> Command {
    let mut cmd = Command::new(insta_cmd::get_cargo_bin("rastair"));
    cmd.env("NO_COLOR", "1");
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
