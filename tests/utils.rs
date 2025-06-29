#![allow(unused_imports, dead_code)]

use std::process::Command;

use color_eyre::eyre::bail;
pub use color_eyre::{Result, eyre::Context as _};
pub use insta::{assert_debug_snapshot, assert_snapshot};
pub use insta_cmd::assert_cmd_snapshot;
pub use tempfile::TempDir;

pub fn rastair() -> Command {
    let mut cmd = Command::new(insta_cmd::get_cargo_bin("rastair2"));
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
        let _bound = settings.bind_to_scope();
    }
}

pub trait ExitStatusResultExt {
    fn is_success(&self) -> Result<()>;
}

impl ExitStatusResultExt for std::process::ExitStatus {
    fn is_success(&self) -> Result<()> {
        if !self.success() {
            bail!("Command failed with status: {}", self)
        }
        Ok(())
    }
}
