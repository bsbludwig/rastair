#![allow(unused_imports, dead_code, reason = "test code")]

pub use color_eyre::eyre::{bail, ensure, eyre};
pub use color_eyre::{Result, eyre::Context as _};
pub use insta::{assert_debug_snapshot, assert_snapshot};
pub use insta_cmd::assert_cmd_snapshot;
pub use std::{collections::BTreeSet, path::Path, process::Command};
pub use tempfile::TempDir;

pub const CALL_TEST_BAM: [&str; 3] =
    ["call", "--fasta-file=tests/data/test.fasta.gz", "tests/data/test.bam"];
pub const CHR19_SMALL: &str = "--region=chr19:6105700-6105800";
pub const NO_ML: &str = "--no-ml"; // disable ML for faster tests

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
        settings.add_filter(r"\[TIME\] INFO rastair: Using experimental seqair backend\.\n", "");
        settings.add_filter(r"duration=[\w.]+", "[DURATION]");
        settings.add_filter(r": close time.*", " [CLOSE]");
        settings.add_filter(r#"file="/.*/*.vcf"#, "file=[PATH]");
        settings.add_filter(r#"file="/.*/*.vcf.gz"#, "file=[PATH]");
        settings.add_filter(r#"file="/.*/*.bcf"#, "file=[PATH]");
        settings.add_filter(r#"file="/.*/*.mpk.lz4"#, "file=[PATH]");
        settings.add_filter(r#"file="/.*/*.bed"#, "file=[PATH]");
        settings.add_filter(r#"/var/.*/*.bam"#, "[PATH]");
        settings.add_filter(r#"/tmp/.*/*.bam"#, "[PATH]");
        settings.add_filter(r#"/var/.*/*.bed"#, "[PATH]");
        settings.add_filter(r#"/tmp/.*/*.gz"#, "[PATH]");
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

/// Writes a plain FASTA holding only `contigs` (one short dummy sequence line each)
/// plus a matching `.fai` index into `dir`, returning the FASTA path. Used to build a
/// reference that deliberately omits a contig the test BAM carries.
pub fn write_fasta_with_contigs(dir: &Path, contigs: &[&str]) -> Result<std::path::PathBuf> {
    const SEQ: &str = "ACGTACGTAC";
    let fasta = dir.join("ref.fa");
    let mut body = String::new();
    let mut fai = String::new();
    for contig in contigs {
        let header = format!(">{contig}\n");
        // Offset of the sequence's first base is the current byte length of the file.
        let offset = body.len() + header.len();
        body.push_str(&header);
        body.push_str(SEQ);
        body.push('\n');
        // fai columns: NAME, LENGTH, OFFSET, LINEBASES, LINEWIDTH (bases + newline).
        fai.push_str(&format!(
            "{contig}\t{}\t{offset}\t{}\t{}\n",
            SEQ.len(),
            SEQ.len(),
            SEQ.len() + 1
        ));
    }
    std::fs::write(&fasta, body).wrap_err("write fasta")?;
    std::fs::write(dir.join("ref.fa.fai"), fai).wrap_err("write fai")?;
    Ok(fasta)
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

pub trait CommandStdioExt {
    fn silent(&mut self) -> &mut Self;
}

impl CommandStdioExt for std::process::Command {
    fn silent(&mut self) -> &mut Self {
        self.stdout(std::process::Stdio::null()).stderr(std::process::Stdio::null())
    }
}
