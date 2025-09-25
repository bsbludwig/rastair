use std::io::Write as _;

use clio::ClioPath;
use color_eyre::{Result, eyre::Context};

use crate::{
    io::mpk::{MessagePackReader, MpkEntry},
    utils::cli,
};

#[derive(Debug, clap::Args)]
pub struct MpkViewParams {
    /// Message Pack file to view
    pub input: ClioPath,

    /// Message Pack file to view
    #[arg(short = 'o', long, default_value = "-")]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub output: ClioPath,
}

pub fn view(params: &MpkViewParams) -> Result<()> {
    let reader = MessagePackReader::new(&params.input)
        .and_then(|r| r.read())
        .wrap_err_with(|| format!("Failed to read input data from {}", params.input))?;
    let mut out = params
        .output
        .clone()
        .create()
        .wrap_err_with(|| format!("Failed to open output {}", params.output))?;

    serde_json::to_writer(&mut out, &MpkEntry::Header(reader.header))?;
    out.write_all(b"\n")?;

    if let Some(vcf_header) = reader.vcf_header {
        serde_json::to_writer(&mut out, &MpkEntry::VcfHeader(vcf_header))?;
        out.write_all(b"\n")?;
    }

    for entry in reader.entries {
        let entry = entry?;
        serde_json::to_writer(&mut out, &entry)?;
        out.write_all(b"\n")?;
    }

    Ok(())
}
