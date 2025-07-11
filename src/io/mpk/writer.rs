use crate::{
    io::mpk::format::{MpkEntry, MpkHeader, MpkVcfHeader},
    vcf,
};
use clio::ClioPath;
use color_eyre::eyre::{Context as _, Result};
use std::{
    borrow::Cow,
    io::{BufWriter, Write},
};
use tracing::instrument;

pub struct MessagePackWriter {
    pub path: ClioPath,
    writer: Box<dyn Write + Send>,
}

impl MessagePackWriter {
    /// Create a new `MessagePackWriter` with the specified output path.
    #[instrument(level = "debug")]
    pub fn new(path: &ClioPath) -> Result<Self> {
        let file =
            path.clone().create().wrap_err_with(|| format!("Failed to create output {path}"))?;

        let writer = lz4::EncoderBuilder::new()
            .level(0)
            .build(file)
            .wrap_err("Failed to create LZ4 encoder")?;
        let mut me = Self { path: path.clone(), writer: Box::new(BufWriter::new(writer)) };
        me.write(&MpkEntry::Header(MpkHeader {
            rastair2_version: env!("CARGO_PKG_VERSION").into(),
        }))?;
        Ok(me)
    }

    pub fn add_metadata(&mut self, data: MpkVcfHeader) -> Result<()> {
        self.write(&MpkEntry::VcfHeader(data.clone()))
            .wrap_err("Failed to write VCF header to Message Pack file")
    }

    /// Write a record to the Message Pack file.
    pub fn add(&mut self, record: &vcf::Record) -> Result<()> {
        self.write(&MpkEntry::Record(Cow::Borrowed(record))).wrap_err("Failed to write record")
    }

    fn write(&mut self, entry: &MpkEntry) -> Result<()> {
        rmp_serde::encode::write(&mut self.writer, entry)
            .wrap_err("Failed to write entry to MessagePack file")
    }
}
