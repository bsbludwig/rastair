use crate::{
    io::mpk::format::{MpkEntry, MpkHeader, MpkVcfHeader},
    metrics::PileupMetrics,
};
use clio::{ClioPath, Output};
use color_eyre::eyre::{Context as _, Result};
use std::{borrow::Cow, io::BufWriter};
use tracing::instrument;

pub struct MessagePackWriter {
    pub path: ClioPath,
    writer: BufWriter<lz4::Encoder<Output>>,
}

impl MessagePackWriter {
    /// Create a new `MessagePackWriter` with the specified output path.
    #[instrument(level = "debug", skip_all)]
    pub fn new(path: &ClioPath) -> Result<Self> {
        let file =
            path.clone().create().wrap_err_with(|| format!("Failed to create output {path}"))?;

        let writer = lz4::EncoderBuilder::new()
            .level(0)
            .block_size(lz4::BlockSize::Max1MB)
            .build(file)
            .wrap_err("Failed to create LZ4 encoder")?;
        let one_mb = 1024 * 1024;
        let mut me = Self { path: path.clone(), writer: BufWriter::with_capacity(one_mb, writer) };
        me.write(&MpkEntry::Header(MpkHeader {
            rastair_version: env!("CARGO_PKG_VERSION").into(),
        }))?;
        Ok(me)
    }

    pub fn add_metadata(&mut self, data: MpkVcfHeader) -> Result<()> {
        self.write(&MpkEntry::VcfHeader(data.clone()))
            .wrap_err("Failed to write VCF header to Message Pack file")
    }

    /// Write a record to the Message Pack file.
    pub fn add(&mut self, record: &PileupMetrics) -> Result<()> {
        self.write(&MpkEntry::Record(Cow::Borrowed(record))).wrap_err("Failed to write record")
    }

    fn write(&mut self, entry: &MpkEntry) -> Result<()> {
        rmp_serde::encode::write(&mut self.writer, entry)
            .wrap_err("Failed to write entry to MessagePack file")
    }

    /// Flush and close the file. Without this, the output is truncated.
    pub fn finish(self) -> Result<()> {
        let encoder = self
            .writer
            .into_inner()
            .map_err(std::io::IntoInnerError::into_error)
            .wrap_err("Failed to flush MessagePack output")?;
        let (output, ended) = encoder.finish();
        ended.wrap_err("Failed to end the LZ4 stream")?;
        output.finish().wrap_err("Failed to close MessagePack output")
    }
}
