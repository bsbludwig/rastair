//! Message Pack format tooling
//!
//! This is just the internal format used by rastair2. We expose it only for
//! debugging.

use crate::vcf;
use clio::ClioPath;
use color_eyre::eyre::{Context as _, Result, eyre};
use serde::de::DeserializeOwned;
use smol_str::SmolStr;
use std::{
    borrow::Cow,
    io::{BufReader, BufWriter, ErrorKind, Read, Write},
    marker::PhantomData,
};
use tracing::{instrument, trace};

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MpkHeader {
    pub rastair2_version: String,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MpkVcfHeader {
    pub contigs: Vec<SmolStr>,
    pub samples: Vec<SmolStr>,
    pub metadata: Vec<String>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[allow(clippy::large_enum_variant)] // all but two entries in a file are `Record`s
pub enum MpkEntry<'r> {
    Header(MpkHeader),
    VcfHeader(MpkVcfHeader),
    Record(Cow<'r, vcf::Record>),
}

pub struct MessagePackWriter {
    pub path: ClioPath,
    writer: Box<dyn Write + Send>,
}

impl MessagePackWriter {
    /// Create a new `MessagePackWriter` with the specified output path.
    #[instrument(level = "debug")]
    pub fn new(path: &ClioPath) -> Result<Self> {
        let file = BufWriter::new(
            path.clone().create().wrap_err_with(|| format!("Failed to create output {path}"))?,
        );

        let writer = lz4::EncoderBuilder::new()
            .level(0)
            .build(file)
            .wrap_err("Failed to create LZ4 encoder")?;
        let mut me = Self { path: path.clone(), writer: Box::new(writer) };
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

pub struct MessagePackReader {
    pub path: ClioPath,
    reader: Box<dyn Read + Send>,
}

pub struct MpkFile {
    pub header: MpkHeader,
    pub vcf_header: Option<MpkVcfHeader>,
    pub entries: Box<dyn Iterator<Item = Result<MpkEntry<'static>>> + Send>,
}

impl MessagePackReader {
    /// Create a new `MessagePackReader` with the specified input path.
    #[instrument(level = "debug")]
    pub fn new(path: &ClioPath) -> Result<Self> {
        let file = path.clone().open().wrap_err_with(|| format!("Failed to open {path}"))?;
        let reader =
            lz4::Decoder::new(BufReader::new(file)).wrap_err("Failed to create LZ4 decoder")?;
        Ok(Self { path: path.clone(), reader: Box::new(reader) })
    }

    pub fn read(self) -> Result<MpkFile> {
        let mut entries = self.read_entry().peekable();
        let header = match entries.next() {
            Some(Ok(MpkEntry::Header(header))) => header,
            Some(Ok(_)) => {
                return Err(eyre!("Expected header entry but found something else"));
            }
            Some(Err(e)) => return Err(e.wrap_err("Failed to read header from Message Pack file")),
            None => return Err(eyre!("No entries found in Message Pack file")),
        };
        let vcf_header = if let Some(Ok(MpkEntry::VcfHeader(header))) = entries.peek() {
            Some(header.clone())
        } else {
            None
        };
        if vcf_header.is_some() {
            entries.next(); // Consume the VCF header entry
        }

        Ok(MpkFile { header, vcf_header, entries: Box::new(entries) })
    }

    /// Read the entries from the Message Pack file
    fn read_entry(self) -> impl Iterator<Item = Result<MpkEntry<'static>>> {
        // Use a streaming deserializer to read entries one by one
        StreamingDeserializer { reader: self.reader, item: PhantomData }
    }
}

struct StreamingDeserializer<R: Read, T: DeserializeOwned> {
    reader: R,
    item: PhantomData<T>,
}

impl<R: Read, T: DeserializeOwned> Iterator for StreamingDeserializer<R, T> {
    type Item = Result<T>;

    fn next(&mut self) -> Option<Self::Item> {
        match rmp_serde::decode::from_read(&mut self.reader) {
            Ok(entry) => Some(Ok(entry)),
            Err(rmp_serde::decode::Error::InvalidMarkerRead(e))
            | Err(rmp_serde::decode::Error::InvalidDataRead(e))
                if e.kind() == ErrorKind::UnexpectedEof =>
            {
                // If we hit EOF, we just return None to end the iteration
                trace!(%e, "reached end of Message Pack file");
                None
            }
            error => Some(error.wrap_err("Failed to decode structure from MessagePack")),
        }
    }
}
