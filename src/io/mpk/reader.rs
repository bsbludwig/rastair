use crate::io::mpk::format::{MpkEntry, MpkHeader, MpkVcfHeader};
use clio::ClioPath;
use color_eyre::eyre::{Context as _, Result, eyre};
use serde::de::DeserializeOwned;
use std::{
    io::{BufReader, ErrorKind, Read},
    marker::PhantomData,
};
use tracing::{instrument, trace};

pub struct MessagePackReader {
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
        Ok(Self { reader: Box::new(reader) })
    }

    #[instrument(level = "debug", skip(self))]
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
