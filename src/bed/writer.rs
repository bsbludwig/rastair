use crate::bed::{BedFormat, BedRecord};
use clio::ClioPath;
use color_eyre::{Result, eyre::Context as _};
use rust_htslib::tbx;
use std::{
    any::type_name,
    fmt,
    io::{BufWriter, Write},
    marker::PhantomData,
};
use tracing::{debug, info, instrument, warn};

pub struct BedWriter<R: BedRecord> {
    pub path: ClioPath,
    pub format: BedFormat,
    writer: Writer,
    record_type: PhantomData<R>,
}

type BoxedWriter = Box<dyn Write + Send + Sync>;

enum Writer {
    /// BED as text files
    Bed(BoxedWriter),
    /// BGZF-compressed BED files
    BedGz(bgzf::Writer<BoxedWriter>),
}

impl Write for Writer {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Writer::Bed(writer) => writer.write(buf),
            Writer::BedGz(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Writer::Bed(writer) => writer.flush(),
            Writer::BedGz(writer) => writer.flush(),
        }
    }
}

impl<R: BedRecord> BedWriter<R> {
    #[instrument(level = "info", name = "BedWriter", skip(path), fields(path = %path))]
    pub fn new(path: &ClioPath, format: BedFormat) -> Result<Self> {
        let writer = path.clone().create().wrap_err("Failed to create output")?;
        let writer: Box<dyn Write + Send + Sync> = Box::new(BufWriter::new(writer));
        let mut writer: Writer = match format {
            BedFormat::BedGz => {
                let compression_level = bgzf::CompressionLevel::try_from(6)
                    .wrap_err("Failed to create compression level")?;
                Writer::BedGz(bgzf::Writer::new(writer, compression_level))
            }
            BedFormat::Bed => Writer::Bed(writer),
        };
        writeln!(&mut writer, "{}", R::HEADER).wrap_err("Failed to write header")?;
        debug!("Writing reads to BED");
        Ok(Self { path: path.clone(), format, writer, record_type: PhantomData })
    }

    pub fn write_record(&mut self, record: &R) -> Result<()> {
        record.write(&mut self.writer).wrap_err("Failed to write record")?;
        Ok(())
    }

    /// Close the writer and create a tabix index for BGZF files
    #[instrument(level = "debug", skip(self))]
    pub fn close(self) -> Result<()> {
        match self.writer {
            Writer::Bed(mut writer) => writer.flush().wrap_err("Failed to flush BED writer")?,
            Writer::BedGz(writer) => {
                writer.finish().wrap_err("Failed to finish BGZF writer")?;
                if self.path.is_file() {
                    let path = self.path.path();
                    match tbx::build_index(path, tbx::TabixFormat::Bed, 0) {
                        Ok(()) => info!(path = %path.display(), "Created tabix index"),
                        Err(error) => warn!(
                            %error,
                            path = %path.display(),
                            "Failed to create tabix index. \
                             You can create it manually with: tabix {}",
                            path.display()
                        ),
                    }
                }
            }
        }
        Ok(())
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl<R: BedRecord> fmt::Debug for BedWriter<R> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("BedWriter")
            .field("type", &type_name::<R>())
            .field("path", &self.path)
            .field("format", &self.format)
            .finish()
    }
}
