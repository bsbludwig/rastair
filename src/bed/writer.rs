use crate::bed::{BedFormat, BedRecord};
use clio::ClioPath;
use color_eyre::{Result, eyre::Context as _};
use std::{
    any::type_name,
    fmt,
    fs::File,
    io::{BufWriter, Write},
    marker::PhantomData,
    path::Path,
};
use tracing::instrument;

pub struct BedWriter<R: BedRecord> {
    pub path: ClioPath,
    pub format: BedFormat,
    writer: Writer,
    record_type: PhantomData<R>,
}

enum Writer {
    Bed(Box<dyn Write + Send + Sync>),
    BedGz(bgzip::BGZFWriter<Box<dyn Write + Send + Sync>>),
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
    #[instrument(level = "debug")]
    pub fn new(path: &ClioPath, format: BedFormat) -> Result<Self> {
        let writer = path.clone().create().wrap_err("Failed to create output")?;
        let writer: Box<dyn Write + Send + Sync> = Box::new(BufWriter::new(writer));
        let mut writer: Writer = match format {
            BedFormat::BedGz => {
                let writer = bgzip::BGZFWriter::with_compress_unit_size(
                    writer,
                    bgzip::Compression::fast(),
                    bgzip::write::DEFAULT_COMPRESS_UNIT_SIZE,
                    // Write index if the path is a local file
                    path.is_local(),
                )
                .wrap_err("Failed to create BGZF writer")?;
                Writer::BedGz(writer)
            }
            BedFormat::Bed => Writer::Bed(writer),
        };
        writeln!(&mut writer, "{}", R::HEADER).wrap_err("Failed to write header")?;
        Ok(Self { path: path.clone(), format, writer, record_type: PhantomData })
    }

    pub fn write_record(&mut self, record: &R) -> Result<()> {
        record.write(&mut self.writer)?;
        writeln!(self.writer)?;
        Ok(())
    }

    #[instrument(level = "debug")]
    pub fn close(mut self) -> Result<()> {
        self.writer.flush().wrap_err("Failed to flush writer")?;
        if let Writer::BedGz(bgzfwriter) = self.writer
            && let Some(index) = bgzfwriter.close().wrap_err("Failed to close BGZF writer")?
            && self.path.is_local()
        {
            write_index(self.path.path(), index)
                .wrap_err_with(|| format!("Failed to write index for `{}`", self.path.display()))?;
        }
        Ok(())
    }
}

impl<R: BedRecord> fmt::Debug for BedWriter<R> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("BedWriter")
            .field("type", &type_name::<R>())
            .field("path", &self.path)
            .field("format", &self.format)
            .finish()
    }
}

#[instrument(level = "debug")]
fn write_index(original_path: &Path, index: bgzip::index::BGZFIndex) -> Result<()> {
    let index_path = original_path.with_extension("gz.gzi");
    let mut index_file = File::create(&index_path)
        .wrap_err_with(|| format!("Failed to create index file `{}`", index_path.display()))?;
    index
        .write(&mut index_file)
        .wrap_err_with(|| format!("Failed to write index to `{}`", index_path.display()))?;
    Ok(())
}
