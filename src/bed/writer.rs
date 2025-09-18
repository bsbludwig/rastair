use crate::{
    bed::{BedFormat, BedRecord},
    utils::logging::ThisIsABug,
};
use clio::ClioPath;
use color_eyre::{
    Result,
    eyre::{Context as _, ContextCompat},
};
use noodles::{
    bgzf::{self, io::writer::CompressionLevel},
    core::Position,
    csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk},
    tabix,
};
use std::{
    any::type_name,
    fmt,
    fs::File,
    io::{BufWriter, Write},
    marker::PhantomData,
    path::PathBuf,
};
use tracing::{debug, info, instrument};

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
    /// BGZF-compressed BED files, will write both tbi index and gzi index files
    /// if the path is a local file
    BedGz { bed: bgzf::io::Writer<BoxedWriter>, tabix: Option<TabixWriter> },
}

struct TabixWriter {
    indexer: tabix::index::Indexer,
    output_path: PathBuf,
    output: tabix::io::Writer<BoxedWriter>,
}

impl Write for Writer {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Writer::Bed(writer) => writer.write(buf),
            Writer::BedGz { bed: writer, .. } => writer.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Writer::Bed(writer) => writer.flush(),
            Writer::BedGz { bed: writer, .. } => writer.flush(),
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
                let writer = bgzf::io::writer::Builder::default()
                    .set_compression_level(CompressionLevel::default())
                    .build_from_writer(writer);
                let tabix = if path.is_file() {
                    let mut indexer = tabix::index::Indexer::default();
                    indexer.set_header(csi::binning_index::index::header::Builder::bed().build());

                    let output_path = path.with_file_name({
                        let mut p = path
                            .path()
                            .file_name()
                            .wrap_err_with(|| {
                                format!(
                                    "Failed to get file name from path `{}`",
                                    path.path().display()
                                )
                            })?
                            .to_os_string();
                        p.push(".tbi");
                        p
                    });
                    let output_writer: BoxedWriter = Box::new(BufWriter::new(
                        File::create(&output_path).wrap_err_with(|| {
                            format!(
                                "Failed to create tabix index file for `{}`",
                                path.path().display()
                            )
                        })?,
                    ));
                    let output = tabix::io::Writer::new(output_writer);
                    debug!(path=?path.path(), "Writing tabix index");
                    Some(TabixWriter { indexer, output_path, output })
                } else {
                    None
                };
                Writer::BedGz { bed: writer, tabix }
            }
            BedFormat::Bed => Writer::Bed(writer),
        };
        writeln!(&mut writer, "{}", R::HEADER).wrap_err("Failed to write header")?;
        debug!("Writing reads to BED");
        Ok(Self { path: path.clone(), format, writer, record_type: PhantomData })
    }

    pub fn write_record(&mut self, record: &R) -> Result<()> {
        if let Writer::BedGz { bed, tabix } = &mut self.writer
            && let Some(tabix) = tabix.as_mut()
        {
            let start_position = bed.virtual_position();

            record.write(bed).wrap_err("Failed to write record")?;

            let end_position = bed.virtual_position();
            let chunk = Chunk::new(start_position, end_position);

            tabix.indexer.add_record(
                record.chr(),
                Position::new(record.start())
                    .wrap_err("Failed to convert start to position")
                    .this_is_a_bug()?,
                Position::new(record.end())
                    .wrap_err("Failed to convert end to position")
                    .this_is_a_bug()?,
                chunk,
            )?;
        } else {
            record.write(&mut self.writer).wrap_err("Failed to write record")?;
        }
        Ok(())
    }

    /// Close the writer and write the index if applicable
    #[instrument(level = "debug", skip(self))]
    pub fn close(self) -> Result<()> {
        match self.writer {
            Writer::Bed(mut writer) => writer.flush().wrap_err("Failed to flush BED writer")?,
            Writer::BedGz { bed, tabix } => {
                bed.finish().wrap_err("Failed to finish BGZF writer")?;
                if let Some(mut tabix) = tabix {
                    debug!("Writing tabix index");
                    let index = tabix.indexer.build();
                    tabix.output.write_index(&index).wrap_err_with(|| {
                        format!("Failed to write tabix index to `{}`", tabix.output_path.display())
                    })?;
                    info!(path=?tabix.output_path, "Wrote tabix index");
                }
            }
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
