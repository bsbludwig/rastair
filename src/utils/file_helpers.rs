use bio::io::fasta::IndexedReader;
use color_eyre::{
    Section,
    eyre::{Context, ContextCompat, Result, eyre},
};
use noodles_bgzf as bgzf;
use std::{
    fmt,
    fs::File,
    io::{Read, Seek},
    ops::{Deref, DerefMut},
    path::{Path, PathBuf},
};
use tracing::{debug, instrument};

pub struct FastaReader {
    fasta_file: PathBuf,
    index_file: PathBuf,
    reader: IndexedReader<Box<dyn ReadAndSeek>>,
}

impl fmt::Debug for FastaReader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FastaReader")
            .field("fasta", &self.fasta_file)
            .field("index", &self.index_file)
            .finish()
    }
}

impl Deref for FastaReader {
    type Target = IndexedReader<Box<dyn ReadAndSeek>>;

    fn deref(&self) -> &Self::Target {
        &self.reader
    }
}

impl DerefMut for FastaReader {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.reader
    }
}

/// Combines `Read` and `Seek` traits.
pub trait ReadAndSeek: Read + Seek {}

/// This trait is implemented for any type that implements both `Read` and `Seek`.
impl<R: Read + Seek> ReadAndSeek for R {}

/// Open a FASTA file with a FAI index.
///
/// Assumes that next to the FASTA file, there is a corresponding `.fai` index file.
#[instrument(level = "debug")]
pub fn open_fasta(fasta_path: &Path) -> Result<FastaReader> {
    let fasta_file = open_maybe_bgzip(fasta_path)
        .wrap_err_with(|| format!("Failed to open FASTA file {fasta_path:?}"))?;
    let possible_index_files = [
        fasta_path.with_file_name({
            let mut name = fasta_path
                .file_name()
                .wrap_err("FASTA file path does not have a name")
                .note("This happened when looking for the index file")
                .note("Rastair already opened the FASTA file successfully")?
                .to_os_string();
            name.push(".fai");
            name
        }),
        fasta_path.with_extension("fai"),
        fasta_path.with_extension("fa.fai"),
        fasta_path.with_extension("gz.fai"),
    ];
    let index_path = possible_index_files.iter().find(|p| p.exists()).ok_or_else(|| {
        eyre!("No index file found for FASTA input {fasta_path:?}")
            .with_note(|| format!("Expected index file to be one of: {possible_index_files:?}"))
            .with_suggestion(|| format!("Create FASTA index with `samtools faidx {fasta_path:?}`"))
    })?;

    let fasta_index = bio::io::fasta::Index::from_file(&index_path)
        .map_err(|err| eyre!(Box::new(err)))
        .wrap_err_with(|| format!("Failed to read FASTA index file {index_path:?}"))
        .with_suggestion(|| {
            format!("You can recreate the FASTA index using `samtools faidx {fasta_path:?}`.")
        })?;
    Ok(FastaReader {
        reader: IndexedReader::with_index(fasta_file.into_reader(), fasta_index),
        fasta_file: fasta_path.to_path_buf(),
        index_file: index_path.to_path_buf(),
    })
}

#[cfg(test)]
mod fasta_tests {
    use super::*;
    use std::{io::Write, path::PathBuf};
    use tempfile::TempDir;

    #[test]
    fn test_open_fasta_success() -> Result<()> {
        let dir = TempDir::new()?;
        let fasta_path = dir.path().join("test.fa");
        let index_path = dir.path().join("test.fai");

        // Create a simple FASTA file
        let mut fasta_file = File::create(&fasta_path)?;
        fasta_file.write_all(b">seq1\nACGT\n>seq2\nGTCA\n")?;

        // Create a simple FAI index
        let mut index_file = File::create(&index_path)?;
        // Format: name, length, offset, line_bases, line_width
        index_file.write_all(b"seq1\t4\t6\t4\t5\nseq2\t4\t16\t4\t5\n")?;

        open_fasta(&fasta_path)?;

        Ok(())
    }

    #[test]
    fn test_open_fasta_missing_file() {
        let non_existent_path = PathBuf::from("/path/that/does/not/exist.fa");
        let result = open_fasta(&non_existent_path);
        assert!(result.is_err());
    }

    #[test]
    fn test_open_fasta_missing_index() -> Result<()> {
        let dir = TempDir::new()?;
        let fasta_path = dir.path().join("test.fa");

        // Create FASTA file without index
        let mut fasta_file = File::create(&fasta_path)?;
        fasta_file.write_all(b">seq1\nACGT\n")?;

        let result = open_fasta(&fasta_path);
        assert!(result.is_err());

        Ok(())
    }

    #[test]
    fn test_open_fasta_with_different_index_names() -> Result<()> {
        let dir = TempDir::new()?;
        let fasta_path = dir.path().join("test.fasta");
        let index_path = dir.path().join("test.fasta.fai");

        // Create a simple FASTA file
        let mut fasta_file = File::create(&fasta_path)?;
        fasta_file.write_all(b">seq1\nACGT\n>seq2\nGTCA\n")?;

        // Create a simple FAI index with .fai extension
        let mut index_file = File::create(&index_path)?;
        index_file.write_all(b"seq1\t4\t6\t4\t5\nseq2\t4\t16\t4\t5\n")?;

        open_fasta(&fasta_path)?;

        Ok(())
    }

    #[test]
    fn test_open_fasta_with_bgzip() -> Result<()> {
        let dir = TempDir::new()?;
        let fasta_path = dir.path().join("test.fa.gz");
        let fasta_index_path = dir.path().join("test.fa.fai");
        let bgzip_index_path = dir.path().join("test.fa.gz.gzi");

        // Create a bgzipped FASTA file with gzi index
        let mut gz = bgzf::io::Writer::new(File::create(&fasta_path)?);
        gz.write_all(b">seq1\nACGT\n>seq2\nGTCA\n")?;
        gz.finish()?;
        bgzf::gzi::fs::write(&bgzip_index_path, &bgzf::gzi::Index::default())?;

        // Create FAI index
        let mut index_file = File::create(&fasta_index_path)?;
        index_file.write_all(b"seq1\t4\t6\t4\t5\nseq2\t4\t16\t4\t5\n")?;

        open_fasta(&fasta_path)?;

        Ok(())
    }

    #[test]
    fn test_open_fasta_corrupt_index() -> Result<()> {
        let dir = TempDir::new()?;
        let fasta_path = dir.path().join("test.fa");
        let index_path = dir.path().join("test.fa.fai");

        // Create FASTA file
        let mut fasta_file = File::create(&fasta_path)?;
        fasta_file.write_all(b">seq1\nACGT\n")?;

        // Create corrupt FAI index
        let mut index_file = File::create(&index_path)?;
        index_file.write_all(b"corrupted_content")?;

        let result = open_fasta(&fasta_path);
        assert!(result.is_err());

        Ok(())
    }
}

enum BgzipReader {
    Gz { gz_file: PathBuf, index_file: PathBuf, reader: Box<dyn ReadAndSeek> },
    Uncompressed { file: PathBuf, reader: Box<dyn ReadAndSeek> },
}

impl BgzipReader {
    pub fn into_reader(self) -> Box<dyn ReadAndSeek> {
        match self {
            BgzipReader::Gz { reader, .. } => reader,
            BgzipReader::Uncompressed { reader, .. } => reader,
        }
    }
}

impl fmt::Debug for BgzipReader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BgzipReader::Gz { gz_file, index_file, .. } => f
                .debug_struct("BgzipReader")
                .field("gz_file", gz_file)
                .field("index_file", index_file)
                .finish(),
            BgzipReader::Uncompressed { file, .. } => {
                f.debug_struct("BgzipReader").field("file", file).finish()
            }
        }
    }
}

impl Deref for BgzipReader {
    type Target = Box<dyn ReadAndSeek>;

    fn deref(&self) -> &Self::Target {
        match self {
            BgzipReader::Gz { reader, .. } => reader,
            BgzipReader::Uncompressed { reader, .. } => reader,
        }
    }
}

impl DerefMut for BgzipReader {
    fn deref_mut(&mut self) -> &mut Self::Target {
        match self {
            BgzipReader::Gz { reader, .. } => reader,
            BgzipReader::Uncompressed { reader, .. } => reader,
        }
    }
}

/// Open the file at the path. If the file has a `.gz` extension, it is assumed
/// to be a bgzip-compressed file and that there is a corresponding `.gz.gzi`
/// index file.
#[instrument(level = "debug", skip_all)]
fn open_maybe_bgzip<P: AsRef<Path> + std::fmt::Debug>(path: P) -> Result<BgzipReader> {
    let path = path.as_ref();

    if path.extension().unwrap_or_default() == "gz" {
        let possible_index_files = [path.with_extension("gzi"), path.with_extension("gz.gzi")];
        let index_path = possible_index_files.iter().find(|p| p.exists()).ok_or_else(|| {
            eyre!("No index file found for bgzip input {path:?}")
                .with_note(|| format!("Expected index file to be one of: {possible_index_files:?}"))
                .with_suggestion(|| format!("Create bgzip index with `bgzip -r {path:?}`"))
        })?;
        let index = bgzf::gzi::io::Reader::new(
            open(index_path)
                .wrap_err_with(|| format!("Failed to open bgzip index file {index_path:?}"))?,
        )
        .read_index()
        .wrap_err_with(|| format!("Failed to parse bgzip index file {index_path:?}"))?;
        let bgzf_file =
            open(path).wrap_err_with(|| format!("Failed to open bgzip file {path:?}"))?;
        let in_file = bgzf::io::IndexedReader::new(bgzf_file, index);
        debug!(?index_path, "Opened bgzip file");
        Ok(BgzipReader::Gz {
            gz_file: path.to_path_buf(),
            index_file: index_path.to_path_buf(),
            reader: Box::new(in_file),
        })
    } else {
        let in_file = open(path)?;
        debug!("Opened uncompressed file");
        Ok(BgzipReader::Uncompressed { file: path.to_path_buf(), reader: in_file })
    }
}

fn open(path: &Path) -> Result<Box<dyn ReadAndSeek>> {
    let file = File::open(path).wrap_err_with(|| format!("Failed to open file {path:?}"))?;
    let file = std::io::BufReader::new(file);
    Ok(Box::new(file))
}

#[cfg(test)]
mod open_maybe_bgzip {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    fn write_dummy_file(path: &Path) -> Result<()> {
        let mut file = File::create(path)?;
        file.write_all(b"compressed content")?;
        Ok(())
    }

    #[test]
    fn test_open_uncompressed() -> Result<()> {
        let dir = TempDir::new()?;
        let file_path = dir.path().join("test_file.txt");

        write_dummy_file(&file_path)?;

        let _ = open_maybe_bgzip(&file_path)?;
        Ok(())
    }

    #[test]
    fn test_open_bgzip_with_index() -> Result<()> {
        let dir = TempDir::new()?;

        let file_path = dir.path().join("test_file.gz");
        let index_path = dir.path().join("test_file.gz.gzi");

        // create a dummy bgzip file with an index
        let mut gz = bgzf::io::Writer::new(File::create(&file_path)?);
        gz.write_all(b"compressed content :)")?;
        gz.finish()?;
        bgzf::gzi::fs::write(&index_path, &bgzf::gzi::Index::default())?;

        let _ = open_maybe_bgzip(&file_path)?;

        Ok(())
    }

    #[test]
    fn test_open_bgzip_fails_without_index() -> Result<()> {
        let dir = TempDir::new()?;

        let file_path = dir.path().join("test_file.gz");
        write_dummy_file(&file_path)?;

        // no index!

        let result = open_maybe_bgzip(&file_path);
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn test_open_file_error() {
        let path = PathBuf::from("/nonexistent/file.txt");
        let result = open_maybe_bgzip(&path);
        assert!(result.is_err());
    }

    #[test]
    fn test_no_index_file_error() -> Result<()> {
        let dir = TempDir::new()?;
        let file_path = dir.path().join("test.gz");

        // Create a file with .gz extension but no index
        let mut file = File::create(&file_path)?;
        file.write_all(b"some content")?;

        let result = open_maybe_bgzip(&file_path);
        assert!(result.is_err());

        Ok(())
    }

    #[test]
    fn test_read_index_error() -> Result<()> {
        let dir = TempDir::new()?;
        let file_path = dir.path().join("test.gz");
        let index_path = dir.path().join("test.gz.gzi");

        // Create the main file
        let mut file = File::create(&file_path)?;
        file.write_all(b"some content")?;

        // Create a corrupted/invalid index file
        let mut index_file = File::create(&index_path)?;
        index_file.write_all(b"this is not a valid bgzip index")?;

        let result = open_maybe_bgzip(&file_path);
        assert!(result.is_err());

        Ok(())
    }

    #[test]
    fn test_read_bgzip_error() -> Result<()> {
        let dir = TempDir::new()?;
        let file_path = dir.path().join("test.gz");
        let index_path = dir.path().join("test.gz.gzi");

        // Create a file that's not a valid BGZIP file
        let mut file = File::create(&file_path)?;
        file.write_all(b"this is not a valid bgzip file")?;

        // Create an empty GZI index so we get past the index check
        bgzf::gzi::fs::write(&index_path, &bgzf::gzi::Index::default())?;

        // Construction succeeds but reading fails on invalid BGZF data
        let mut reader = open_maybe_bgzip(&file_path)?.into_reader();
        let mut buf = [0u8; 1];
        assert!(reader.read(&mut buf).is_err());

        Ok(())
    }

    #[test]
    fn test_indexed_bgzip_reader_error() -> Result<()> {
        let dir = TempDir::new()?;
        let file_path = dir.path().join("test.gz");
        let index_path = dir.path().join("test.gz.gzi");

        // Create a valid BGZIP file
        let mut gz = bgzf::io::Writer::new(File::create(&file_path)?);
        gz.write_all(b"content")?;
        gz.finish()?;

        // Create an index that doesn't match the file
        let mut index = File::create(&index_path)?;
        write!(&mut index, "no thanks")?;

        let result = open_maybe_bgzip(&file_path);
        assert!(result.is_err());

        Ok(())
    }
}
