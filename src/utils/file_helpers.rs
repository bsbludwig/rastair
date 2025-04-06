use bgzip::{BGZFReader, index::BGZFIndex, read::IndexedBGZFReader};
use bio::io::fasta::IndexedReader;
use color_eyre::eyre::{Result, bail, eyre};
use std::{
    fs::File,
    io::{Read, Seek},
    path::{Path, PathBuf},
};
use tracing::{debug, instrument};

/// Combines `Read` and `Seek` traits.
pub trait ReadAndSeek: Read + Seek {}

/// This trait is implemented for any type that implements both `Read` and `Seek`.
impl<R: Read + Seek> ReadAndSeek for R {}

/// Open a FASTA file with a FAI index.
///
/// Assumes that next to the FASTA file, there is a corresponding `.fai` index file.
/// These can be created using `samtools faidx` or `bgzip -r` commands.
#[instrument(level = "debug")]
pub fn open_fasta(path: &Path) -> Result<IndexedReader<Box<dyn ReadAndSeek>>, OpenFastaError> {
    let fasta_file = open_maybe_bgzip(path)?;
    let index_path = path.with_extension("fai");
    let fasta_index = bio::io::fasta::Index::from_file(&index_path).map_err(|err| {
        OpenFastaError::OpenFastaIndex { index_path, source: eyre!(Box::new(err)) }
    })?;
    Ok(IndexedReader::with_index(fasta_file, fasta_index))
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum OpenFastaError {
    #[error("Failed to open file {source}")]
    OpenFile {
        #[from]
        source: OpenMaybeBgzipError,
    },
    #[error("Failed to read index file `{index_path}`: {source}")]
    OpenFastaIndex { index_path: PathBuf, source: color_eyre::Report },
}

/// Open the file at the path. If the file has a `.gz` extension, it is assumed
/// to be a bgzip-compressed file and that there is a corresponding `.gz.gzi`
/// index file.
#[instrument(level = "debug")]
fn open_maybe_bgzip<P: AsRef<Path> + std::fmt::Debug>(
    path: P,
) -> Result<Box<dyn ReadAndSeek>, OpenMaybeBgzipError> {
    let path = path.as_ref();

    fn open(path: &Path) -> Result<Box<dyn ReadAndSeek>, OpenMaybeBgzipError> {
        let file = File::open(path)
            .map_err(|source| OpenMaybeBgzipError::OpenFile { path: path.to_owned(), source })?;
        let buffered = std::io::BufReader::new(file);
        Ok(Box::new(buffered))
    }

    if path.extension().unwrap_or_default() == "gz" {
        let mut index_path = Path::new(path).to_owned();
        index_path.set_extension("gz.gzi");
        if !index_path.exists() {
            return Err(OpenMaybeBgzipError::NoIndexFile { path: path.to_owned(), index_path });
        }
        let index = BGZFIndex::from_reader(open(&index_path)?).map_err(|source| {
            OpenMaybeBgzipError::ReadIndex { path: index_path.to_owned(), source }
        })?;
        let gzreader = BGZFReader::new(open(path)?)
            .map_err(|source| OpenMaybeBgzipError::ReadBgzip { path: path.to_owned(), source })?;
        let in_file = IndexedBGZFReader::new(gzreader, index).map_err(|source| {
            OpenMaybeBgzipError::IndexedBgzipReader { path: path.to_owned(), source }
        })?;
        debug!(?index_path, "Opened bgzip file");
        Ok(Box::new(in_file))
    } else {
        let in_file = open(path)?;
        debug!("Opened uncompressed file");
        Ok(Box::new(in_file))
    }
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
enum OpenMaybeBgzipError {
    #[error("Failed to open file `{path}`: {source}")]
    OpenFile { path: PathBuf, source: std::io::Error },
    #[error("{index_path:?} does not exist. bgzip input must be indexed (try bgzip -r {path:?})")]
    NoIndexFile { path: PathBuf, index_path: PathBuf },
    #[error("Failed to read index file `{path}`: {source}")]
    ReadIndex { path: PathBuf, source: std::io::Error },
    #[error("Failed to read bgzip file `{path}`: {source}")]
    ReadBgzip { path: PathBuf, source: bgzip::BGZFError },
    #[error("Failed to read indexed bgzip file `{path}`: {source}")]
    IndexedBgzipReader { path: PathBuf, source: bgzip::BGZFError },
}

#[cfg(test)]
mod tests {
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
        let mut gz = bgzip::BGZFWriter::with_compress_unit_size(
            File::create(&file_path)?,
            bgzip::Compression::fast(),
            16, // extra small chunk size so we get index entries!
            true,
        )?;
        gz.write_all(b"compressed content :)")?;
        let index = gz.close()?.expect("index");
        index.write(&mut File::create(&index_path)?)?;

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
}
