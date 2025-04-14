use bgzip::{BGZFReader, index::BGZFIndex, read::IndexedBGZFReader};
use bio::io::fasta::IndexedReader;
use color_eyre::eyre::{Result, eyre};
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

#[cfg(test)]
mod fasta_tests {
    use super::*;
    use std::io::Write;
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

        assert!(matches!(result, Err(OpenFastaError::OpenFile { .. })));
    }

    #[test]
    fn test_open_fasta_missing_index() -> Result<()> {
        let dir = TempDir::new()?;
        let fasta_path = dir.path().join("test.fa");

        // Create FASTA file without index
        let mut fasta_file = File::create(&fasta_path)?;
        fasta_file.write_all(b">seq1\nACGT\n")?;

        let result = open_fasta(&fasta_path);

        assert!(matches!(result, Err(OpenFastaError::OpenFastaIndex { .. })));

        Ok(())
    }

    #[test]
    fn test_open_fasta_with_bgzip() -> Result<()> {
        let dir = TempDir::new()?;
        let fasta_path = dir.path().join("test.fa.gz");
        let fasta_index_path = dir.path().join("test.fa.fai");
        let bgzip_index_path = dir.path().join("test.fa.gz.gzi");

        // Create a bgzipped FASTA file with gzi index
        let mut gz = bgzip::BGZFWriter::with_compress_unit_size(
            File::create(&fasta_path)?,
            bgzip::Compression::fast(),
            16, // extra small chunk size so we get index entries!
            true,
        )?;
        gz.write_all(b">seq1\nACGT\n>seq2\nGTCA\n")?;
        let index = gz.close()?.expect("index");
        index.write(&mut File::create(&bgzip_index_path)?)?;

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

        assert!(matches!(result, Err(OpenFastaError::OpenFastaIndex { .. })));

        Ok(())
    }
}

/// Open the file at the path. If the file has a `.gz` extension, it is assumed
/// to be a bgzip-compressed file and that there is a corresponding `.gz.gzi`
/// index file.
#[instrument(level = "debug")]
fn open_maybe_bgzip<P: AsRef<Path> + std::fmt::Debug>(
    path: P,
) -> Result<Box<dyn ReadAndSeek>, OpenMaybeBgzipError> {
    let path = path.as_ref();

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

fn open(path: &Path) -> Result<Box<dyn ReadAndSeek>, OpenMaybeBgzipError> {
    let file = File::open(path)
        .map_err(|source| OpenMaybeBgzipError::OpenFile { path: path.to_owned(), source })?;
    let buffered = std::io::BufReader::new(file);
    Ok(Box::new(buffered))
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum OpenMaybeBgzipError {
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

    #[test]
    fn test_open_file_error() {
        let path = PathBuf::from("/nonexistent/file.txt");
        let result = open_maybe_bgzip(&path);

        assert!(matches!(result, Err(OpenMaybeBgzipError::OpenFile { path: p, .. }) if p == path));
    }

    #[test]
    fn test_no_index_file_error() -> Result<()> {
        let dir = TempDir::new()?;
        let file_path = dir.path().join("test.gz");

        // Create a file with .gz extension but no index
        let mut file = File::create(&file_path)?;
        file.write_all(b"some content")?;

        let result = open_maybe_bgzip(&file_path);

        assert!(matches!(
            result,
            Err(OpenMaybeBgzipError::NoIndexFile { path, index_path })
            if path == file_path && index_path == file_path.with_extension("gz.gzi")
        ));

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

        assert!(matches!(
            result,
            Err(OpenMaybeBgzipError::ReadIndex { path, .. })
            if path == index_path
        ));

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

        // Create a dummy index file so we get past the index check
        let index = BGZFIndex::default();
        index.write(&mut File::create(&index_path)?)?;

        let result = open_maybe_bgzip(&file_path);

        assert!(matches!(
            result,
            Err(OpenMaybeBgzipError::ReadBgzip { path, .. })
            if path == file_path
        ));

        Ok(())
    }

    #[test]
    fn test_indexed_bgzip_reader_error() -> Result<()> {
        let dir = TempDir::new()?;
        let file_path = dir.path().join("test.gz");
        let index_path = dir.path().join("test.gz.gzi");

        // Create a valid BGZIP file
        let mut gz = bgzip::BGZFWriter::new(File::create(&file_path)?, bgzip::Compression::fast());
        gz.write_all(b"content")?;
        gz.close()?;

        // Create an index that doesn't match the file
        let mut index = File::create(&index_path)?;
        write!(&mut index, "no thanks")?;

        let result = open_maybe_bgzip(&file_path);

        assert!(matches!(result, Err(OpenMaybeBgzipError::ReadIndex { .. })));

        Ok(())
    }
}
