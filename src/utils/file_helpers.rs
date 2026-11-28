use color_eyre::{
    Section,
    eyre::{Context, Result},
};
use seqair::fasta::IndexedFastaReader;
use seqair_types::Pos0;
use std::{
    collections::HashSet,
    fmt,
    path::{Path, PathBuf},
};
use tracing::instrument;

pub struct FastaReader {
    fasta_file: PathBuf,
    reader: IndexedFastaReader,
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for FastaReader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FastaReader").field("fasta", &self.fasta_file).finish()
    }
}

impl FastaReader {
    /// Names of every sequence in the FASTA index.
    pub fn sequence_names(&self) -> HashSet<String> {
        self.reader.index().sequence_names().into_iter().map(String::from).collect()
    }

    /// The bases of `name` over the 0-based half-open `[start, stop)`, in
    /// upper case; nothing for an empty stretch.
    pub fn fetch_seq(&mut self, name: &str, start: u64, stop: u64) -> Result<Vec<u8>> {
        let pos = |at: u64| u32::try_from(at).ok().and_then(Pos0::new);
        let last = stop.checked_sub(1).filter(|last| *last >= start);
        let (Some(first), Some(last)) = (pos(start), last.and_then(pos)) else {
            return Ok(Vec::new());
        };
        self.reader.fetch_seq(name, (first..=last).into()).wrap_err_with(|| {
            format!("Failed to fetch sequence {name}:{start}-{stop} from {:?}", self.fasta_file)
        })
    }
}

/// Open a FASTA file with its FAI index (and GZI index, if it is
/// bgzip-compressed).
#[instrument(level = "debug", skip_all)]
pub fn open_fasta(fasta_path: &Path) -> Result<FastaReader> {
    let reader = IndexedFastaReader::open(fasta_path)
        .wrap_err_with(|| format!("Failed to open FASTA file {fasta_path:?}"))
        .with_suggestion(|| {
            format!(
                "Ensure the FASTA file exists and has an index: `samtools faidx {fasta_path:?}`"
            )
        })?;

    Ok(FastaReader { fasta_file: fasta_path.to_path_buf(), reader })
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
        let index_path = dir.path().join("test.fa.fai");

        // Create a simple FASTA file
        let mut fasta_file = std::fs::File::create(&fasta_path)?;
        fasta_file.write_all(b">seq1\nACGT\n>seq2\nGTCA\n")?;

        // Create a simple FAI index
        let mut index_file = std::fs::File::create(&index_path)?;
        // Format: name, length, offset, line_bases, line_width
        index_file.write_all(b"seq1\t4\t6\t4\t5\nseq2\t4\t16\t4\t5\n")?;

        let mut reader = open_fasta(&fasta_path)?;
        let seq = reader.fetch_seq("seq1", 0, 4)?;
        assert_eq!(seq, b"ACGT");

        Ok(())
    }

    #[test]
    fn test_open_fasta_missing_file() {
        let non_existent_path = PathBuf::from("/path/that/does/not/exist.fa");
        let result = open_fasta(&non_existent_path);
        assert!(result.is_err());
    }

    #[test]
    fn test_open_fasta_fetch_region() -> Result<()> {
        let dir = TempDir::new()?;
        let fasta_path = dir.path().join("test.fa");
        let index_path = dir.path().join("test.fa.fai");

        let mut fasta_file = std::fs::File::create(&fasta_path)?;
        fasta_file.write_all(b">seq1\nACGTACGT\n")?;

        let mut index_file = std::fs::File::create(&index_path)?;
        index_file.write_all(b"seq1\t8\t6\t8\t9\n")?;

        let mut reader = open_fasta(&fasta_path)?;

        // Fetch a sub-region using 0-based half-open coords
        let seq = reader.fetch_seq("seq1", 2, 6)?;
        assert_eq!(seq, b"GTAC");

        Ok(())
    }

    /// The index is not built on the fly: a missing one is an error that
    /// says how to make it.
    #[test]
    fn test_open_fasta_missing_index_is_an_error() -> Result<()> {
        let dir = TempDir::new()?;
        let fasta_path = dir.path().join("test.fa");

        let mut fasta_file = std::fs::File::create(&fasta_path)?;
        fasta_file.write_all(b">seq1\nACGT\n")?;

        let Err(error) = open_fasta(&fasta_path) else {
            color_eyre::eyre::bail!("a FASTA without an index was opened");
        };
        assert!(format!("{error:?}").contains("samtools faidx"), "{error:?}");
        Ok(())
    }

    #[test]
    fn test_open_fasta_with_different_index_names() -> Result<()> {
        let dir = TempDir::new()?;
        let fasta_path = dir.path().join("test.fasta");
        let index_path = dir.path().join("test.fasta.fai");

        let mut fasta_file = std::fs::File::create(&fasta_path)?;
        fasta_file.write_all(b">seq1\nACGT\n>seq2\nGTCA\n")?;

        let mut index_file = std::fs::File::create(&index_path)?;
        index_file.write_all(b"seq1\t4\t6\t4\t5\nseq2\t4\t16\t4\t5\n")?;

        open_fasta(&fasta_path)?;

        Ok(())
    }

    #[test]
    fn test_open_fasta_with_bgzip() -> Result<()> {
        let dir = TempDir::new()?;
        let fasta_path = dir.path().join("test.fa.gz");
        let fasta_index_path = dir.path().join("test.fa.gz.fai");
        let bgzip_index_path = dir.path().join("test.fa.gz.gzi");

        // Create a bgzipped FASTA file with gzi index
        let mut gz = noodles_bgzf::io::Writer::new(std::fs::File::create(&fasta_path)?);
        gz.write_all(b">seq1\nACGT\n>seq2\nGTCA\n")?;
        gz.finish()?;
        noodles_bgzf::gzi::fs::write(&bgzip_index_path, &noodles_bgzf::gzi::Index::default())?;

        // Create FAI index (htslib expects .fa.gz.fai for bgzipped FASTA)
        let mut index_file = std::fs::File::create(&fasta_index_path)?;
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
        let mut fasta_file = std::fs::File::create(&fasta_path)?;
        fasta_file.write_all(b">seq1\nACGT\n")?;

        // Create corrupt FAI index
        let mut index_file = std::fs::File::create(&index_path)?;
        index_file.write_all(b"corrupted_content")?;

        let result = open_fasta(&fasta_path);
        assert!(result.is_err());

        Ok(())
    }

    #[test]
    fn test_fetch_nonexistent_sequence() -> Result<()> {
        let dir = TempDir::new()?;
        let fasta_path = dir.path().join("test.fa");
        let index_path = dir.path().join("test.fa.fai");

        let mut fasta_file = std::fs::File::create(&fasta_path)?;
        fasta_file.write_all(b">seq1\nACGT\n")?;

        let mut index_file = std::fs::File::create(&index_path)?;
        index_file.write_all(b"seq1\t4\t6\t4\t5\n")?;

        let mut reader = open_fasta(&fasta_path)?;
        let result = reader.fetch_seq("nonexistent", 0, 4);
        assert!(result.is_err());

        Ok(())
    }
}
