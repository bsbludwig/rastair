use bgzip::{BGZFReader, index::BGZFIndex, read::IndexedBGZFReader};
use color_eyre::eyre::{Result, bail};
use std::{
    fs::File,
    io::{Read, Seek},
    path::Path,
};
use tracing::{debug, instrument};

pub trait ReadAndSeek: Read + Seek {}
impl<R: Read + Seek> ReadAndSeek for R {}

#[instrument(level = "debug")]
pub fn open_maybe_bgzip<P: AsRef<Path> + std::fmt::Debug>(path: P) -> Result<Box<dyn ReadAndSeek>> {
    if path.as_ref().extension().unwrap_or_default() == "gz" {
        let mut index_path = Path::new(path.as_ref()).to_owned();
        index_path.set_extension("gz.gzi");
        if !index_path.exists() {
            bail!(
                "{:?} does not exist. bgzip input must be indexed (try bgzip -r {:?})",
                index_path,
                path.as_ref()
            );
        }
        let index = BGZFIndex::from_reader(File::open(&index_path)?)?;
        let gzreader = BGZFReader::new(File::open(&path)?)?;
        let in_file = IndexedBGZFReader::new(gzreader, index)?;
        debug!(?index_path, "Opened bgzip file");
        Ok(Box::new(in_file))
    } else {
        let in_file = File::open(&path)?;
        debug!("Opened uncompressed file");
        Ok(Box::new(in_file))
    }
}
