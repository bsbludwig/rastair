use std::{io::{Read, Seek}, path::Path, fmt::Debug, fs::File};
use anyhow::{bail, Result};
use bgzip::{index::BGZFIndex, read::IndexedBGZFReader, BGZFReader};

pub trait ReadAndSeek: Read+Seek {}
impl <R: Read+Seek> ReadAndSeek for R {}

pub fn open_file<P: AsRef<Path> + Debug>(path: P) -> Result<Box<dyn ReadAndSeek>>
{
    if path.as_ref().extension().unwrap_or_default() == "gz"
    {
        let mut index_path = Path::new(path.as_ref()).to_owned();
        index_path.set_extension("gz.gzi");
        if !index_path.exists()
        {
            bail!("{} does not exist. bgzip input must be indexed (try bgzip -r {})", index_path.to_str().unwrap_or_default(), path.as_ref().to_str().unwrap_or_default());
        }
        let index = BGZFIndex::from_reader(File::open(index_path)?)?;
        let gzreader = BGZFReader::new(File::open(path)?)?;
        let in_file = IndexedBGZFReader::new(gzreader, index)?;
        Ok(Box::new(in_file))
    }
    else
    {
        let in_file = File::open(path)?;
        Ok(Box::new(in_file))
    }
}
