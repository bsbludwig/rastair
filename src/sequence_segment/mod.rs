use std::path::{Path, PathBuf};
use std::io::{Read, Seek};
use std::fs;
use log::{trace, debug, error};
use anyhow::{bail, Result};

use bio::io::fasta::{IndexedReader, Sequence};
use bio::utils::Text;

const DEFAULT_STEP_SIZE: u64 = 100000;

pub struct SequenceSegment {
    pub sequence: Text,
    pub contig: String,
    pub start:	u64,
    pub stop: u64,
}

pub struct SequenceSegmentIterator<R: Read + Seek> {
    reader: IndexedReader<R>,
    sequences: Vec<Sequence>,
    index_pos: usize,
    pos: u64,
    step_size: u64
}

impl SequenceSegmentIterator<fs::File> {
    pub fn with_file<P: AsRef<Path> + std::fmt::Debug>(fasta_path: P) -> Result<Self> {
        let reader = IndexedReader::from_file(&fasta_path)?;
        // Get the first contig in the index. If there is none, return an Error
        Self::with_reader(reader)
    }

    pub fn with_file_and_stepsize(fasta_path: &PathBuf, step_size: u64) -> Result<Self> {
        let mut new_seq_seg = Self::with_file(fasta_path)?;
        new_seq_seg.step_size = step_size;
        Ok(new_seq_seg)
    }
}

impl <R> SequenceSegmentIterator<R> 
where 
    R: Read + Seek
{
    pub fn with_reader(reader:IndexedReader<R>) -> Result<Self> {
        let sequences = reader.index.sequences();

        debug!("Read index with {} sequences", sequences.len());

        if sequences.is_empty() {
            bail!("No sequences in FASTA file");
        }
      
        let new_seq_seg = SequenceSegmentIterator {
            reader,
            sequences,
            index_pos: 0,
            pos: 0,
            step_size: DEFAULT_STEP_SIZE
        };
        Ok(new_seq_seg)
    }

    pub fn reached_end(&self) -> bool {
        self.index_pos >= self.sequences.len()
    }
}

impl <R> Iterator for SequenceSegmentIterator<R> 
where 
    R: Read + Seek
{
    type Item = SequenceSegment;
    fn next(&mut self) -> Option<Self::Item> {
        // Check if we've reached the end already
        if self.reached_end() {
            return None;
        }

        let start = self.pos;
        let seq_info = &self.sequences[self.index_pos];
        let stop = {
            if start + self.step_size > seq_info.len { 
                trace!("Reached end of {}, clipping to {}", &seq_info.name, seq_info.len);
                seq_info.len
            } else { 
                start + self.step_size 
            }
        };

        trace!("Moving cursor in fasta file to {}:{}-{}", &seq_info.name, start, stop);
        match self.reader.fetch(&seq_info.name, start, stop) {
            Err(e) => {
                error!("Error fetching sequence: {}", e);
                return None;
            },
            Ok(_) => {},
        };
        // Increment internal pointer
        if stop >= seq_info.len {
            trace!("Reached end of {}", &seq_info.name);
            self.index_pos = self.index_pos + 1;
            if !self.reached_end() {
                self.pos = 0;
            }
        } else {
            self.pos = stop;
        }
        // "Allocate" a sufficiently large chunk of memory
        let mut sequence: Text = vec![0 as u8; (stop-start) as usize];
        match self.reader.read(&mut sequence) {
            Err(e) => {
                error!("Error reading from fasta ({} from {} to {}): {}", &seq_info.name, start, stop, e);
                return None;
            },
            Ok(_)   => {}
        };
        let segment = SequenceSegment {
            sequence,
            contig: seq_info.name.clone(),
            start:  start,
            stop:   stop,
        };
        Some(segment)
    }
}

/***************************************************** 
 * Unit Tests
*****************************************************/
#[cfg(test)]
mod tests {
    // For testing
    use super::*;
    use std::io::Cursor;

    const FASTA_FILE: &[u8] = b">id desc
ACCGTAGGCTGA
CCGTAGGCTGAA
CGTAGGCTGAAA
GTAGGCTGAAAA
CCCC
>id2
ATTGTTGTTTTA
ATTGTTGTTTTA
ATTGTTGTTTTA
GGGG
";
    const FAI_FILE: &[u8] = b"id\t52\t9\t12\t13
id2\t40\t71\t12\t13
";

    fn new_seg_iter() -> Result<SequenceSegmentIterator<Cursor<&'static[u8]>>> {
        let reader = IndexedReader::new(Cursor::new(FASTA_FILE), FAI_FILE)?;
        let seg_iter = SequenceSegmentIterator::with_reader(reader)?;
        Ok(seg_iter)
    }

    #[test]
    fn can_create_segment_iterator() -> Result<()> {
        let seg_iter = new_seg_iter()?;

        assert_eq!(seg_iter.index_pos, 0);
        assert_eq!(seg_iter.pos, 0);
        Ok(())
    }
    
    #[test]
    fn can_check_end() -> Result<()>{
        let mut seg_iter = new_seg_iter()?;

        assert!(!seg_iter.reached_end());
        seg_iter.index_pos = 2;
        assert!(seg_iter.reached_end());
        Ok(())
    }

    #[test]
    fn can_iterate() -> Result<()> {
        let mut seg_iter = new_seg_iter()?;
        seg_iter.step_size = 30;
        let seq_info = seg_iter.next().unwrap();

        assert_eq!(&seq_info.contig, "id");
        assert_eq!(seq_info.start, 0);
        assert_eq!(seq_info.stop, 30);
        assert_eq!(&seq_info.sequence, b"ACCGTAGGCTGACCGTAGGCTGAACGTAGG");

        let seq_info2 = seg_iter.next().unwrap();
        assert_eq!(&seq_info2.contig, "id");
        assert_eq!(seq_info2.start, 30);
        assert_eq!(seq_info2.stop, 52);
        assert_eq!(&seq_info2.sequence, b"CTGAAAGTAGGCTGAAAACCCC");

        let seq_info3 = seg_iter.next().unwrap();
        assert_eq!(&seq_info3.contig, "id2");
        assert_eq!(seq_info3.start, 0);
        assert_eq!(seq_info3.stop, 30);
        assert_eq!(&seq_info3.sequence, b"ATTGTTGTTTTAATTGTTGTTTTAATTGTT");

        let seq_info4 = seg_iter.next().unwrap();
        assert_eq!(&seq_info4.contig, "id2");
        assert_eq!(seq_info4.start, 30);
        assert_eq!(seq_info4.stop, 40);
        assert_eq!(&seq_info4.sequence, b"GTTTTAGGGG");

        assert!(seg_iter.reached_end());

        assert!(seg_iter.next().is_none());

        Ok(())
    }
}