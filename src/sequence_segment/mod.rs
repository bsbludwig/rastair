use std::path::PathBuf;

use bio::io::fasta::{IndexedReader, Sequence};
use bio::utils::Text;

use log::{trace, debug, error};
use std::fs;

use anyhow::{bail, Result};

const DEFAULT_STEP_SIZE: u64 = 100000;

pub struct SequenceSegment {
    pub sequence: Text,
    pub contig: String,
    pub start:	u64,
    pub stop: u64,
}

pub struct SequenceSegmentIterator {
    reader: IndexedReader<fs::File>,
    sequences: Vec<Sequence>,
    index_pos: usize,
    pos: u64,
    step_size: u64
}

impl SequenceSegmentIterator{
    pub fn from_file(fasta_path: &PathBuf) -> Result<Self> {
        let reader = IndexedReader::from_file(fasta_path)?;
        // Get the first contig in the index. If there is none, return an Error
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

    pub fn from_file_with_stepsize(fasta_path: &PathBuf, step_size: u64) -> Result<Self> {
        let mut new_seq_seg = Self::from_file(fasta_path)?;
        new_seq_seg.step_size = step_size;
        Ok(new_seq_seg)
    }

    pub fn reached_end(&self) -> bool {
        self.index_pos >= self.sequences.len()
    }
}

impl Iterator for SequenceSegmentIterator {
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