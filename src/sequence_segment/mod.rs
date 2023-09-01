use std::fmt::{Display, Formatter};
use std::path::Path;
use std::io::{Read, Seek};
use std::fs;
use log::{trace, debug, error};
use anyhow::{bail, Result};

use bio::io::fasta::IndexedReader;
use bio::utils::Text;

const DEFAULT_STEP_SIZE: u64 = 100000;

/// A genomic position, represented by its base and position, and the position
/// relative to the (arbitrary) segment slice it belongs to.
pub struct ContigPosition<'a>
{
    pub pos_in_segment: usize,
    segment: &'a SequenceSegment
}

impl<'a> ContigPosition<'a>
{
    /// The base at the represented position
    pub fn base(& self) -> u8
    {
        self.segment.sequence[self.pos_in_segment]
    }

    /// The position of the represented position in the overall contig
    pub fn pos_in_contig(& self) -> u64
    {
        (self.pos_in_segment as u64) + self.segment.start
    }
}

impl <'a> Display for ContigPosition<'a>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result
    {
        let char = char::from_u32(self.base() as u32).unwrap_or_default();
        write!(f, "{}:{} ({})", self.segment.contig, self.pos_in_contig(), char)
    }
}

/// A segment of (DNA) sequence with the associated sequence included
pub struct SequenceSegment
{
    pub sequence: Text,
    pub contig: String,
    pub start:	u64,
    pub stop: u64,
}

impl SequenceSegment
{
    /// Find the positions (relative to the contig coordinates) of all CpGs
    /// in the current segment
    pub fn find_cpgs(&self) -> Option<Vec<ContigPosition>>
    {
        let mut positions: Vec<ContigPosition> = Vec::new();
        if self.sequence.len() == 0
        {
            return None
        }

        for i in 1..(self.sequence.len())
        {
            let sub_seq = &self.sequence[(i-1)..(i+1)];
            if sub_seq == b"CG" {

                positions.push(ContigPosition{segment: self,
                                              pos_in_segment: i-1});
                positions.push(ContigPosition{segment: self,
                                              pos_in_segment: i});
            }
        }
        Some(positions)
    }
}

impl Display for SequenceSegment
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result
    {
        write!(f, "{}:{}-{}", self.contig, self.start, self.stop)
    }
}

/// An iterator over a fasta file that tokenises sequences
/// into shorter, more manageable chunks
pub struct SequenceSegmentIterator<R: Read + Seek>
{
    reader: IndexedReader<R>,
    sequences: Vec<(String, u64, u64)>,
    index_pos: usize,
    pos: u64,
    step_size: u64,
    tiling: usize,
}

impl SequenceSegmentIterator<fs::File>
{
    /// Create a new iterator with a reference to a FASTA file
    pub fn with_file<P: AsRef<Path> + std::fmt::Debug>(fasta_path: P) -> Result<Self>
    {
        let reader = IndexedReader::from_file(&fasta_path)?;
        // Get the first contig in the index. If there is none, return an Error
        Self::with_reader(reader)
    }

    pub fn with_file_and_stepsize<P: AsRef<Path> + std::fmt::Debug>(fasta_path: P, step_size: u64) -> Result<Self>
    {
        let mut new_seq_seg = Self::with_file(fasta_path)?;
        new_seq_seg.step_size = step_size;
        Ok(new_seq_seg)
    }
}

impl <R> SequenceSegmentIterator<R>
where
    R: Read + Seek
{
    pub fn with_reader(reader:IndexedReader<R>) -> Result<Self>
    {
        let sequences: Vec<(String, u64, u64)> = reader.index.sequences().iter().map(|seq| (seq.name.clone(), 0, seq.len)).collect();

        debug!("Read index with {} sequences", sequences.len());

        if sequences.is_empty()
        {
            bail!("No sequences in FASTA file");
        }

        let new_seq_seg = SequenceSegmentIterator
        {
            reader,
            sequences,
            index_pos: 0,
            pos: 0,
            step_size: DEFAULT_STEP_SIZE,
            tiling: 1
        };
        Ok(new_seq_seg)
    }

    fn reached_end(&self) -> bool
    {
        self.index_pos >= self.sequences.len()
    }

    /// Intersect the sequences in the fasta file with sequences from e.g. the bam index.
    /// This will reset the iterator, so the next call to `next()` will start form the leftmost interval
    /// in the first sequence again.
    pub fn subset_to_intervals(&mut self, intervals: &[(Vec<u8>, u64, u64)]) -> Option<()> 
    {
        let mut new_sequences: Vec<(String, u64, u64)> = Vec::with_capacity(self.sequences.len());
        for interval in intervals
        {
            let interval_id = std::str::from_utf8(&interval.0).unwrap_or_default();
            if let Some(sequence) = self.sequences.iter().find(|&seq| &seq.0 == &interval_id) 
            {
                new_sequences.push((sequence.0.clone(), std::cmp::max(sequence.1, interval.1), std::cmp::min(sequence.2, interval.2)));
            }
        }
        if new_sequences.len() == 0
        {
            return None;
        }
        self.sequences = new_sequences;
        self.rewind()
    }

    fn rewind(&mut self) -> Option<()> {
        self.index_pos = 0;
        if self.sequences.len() > 0 {
            let seq = &self.sequences[self.index_pos];
            self.pos = seq.1; // initialise as min start position
            return Some(());
        } 
        else 
        {
            self.pos = 0;
            return None;
        }
    }
}

impl <R> Iterator for SequenceSegmentIterator<R>
where
    R: Read + Seek
{
    type Item = SequenceSegment;
    fn next(&mut self) -> Option<Self::Item>
    {
        // Check if we've reached the end already
        if self.reached_end()
        {
            return None;
        }

        let start = self.pos;
        let seq_info = &self.sequences[self.index_pos];
        let stop =
        {
            if start + self.step_size > seq_info.2
            {
                trace!("Reached end of {}, clipping to {}", &seq_info.0, seq_info.2);
                seq_info.2
            }
            else
            {
                start + self.step_size
            }
        };

        trace!("Moving cursor in fasta file to {}:{}-{}", &seq_info.0, start, stop);
        match self.reader.fetch(&seq_info.0, start, stop)
        {
            Err(e) =>
            {
                error!("Error fetching sequence: {}", e);
                return None;
            },
            Ok(_) => {},
        };

        // Increment internal pointer
        if stop >= seq_info.2
        {
            trace!("Reached end of {}", &seq_info.0);
            self.index_pos = self.index_pos + 1;
            if !self.reached_end()
            {
                self.pos = 0;
            }
        }
        else
        {
            self.pos = stop - (self.tiling as u64);
        }

        // "Allocate" a sufficiently large chunk of memory
        let mut sequence: Text = vec![0 as u8; (stop-start) as usize];
        match self.reader.read(&mut sequence)
        {
            Err(e) =>
            {
                error!("Error reading from fasta ({} from {} to {}): {}", &seq_info.0, start, stop, e);
                return None;
            },
            Ok(_)   => {}
        };
        let segment = SequenceSegment
        {
            sequence,
            contig: seq_info.0.clone(),
            start:  start,
            stop:   stop,
        };
        Some(segment)
    }
}

/*====================================================
 = Unit Tests
====================================================*/
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

    fn new_seg_iter() -> Result<SequenceSegmentIterator<Cursor<&'static[u8]>>>
    {
        let reader = IndexedReader::new(Cursor::new(FASTA_FILE), FAI_FILE)?;
        let seg_iter = SequenceSegmentIterator::with_reader(reader)?;
        Ok(seg_iter)
    }

    #[test]
    fn can_create_segment_iterator() -> Result<()>
    {
        let seg_iter = new_seg_iter()?;

        assert_eq!(seg_iter.index_pos, 0);
        assert_eq!(seg_iter.pos, 0);
        Ok(())
    }

    #[test]
    fn can_check_end() -> Result<()>
    {
        let mut seg_iter = new_seg_iter()?;

        assert!(!seg_iter.reached_end());
        seg_iter.index_pos = 2;
        assert!(seg_iter.reached_end());
        Ok(())
    }

    #[test]
    fn can_iterate() -> Result<()>
    {
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

    #[test]
    fn can_find_cpgs() -> Result<()>
    {
        let mut seg_iter = new_seg_iter()?;
        seg_iter.step_size = 12;
        let seq_info = seg_iter.next().unwrap();

        let cpg_pos = seq_info.find_cpgs().unwrap();
        assert_eq!(cpg_pos.len(), 2);
        assert_eq!(cpg_pos[0].pos_in_contig(), 2);
        assert_eq!(cpg_pos[1].pos_in_contig(), 3);
        assert_eq!(cpg_pos[0].base(), b'C');
        assert_eq!(cpg_pos[1].base(), b'G');

        let seq_info2 = seg_iter.next().unwrap();
        let cpg_pos2 = seq_info2.find_cpgs().unwrap();
        assert_eq!(cpg_pos2.len(), 2);
        assert_eq!(cpg_pos2[0].pos_in_contig(), 13);
        assert_eq!(cpg_pos2[0].pos_in_segment, 1);
        Ok(())
    }
}