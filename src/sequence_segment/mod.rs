use std::fmt::{Display, Formatter};
use std::path::{Path, PathBuf};
use std::io::{stdout, Write, Read, Seek};
use std::fs;
use log::{trace, debug, error};

use anyhow::{bail, Result, anyhow};
// Very fast substr search
use memchr::memmem::find_iter;

use bio::io::fasta::IndexedReader;
use bio::utils::Text;
use rust_htslib::bam::FetchDefinition;

use crate::utils::FetchDefinitionExt;

const DEFAULT_STEP_SIZE: usize = 10000;
const DEFAULT_TILING: usize = 1;

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
    pub is_last: bool,
}

impl SequenceSegment
{
    /// Find the positions (relative to the contig coordinates) of all CpGs
    /// in the current segment. Return None if no CpGs are found
    pub fn find_cpgs(&self) -> Option<Vec<ContigPosition<'_>>>
    {
        if self.sequence.len() == 0
        {
            return None
        }
        let start_positions: Vec<usize> =
            find_iter(&self.sequence, b"CG")
            .collect();
        let results: Vec<ContigPosition> = start_positions
            .iter()
            .map(|pos|
                vec![ContigPosition { pos_in_segment: *pos, segment: self },
                     ContigPosition { pos_in_segment: *pos+1, segment: self }])
            .flatten()
            .collect();
        match results.len() {
            0   => None,
            _   => Some(results)
        }
    }
}

impl Display for SequenceSegment
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result
    {
        write!(f, "{}:{}-{}", self.contig, self.start, self.stop)
    }
}

pub struct GenomicRegion
{
    contig: String,
    start: u64,
    end: u64
}
impl Clone for GenomicRegion
{
    fn clone(&self) -> Self {
        GenomicRegion {
            contig: self.contig.clone(),
            start: self.start,
            end: self.end
        }
    }
}
/// An iterator over a fasta file that tokenises sequences
/// into shorter, more manageable chunks
pub struct SequenceSegmentIterator<R: Read + Seek>
{
    /// FASTA reader
    reader: IndexedReader<R>,
    /// Array of target regions to iterate over (contigId/from/to)
    sequences: Vec<GenomicRegion>,
    /// current position in the `sequences`` index
    index_pos: usize,
    /// position in current contig
    pos: u64,
    /// size of the interval to process at each step
    step_size: usize,
    /// overlap of the next window with the previous one.
    /// Defaults to 1, ie include the last base of the previous window in the same contig
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

    pub fn with_file_and_stepsize<P: AsRef<Path> + std::fmt::Debug>(fasta_path: P, step_size: usize) -> Result<Self>
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
        let sequences: Vec<GenomicRegion> =
            reader
                .index
                .sequences()
                .iter()
                .map(|seq| GenomicRegion {contig: seq.name.clone(), start: 0, end: seq.len})
                .collect();

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
            tiling: DEFAULT_TILING
        };

        Ok(new_seq_seg)
    }

    fn reached_end(&self) -> bool
    {
        self.index_pos >= self.sequences.len()
    }
    pub fn set_tiling(&mut self, new_tiling: usize) -> Result<()>
    {
        if new_tiling == 0 || new_tiling > self.step_size
        {
            return Err(anyhow!("Incorrect tiling setting"));
        }
        self.tiling = new_tiling;
        Ok(())
    }

    pub fn subset_to_region(&mut self, region: &String) -> Result<()>
    {
        match FetchDefinition::from_region_string(region)? {
            FetchDefinition::RegionString(chr_bytes, start, end) => {
                let chr = std::str::from_utf8(chr_bytes)?;
                if let Some(sequence) = self.sequences.iter().find(|&seq| &seq.contig == chr)
                {
                    self.sequences = Vec::from([GenomicRegion{ contig: chr.to_owned(),
                                                               start: std::cmp::max(sequence.start, start as u64),
                                                               end: std::cmp::min(sequence.end, end as u64)
                                                             }
                                                ]);
                }
                else
                {
                    bail!("Could not find {}", region);
                }
            },
            FetchDefinition::String(chr_bytes) => {
                let chr = std::str::from_utf8(chr_bytes)?;
                if let Some(sequence) = self.sequences.iter().find(|&seq| &seq.contig == chr)
                {
                    self.sequences = Vec::from([sequence.clone()]);
                }
                else
                {
                    bail!("Could not find {}", region);
                }
            },
            FetchDefinition::All => return Ok(()),
            _   => bail!("Error subsetting to region: {}", region)
        };
        Ok(())
    }

    /// Intersect the sequences in the fasta file with sequences from e.g. the bam index.
    /// This will reset the iterator, so the next call to `next()` will start form the leftmost interval
    /// in the first sequence again.
    pub fn subset_to_intervals(&mut self, intervals: &[(Vec<u8>, u64, u64)]) -> Result<()>
    {
        if intervals.len() == 0
        {
            bail!("Empty intervals provided");
        }

        let mut new_sequences: Vec<GenomicRegion> = Vec::with_capacity(intervals.len());
        for interval in intervals
        {
            let interval_id = std::str::from_utf8(&interval.0).unwrap_or_default();
            if let Some(sequence) = self.sequences.iter().find(|&seq| &seq.contig == &interval_id)
            {
                new_sequences.push(GenomicRegion{
                    contig: sequence.contig.clone(),
                    start: std::cmp::max(sequence.start, interval.1),
                    end: std::cmp::min(sequence.end, interval.2)
                });
            }
        }
        if new_sequences.len() == 0
        {
            bail!("No overlap between fasta and intervals");
        }
        self.sequences = new_sequences;
        self.rewind()
    }
    fn rewind(&mut self) -> Result<()> {
        self.index_pos = 0;
        if self.sequences.len() > 0 {
            let seq = &self.sequences[self.index_pos];
            self.pos = seq.start; // initialise as min start position
        }
        else
        {
            self.pos = 0;
        }
        Ok(())
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

        let mut start = self.pos;
        if start >= self.tiling as u64
        {
            start = start - self.tiling as u64;
        }
        let seq_info = &self.sequences[self.index_pos];
        let stop =
        {
            if start + (self.step_size as u64) > seq_info.end
            {
                trace!("Reached end of {}, clipping to {}", &seq_info.contig, seq_info.end);
                seq_info.end
            }
            else
            {
                start + (self.step_size as u64)
            }
        };

        trace!("Moving cursor in fasta file to {}:{}-{}", &seq_info.contig, start, stop);
        match self.reader.fetch(&seq_info.contig, start, stop)
        {
            Err(e) =>
            {
                error!("Error fetching sequence: {}", e);
                return None;
            },
            Ok(_) => {},
        };

        // Increment internal pointer
        if stop >= seq_info.end
        {
            trace!("Reached end of {}", &seq_info.contig);
            self.index_pos = self.index_pos + 1;
            if !self.reached_end()
            {
                self.pos = 0;
            }
        }
        else
        {
            self.pos = stop;
        }

        // "Allocate" a sufficiently large chunk of memory
        let mut sequence: Text = vec![0 as u8; (stop-start) as usize];
        match self.reader.read(&mut sequence)
        {
            Err(e) =>
            {
                error!("Error reading from fasta ({} from {} to {}): {}", &seq_info.contig, start, stop, e);
                return None;
            },
            Ok(_)   => {}
        };
        // Make sequence all uppercase
        sequence = sequence.to_ascii_uppercase();
        let segment = SequenceSegment
        {
            sequence,
            contig: seq_info.contig.clone(),
            start:  start,
            stop:   stop,
            is_last: stop >= seq_info.end
        };
        Some(segment)
    }
}

pub fn run_finder (file_path: &PathBuf, step_size: usize) -> Result<()>
{
    let seg_iter = if step_size > 0
    {
        debug!("Using step size {}", step_size);
        SequenceSegmentIterator::with_file_and_stepsize(file_path, step_size)?
    }
    else
    {
        SequenceSegmentIterator::with_file(file_path)?
    };
    let mut lock = stdout().lock();
    for segment in seg_iter
    {
        segment
            .find_cpgs()
            .unwrap_or_default()
            .iter()
            .for_each(|position|
                {
                    let char = char::from_u32(position.base() as u32).unwrap_or_default();
                    let seg_id = &position.segment.contig;
                    writeln!(lock, "{}\t{}\t{}", seg_id, position.pos_in_contig(), char).expect("Error writing to stdout");
                });
    }
    Ok(())
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
>id3
ATCGATCGATCG
AATCGATCGATC
GATCGATCGATC
GGGGG
>id4
ATCGATCGATcG
AATCGATCgATC
gATCGATcGATc
gGGcg
";
    const FAI_FILE: &[u8] = b"id\t52\t9\t12\t13
id2\t40\t71\t12\t13
id3\t41\t120\t12\t13
id4\t41\t170\t12\t13
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
        seg_iter.index_pos = 4;
        assert!(seg_iter.reached_end());
        assert!(seg_iter.next().is_none());
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
        assert_eq!(seq_info.stop, seq_info.start + seq_info.sequence.len() as u64);

        let seq_info2 = seg_iter.next().unwrap();
        assert_eq!(&seq_info2.contig, "id");
        assert_eq!(seq_info2.start, 29);
        assert_eq!(seq_info2.stop, 52);
        assert_eq!(&seq_info2.sequence, b"GCTGAAAGTAGGCTGAAAACCCC");
        assert_eq!(seq_info2.stop, seq_info2.start + seq_info2.sequence.len() as u64);

        let seq_info3 = seg_iter.next().unwrap();
        assert_eq!(&seq_info3.contig, "id2");
        assert_eq!(seq_info3.start, 0);
        assert_eq!(seq_info3.stop, 30);
        assert_eq!(&seq_info3.sequence, b"ATTGTTGTTTTAATTGTTGTTTTAATTGTT");
        assert_eq!(seq_info3.stop, seq_info3.start + seq_info3.sequence.len() as u64);

        let seq_info4 = seg_iter.next().unwrap();
        assert_eq!(&seq_info4.contig, "id2");
        assert_eq!(seq_info4.start, 29);
        assert_eq!(seq_info4.stop, 40);
        assert_eq!(&seq_info4.sequence, b"TGTTTTAGGGG");

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
        assert_eq!(cpg_pos2[0].pos_in_segment, 2);
        Ok(())
    }

    #[test]
    fn can_find_split_cpgs() -> Result<()>
    {
        let mut seg_iter = new_seg_iter()?;
        seg_iter.step_size = 100000;
        let seq_info: Vec<SequenceSegment> = seg_iter.collect();
        assert_eq!(seq_info.len(), 4);
        let cpg_pos = seq_info[2].find_cpgs().unwrap();
        assert_eq!(cpg_pos.len(), 18);
        Ok(())
    }

    #[test]
    fn can_find_mixed_case() -> Result<()>
    {
        let mut seg_iter = new_seg_iter()?;
        seg_iter.step_size = 100000;
        let seq_info: Vec<SequenceSegment> = seg_iter.collect();
        assert_eq!(seq_info.len(), 4);
        let cpg_pos = seq_info[3].find_cpgs().unwrap();
        assert_eq!(cpg_pos.len(), 20);
        Ok(())
    }
}