use bio::io::fasta::IndexedReader;
use bio::bio_types::strand::Strand;

use anyhow::Result;

use std::collections::BTreeMap;
use std::io::{Read, Seek};

use crate::sequence_segment::SequenceSegmentIterator;
pub struct CpgInfo
{
    pub contig  : Vec<u8>,
    pub position: u64,
    pub index   : usize
}

impl CpgInfo {
    pub fn new(contig: Vec<u8>, position: u64, index: usize) -> Self
    {
        CpgInfo{contig, position, index}
    }

    pub fn strand(&self) -> Strand
    {
        if self.index % 2 == 1
        {
            Strand::Reverse
        }
        else
        {
            Strand::Forward
        }
    }
}

pub struct CpgBuffer<R: Read + Seek>
{
    // buffer CpG positions and indices to get info on CpGs between read pairs
    cpg_buffer: BTreeMap<u64, CpgInfo>,
    fasta_reader: SequenceSegmentIterator<R>,
    last_in_segment: bool
}

impl <R> CpgBuffer<R>
where R: Read+Seek
{
    pub fn with_reader(fasta_reader: IndexedReader<R>) -> Result<Self>
    {
        let cpg_buffer = BTreeMap::new();
        let ssi = SequenceSegmentIterator::with_reader(fasta_reader)?;
        Ok(
            Self { cpg_buffer, fasta_reader: ssi, last_in_segment: false }
        )
    }

    pub fn with_reader_and_stepsize(fasta_reader: IndexedReader<R>, step_size: usize) -> Result<Self>
    {
        let cpg_buffer = BTreeMap::new();
        let ssi = SequenceSegmentIterator::with_reader_and_stepsize(fasta_reader, step_size)?;
        Ok(
            Self { cpg_buffer, fasta_reader: ssi, last_in_segment: false }
        )
    }

    fn parse_next_from_file(&mut self) -> Option<()>
    {
        if self.last_in_segment
        {
            return None; // don't progress until explicitly called "progress"
        }

        let last_index =
            if let Some((_, last_cpg)) = self.cpg_buffer.last_key_value()
            {
                last_cpg.index + 1
            }
            else
            {
                0
            };

        let segment = self.fasta_reader.next()?;
        if segment.is_last_in_contig
        {
            self.last_in_segment = true;
        }
        let mut iter: usize = 0;
        for cpg in segment.find_cpgs().unwrap_or_default()
        {
            self.cpg_buffer.insert(cpg.pos_in_contig(), CpgInfo::new(Vec::from(cpg.contig()), cpg.pos_in_contig(), last_index + iter));
            iter += 1;
        }
        Some(())
    }

    pub fn progress_to_contig(&mut self, chr: &[u8]) -> Option<()>
    {
        self.cpg_buffer.clear();
        self.last_in_segment = false;
        self.fasta_reader.move_to_contig(chr).ok()
    }

    pub fn cpgs_in_range(&mut self, chr: &[u8], start: u64, end: u64) -> Option<Vec<&CpgInfo>>
    {
        assert!(end > start, "Inverted slice not allowed");

        if let Some((_, last_cpg)) = self.cpg_buffer.last_key_value()
        {
            if last_cpg.contig != chr
            {
                self.progress_to_contig(chr);
            }
        }

        loop
        {
            if let Some((_, last_cpg)) = self.cpg_buffer.last_key_value()
            {
                if last_cpg.position + 1 < end
                {
                    match self.parse_next_from_file() {
                        None => break,
                        _   => ()
                    }
                }
                else
                {
                    break;
                }
            }
            else
            {
                match self.parse_next_from_file() {
                    None => break,
                    _   => ()
                }
            }
        }
        let slice = self.cpg_buffer.range(start..end).map(|(_, v)| v);
        Some(slice.collect())
    }
}

/*====================================================
 = Unit Tests
====================================================*/
//#[cfg(test)]
// mod tests {
//     use std::io::Cursor;

//     use super::*;
//     use anyhow::{Ok, Result};

//     const COORD_BED: &[u8] = b"chr1\t0\t1\t0\t+
// chr1\t1\t2\t1\t-
// chr1\t6\t7\t2\t+
// chr1\t7\t8\t3\t-
// chr1\t10\t11\t4\t+
// chr1\t11\t12\t5\t-
// chr1\t20\t21\t6\t+
// chr1\t21\t22\t7\t-
// chr1\t100\t101\t8\t+
// chr1\t101\t102\t9\t-
// chr1\t110\t111\t10\t+
// chr1\t111\t112\t11\t-
// chr2\t10\t11\t0\t+
// chr2\t11\t12\t1\t-
// chr2\t16\t17\t2\t+
// chr2\t17\t18\t3\t-
// chr2\t18\t19\t4\t+
// chr2\t19\t20\t5\t-
// ";

//     #[test]
//     fn can_create_buffer() -> Result<()>
//     {
//         let mut buffer = CpgBuffer::with_file(Cursor::new(COORD_BED))?;
//         let mut next_row = buffer.parse_next_from_file().expect("No next row");
//         assert_eq!(next_row.contig, Vec::from("chr1".as_bytes()));
//         assert_eq!(next_row.position, 0);
//         assert_eq!(next_row.index, 0);
//         next_row = buffer.parse_next_from_file().expect("No next row");

//         assert_eq!(next_row.contig, Vec::from("chr1".as_bytes()));
//         assert_eq!(next_row.position, 1);
//         assert_eq!(next_row.index, 1);
//         Ok(())
//     }

//     #[test]
//     fn can_empty_buffer() -> Result<()>
//     {
//         let mut buffer = CpgBuffer::with_file(Cursor::new(COORD_BED))?;
//         for _i in 0..18
//         {
//             buffer.parse_next_from_file().expect("No next row");
//         }
//         assert!(buffer.parse_next_from_file().is_none(), "More data than expected");
//         Ok(())
//     }

//     #[test]
//     fn can_load_until() -> Result<()>
//     {
//         let mut buffer = CpgBuffer::with_file(Cursor::new(COORD_BED))?;
//         buffer.load_data_until("chr1".as_bytes(), 50).expect("Failed to read slice");
//         assert!(buffer.cpg_buffer.len() >= 8, "Buffer too short, must be missing data");
//         Ok(())
//     }

//     #[test]
//     fn can_get_slice() -> Result<()>
//     {
//         let mut buffer = CpgBuffer::with_file(Cursor::new(COORD_BED))?;
//         let mut slice = buffer.cpgs_in_range("chr1".as_bytes().as_ref(), 0, 50).expect("Empty slice");

//         assert_eq!(slice.len(), 8, "Slice length not right");
//         assert_eq!(slice[7].index, 7, "Element index doesn't match");

//         // Shorter slice from a middle start point
//         slice = buffer.cpgs_in_range("chr1".as_bytes().as_ref(), 15, 50).expect("Empty slice");
//         assert_eq!(slice.len(), 2, "Slice length not right");
//         assert_eq!(slice[0].index, 6, "Element index doesn't match");
//         assert_eq!(slice[1].index, 7, "Element index doesn't match");

//         // Go back to the previous slice
//         slice = buffer.cpgs_in_range("chr1".as_bytes().as_ref(), 0, 50).expect("Empty slice");
//         assert_eq!(slice.len(), 8, "Slice length not right");
//         assert_eq!(slice[7].index, 7, "Element index doesn't match");

//         // End slice on position
//         slice = buffer.cpgs_in_range("chr1".as_bytes().as_ref(), 7, 100).expect("Empty slice");
//         assert_eq!(slice.len(), 5, "Slice length not right");
//         assert_eq!(slice.first().expect("empty slice").index, 3, "First element index doesn't match");
//         assert_eq!(slice.last().expect("empty slice").index, 7, "Last element index doesn't match");

//         // Slice beyond max value
//         slice = buffer.cpgs_in_range("chr1".as_bytes().as_ref(), 50, 1000).expect("Empty slice");
//         assert_eq!(slice.len(), 4, "Slice length not right");
//         assert_eq!(slice.first().expect("empty slice").index, 8, "First element index doesn't match");
//         assert_eq!(slice.last().expect("empty slice").index, 11, "Last element index doesn't match");

//         // Slice in second chrom
//         slice = buffer.cpgs_in_range("chr2".as_bytes().as_ref(), 9, 13).expect("Empty slice");
//         assert_eq!(slice.len(), 2, "Slice length not right");
//         assert_eq!(slice.first().expect("empty slice").index, 0, "First element index doesn't match");
//         assert_eq!(slice.last().expect("empty slice").index, 1, "Last element index doesn't match");
//         Ok(())
//     }

//     #[test]
//     fn can_clean_buffer() -> Result<()>
//     {
//         let mut buffer = CpgBuffer::with_file(Cursor::new(COORD_BED))?;
//         buffer.load_data_until("chr1".as_bytes(), 50).expect("Failed to read slice");
//         assert!(buffer.cpg_buffer.back().expect("Empty buffer").index>=7,"Buffer too short, must be missing data");

//         buffer.clear_buffer_until("chr1".as_bytes(), 10);
//         let first_elem = buffer.cpg_buffer.front().expect("Buffer empty");
//         assert_eq!(first_elem.index, 4);

//         buffer.skip_to_contig("chr2".as_bytes());
//         buffer.load_data_until("chr2".as_bytes(), 18).expect("Failed to read slice");
//         assert_eq!(buffer.cpg_buffer.back().expect("Empty buffer").contig, "chr2".as_bytes().as_ref(),"Buffer too short, must be missing data");
//         buffer.clear_buffer_until("chr2".as_bytes(), 15);
//         assert_eq!(buffer.cpg_buffer.front().expect("Buffer empty").position, 16, "Cleaned to wrong element");
//         Ok(())
//     }
// }