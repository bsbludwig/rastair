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
    last_parsed_position: u64,
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
            Self { cpg_buffer, fasta_reader: ssi, last_parsed_position: 0, last_in_segment: false }
        )
    }

    pub fn with_reader_and_stepsize(fasta_reader: IndexedReader<R>, step_size: usize) -> Result<Self>
    {
        let cpg_buffer = BTreeMap::new();
        let ssi = SequenceSegmentIterator::with_reader_and_stepsize(fasta_reader, step_size)?;
        Ok(
            Self { cpg_buffer, fasta_reader: ssi, last_parsed_position: 0, last_in_segment: false }
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
        self.last_parsed_position = segment.region.end-1;
        let mut iter: usize = 0;
        for cpg in segment.find_cpgs().unwrap_or_default()
        {
            let cpg_info = CpgInfo::new(Vec::from(cpg.contig()), cpg.pos_in_contig(), last_index + iter);
            self.cpg_buffer.insert(cpg.pos_in_contig(), cpg_info);
            iter += 1;
        }
        Some(())
    }

    pub fn progress_to_contig(&mut self, chr: &[u8]) -> Option<()>
    {
        self.cpg_buffer.clear();
        self.last_in_segment = false;
        self.last_parsed_position = 0;
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
            if self.last_in_segment || self.last_parsed_position >= end
            {
                break;
            }
            self.parse_next_from_file()?;
        }
        let slice = self.cpg_buffer.range(start..end).map(|(_, v)| v);
        Some(slice.collect())
    }
}

/*====================================================
 = Unit Tests
====================================================*/
#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;
    use anyhow::{Ok, Result};
    use bio::io::fasta::Index;

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

    #[test]
    fn can_create_buffer() -> Result<()>
    {
        let index = Index::new(Cursor::new(FAI_FILE))?;
        let reader = IndexedReader::with_index(Cursor::new(FASTA_FILE), index);
        let mut buffer = CpgBuffer::with_reader_and_stepsize(reader, 12)?;

        let mut rows = buffer.cpgs_in_range("id".as_bytes(), 0, 12).expect("Could not fetch");
        let mut next_row = rows.first().expect("No CpGs found");
        assert_eq!(next_row.contig, Vec::from("id".as_bytes()));
        assert_eq!(next_row.position, 2);
        assert_eq!(next_row.index, 0);

        rows = buffer.cpgs_in_range("id".as_bytes(), 24, 36).expect("Could not fetch");
        next_row = rows.first().expect("No CpGs found");

        assert_eq!(next_row.contig, Vec::from("id".as_bytes()));
        assert_eq!(next_row.position, 24);
        assert_eq!(next_row.index, 4);

        buffer.progress_to_contig("id4".as_bytes());
        rows = buffer.cpgs_in_range("id4".as_bytes(), 0, 12).expect("Could not fetch");
        next_row = rows.first().expect("No CpGs found");

        assert_eq!(next_row.contig, Vec::from("id4".as_bytes()));
        assert_eq!(next_row.position, 2);
        assert_eq!(next_row.index, 0);

        Ok(())
    }
}