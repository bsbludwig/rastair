use bio::io::bed::Reader;
use bio::bio_types::strand::Strand;

use anyhow::Result;

use std::{collections::VecDeque, fs::File, io::{Read, Seek}, path::Path};

const INITIAL_BUFFER_SIZE: usize = 1000; // how many CpGs covered by the average read pair?

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
    cpg_buffer: VecDeque<CpgInfo>,
    coord_reader: Reader<R>
}

impl <R> CpgBuffer<R>
where R: Read+Seek
{
    #[allow(dead_code)] // Used for testing
    pub fn with_reader(coord_reader: Reader<R>) -> Result<Self>
    {
        let cpg_buffer: VecDeque<CpgInfo> = VecDeque::with_capacity(INITIAL_BUFFER_SIZE);
        Ok(
            Self { cpg_buffer: cpg_buffer, coord_reader: coord_reader }
        )
    }
		pub fn with_file(coord_bed: R) -> Result<Self>
    {
        let coord_reader = Reader::new(coord_bed);
        Ok(
            Self::with_reader(coord_reader)?
        )
    }
    fn parse_next_from_file(&mut self) -> Option<CpgInfo>
    {
        if let Ok(record) = self.coord_reader.records().next()?
        {
            let (chr, start, cpg_index_str) = (record.chrom(), record.start(), record.name()?);
            let cpg_index: usize = cpg_index_str.parse().ok()?;
            Some(CpgInfo::new(Vec::from(chr.as_bytes()), start, cpg_index))
        }
        else
        {
            None
        }
    }

    /// Remove CpG info that is before the start of the last read
    pub fn clear_buffer_until(&mut self, chr: &[u8], end: u64) -> Option<()>
    {
        loop
        {
            // Will return None and thus break the loop if buffer is empty
            let cpg = self.cpg_buffer.front()?;
            if cpg.contig == chr && cpg.position < end || cpg.contig != chr
            {
                self.cpg_buffer.pop_front();
            }
            else
            {
                break;
            }
        }
        Some(())
    }

    pub fn skip_to_contig(&mut self, chr: &[u8]) -> Option<()>
    {
        if self.cpg_buffer.len() > 0
        {
            if let Some(last_cpg) = self.cpg_buffer.back()
            {
                if last_cpg.contig == chr
                {
                    // Already there
                    return Some(());
                }
            }
            else
            {
                // clear buffer, we're moving on
                self.cpg_buffer.clear();
            }
        }
        loop {
            let next_cpg = self.parse_next_from_file()?;
            if next_cpg.contig == chr
            {
                self.cpg_buffer.push_back(next_cpg);
                return Some(());
            }
        }
    }
    // Read ahead until at least the end of the current read
    pub fn load_data_until(&mut self, chr: &[u8], end: u64) -> Option<()>
    {
        if let Some(last_cpg) = self.cpg_buffer.back()
        {
            // already loaded to the end of this chromosome
            if last_cpg.contig != chr
            {
                return None;
                // let old_chr = last_cpg.contig.clone();
                // self.load_data_until(&old_chr, std::u64::MAX);
                // self.clear_buffer_until(&old_chr, std::u64::MAX);
            }
        }

        loop
        {
            let next_row = match self.cpg_buffer.back()
            {
                Some(row)   => row,
                None    => {
                    // Special case for the first entry to be read, where the buffer is empty,
                    // I load one row here and push it in the buffer, then return it as the
                    // first row
                    let next_cpg = self.parse_next_from_file()?;
                    self.cpg_buffer.push_back(next_cpg);
                    self.cpg_buffer.back()?
                }
            };

            if next_row.contig != chr
            {
                // Ran into new chromosome
                break;
            }
            if next_row.position < end
            {
                let next_cpg = self.parse_next_from_file()?;
                self.cpg_buffer.push_back(next_cpg);
            }
            else
            {
                break;
            }
        }
        Some(())
    }

    pub fn cpgs_in_range(&mut self, chr: &[u8], start: u64, end: u64) ->Option<&[CpgInfo]>
    {
        assert!(end > start, "Inverted slice not allowed");
        // bring the internal cache up to the current read end
        // Ignore None ouput, just means we're done parsing the file
        self.load_data_until(chr, end);
        // find the indices of the slice that covers the range
        let (mut low_index,mut high_index): (usize, usize) = (std::usize::MAX, 0);
        for (index, cpg) in self.cpg_buffer.iter().enumerate()
        {
            if chr != cpg.contig
            {
                // we already found some data
                if low_index < self.cpg_buffer.len()
                {
                    high_index = index;
                    break;
                }
                else
                {
                    // skip forward until we find a line where the chrom matches
                    continue;
                }
            }

            if low_index > self.cpg_buffer.len() && cpg.position >= start
            {
                low_index = index;
            }

            if cpg.position >= end
            {
                high_index = index;
                break; // done, found the top end of the slice
            }
        }
        if low_index > self.cpg_buffer.len()
        {
            return None;
        }
        // needed to ensure he slice is continuous
        self.cpg_buffer.make_contiguous();
        if low_index > high_index
        {
            Some(&self.cpg_buffer.as_slices().0[low_index..])
        }
        else
        {
            Some(&self.cpg_buffer.as_slices().0[low_index..high_index])
        }
    }
}

impl CpgBuffer<File>
{
    pub fn with_path<P: AsRef<Path> + std::fmt::Debug>(coord_bed: P) -> Result<Self>
    {
        let cpg_buffer: VecDeque<CpgInfo> = VecDeque::with_capacity(INITIAL_BUFFER_SIZE);
        let coord_reader = Reader::from_file(coord_bed)?;
        Ok(
            Self { cpg_buffer: cpg_buffer, coord_reader: coord_reader }
        )
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

    const COORD_BED: &[u8] = b"chr1\t0\t1\t0\t+
chr1\t1\t2\t1\t-
chr1\t6\t7\t2\t+
chr1\t7\t8\t3\t-
chr1\t10\t11\t4\t+
chr1\t11\t12\t5\t-
chr1\t20\t21\t6\t+
chr1\t21\t22\t7\t-
chr1\t100\t101\t8\t+
chr1\t101\t102\t9\t-
chr1\t110\t111\t10\t+
chr1\t111\t112\t11\t-
chr2\t10\t11\t0\t+
chr2\t11\t12\t1\t-
chr2\t16\t17\t2\t+
chr2\t17\t18\t3\t-
chr2\t18\t19\t4\t+
chr2\t19\t20\t5\t-
";

    #[test]
    fn can_create_buffer() -> Result<()>
    {
        let mut buffer = CpgBuffer::with_file(Cursor::new(COORD_BED))?;
        let mut next_row = buffer.parse_next_from_file().expect("No next row");
        assert_eq!(next_row.contig, Vec::from("chr1".as_bytes()));
        assert_eq!(next_row.position, 0);
        assert_eq!(next_row.index, 0);
        next_row = buffer.parse_next_from_file().expect("No next row");

        assert_eq!(next_row.contig, Vec::from("chr1".as_bytes()));
        assert_eq!(next_row.position, 1);
        assert_eq!(next_row.index, 1);
        Ok(())
    }

    #[test]
    fn can_empty_buffer() -> Result<()>
    {
        let mut buffer = CpgBuffer::with_file(Cursor::new(COORD_BED))?;
        for _i in 0..18
        {
            buffer.parse_next_from_file().expect("No next row");
        }
        assert!(buffer.parse_next_from_file().is_none(), "More data than expected");
        Ok(())
    }

    #[test]
    fn can_load_until() -> Result<()>
    {
        let mut buffer = CpgBuffer::with_file(Cursor::new(COORD_BED))?;
        buffer.load_data_until("chr1".as_bytes(), 50).expect("Failed to read slice");
        assert!(buffer.cpg_buffer.len() >= 8, "Buffer too short, must be missing data");
        Ok(())
    }

    #[test]
    fn can_get_slice() -> Result<()>
    {
        let mut buffer = CpgBuffer::with_file(Cursor::new(COORD_BED))?;
        let mut slice = buffer.cpgs_in_range("chr1".as_bytes().as_ref(), 0, 50).expect("Empty slice");

        assert_eq!(slice.len(), 8, "Slice length not right");
        assert_eq!(slice[7].index, 7, "Element index doesn't match");

        // Shorter slice from a middle start point
        slice = buffer.cpgs_in_range("chr1".as_bytes().as_ref(), 15, 50).expect("Empty slice");
        assert_eq!(slice.len(), 2, "Slice length not right");
        assert_eq!(slice[0].index, 6, "Element index doesn't match");
        assert_eq!(slice[1].index, 7, "Element index doesn't match");

        // Go back to the previous slice
        slice = buffer.cpgs_in_range("chr1".as_bytes().as_ref(), 0, 50).expect("Empty slice");
        assert_eq!(slice.len(), 8, "Slice length not right");
        assert_eq!(slice[7].index, 7, "Element index doesn't match");

        // End slice on position
        slice = buffer.cpgs_in_range("chr1".as_bytes().as_ref(), 7, 100).expect("Empty slice");
        assert_eq!(slice.len(), 5, "Slice length not right");
        assert_eq!(slice.first().expect("empty slice").index, 3, "First element index doesn't match");
        assert_eq!(slice.last().expect("empty slice").index, 7, "Last element index doesn't match");

        // Slice beyond max value
        slice = buffer.cpgs_in_range("chr1".as_bytes().as_ref(), 50, 1000).expect("Empty slice");
        assert_eq!(slice.len(), 4, "Slice length not right");
        assert_eq!(slice.first().expect("empty slice").index, 8, "First element index doesn't match");
        assert_eq!(slice.last().expect("empty slice").index, 11, "Last element index doesn't match");

        // Slice in second chrom
        slice = buffer.cpgs_in_range("chr2".as_bytes().as_ref(), 9, 13).expect("Empty slice");
        assert_eq!(slice.len(), 2, "Slice length not right");
        assert_eq!(slice.first().expect("empty slice").index, 0, "First element index doesn't match");
        assert_eq!(slice.last().expect("empty slice").index, 1, "Last element index doesn't match");
        Ok(())
    }

    #[test]
    fn can_clean_buffer() -> Result<()>
    {
        let mut buffer = CpgBuffer::with_file(Cursor::new(COORD_BED))?;
        buffer.load_data_until("chr1".as_bytes(), 50).expect("Failed to read slice");
        assert!(buffer.cpg_buffer.back().expect("Empty buffer").index>=7,"Buffer too short, must be missing data");

        buffer.clear_buffer_until("chr1".as_bytes(), 10);
        let first_elem = buffer.cpg_buffer.front().expect("Buffer empty");
        assert_eq!(first_elem.index, 4);

        buffer.skip_to_contig("chr2".as_bytes());
        buffer.load_data_until("chr2".as_bytes(), 18).expect("Failed to read slice");
        assert_eq!(buffer.cpg_buffer.back().expect("Empty buffer").contig, "chr2".as_bytes().as_ref(),"Buffer too short, must be missing data");
        buffer.clear_buffer_until("chr2".as_bytes(), 15);
        assert_eq!(buffer.cpg_buffer.front().expect("Buffer empty").position, 16, "Cleaned to wrong element");
        Ok(())
    }
}