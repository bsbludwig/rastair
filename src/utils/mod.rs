use log::warn;
use anyhow::{Ok, Result};
use rust_htslib::bam::{IndexedReader, Read, Record};
use bio::bio_types::sequence::SequenceReadPairOrientation::{F1R2, F2R1, R1F2, R2F1, self};

/// Extension of a bam record (ie read)
pub trait RecordExt {
    /// optionally allow ambiguous reads, ie R2F1 reads etc
    fn read_pair_orientation_lenient(&self, exclude_ambiguous: bool) -> SequenceReadPairOrientation;
}

impl RecordExt for Record {
    fn read_pair_orientation_lenient(&self, exclude_ambiguous: bool) -> SequenceReadPairOrientation
    {
        let mut read_pair_orientation = self.read_pair_orientation();
        if ! exclude_ambiguous
        {
            read_pair_orientation = match read_pair_orientation
            {
                F1R2 | R2F1 => F1R2,
                F2R1 | R1F2 => F2R1,
                SequenceReadPairOrientation::None => {
                    warn!("Orientation of {} cannot be unambiguously determined", String::from_utf8(Vec::from(self.qname())).unwrap_or_default());

                    if self.is_first_in_template() && self.is_mate_reverse() ||
                    self.is_last_in_template() && self.is_reverse()
                    {
                        F1R2
                    }
                    // F2R1
                    else if self.is_first_in_template() && self.is_reverse() ||
                            self.is_last_in_template() && self.is_mate_reverse()
                    {
                        F2R1
                    }
                    else {
                        SequenceReadPairOrientation::None
                    }
                },
                _   =>  SequenceReadPairOrientation::None // This should be impossible?
            };
        }
        read_pair_orientation
    }
}

///Extensions on IndexedReader
pub trait IndexedReaderExt {
    /// Return the index as a list of chrName, start, end. Excludes contigs with no mapped reads
    fn expanded_index(&mut self) -> Result<Vec<(Vec<u8>, u64, u64)>>;
}

impl IndexedReaderExt for IndexedReader {
    fn expanded_index(&mut self) -> Result<Vec<(Vec<u8>, u64, u64)>>
    {
        let bam_index: Vec<(Vec<u8>, u64, u64)> =
        self
        .index_stats()
        .unwrap_or(Vec::new())
        .iter()
        .filter(
            |idx| -> bool
            {
                // Exclude contigs with no mapped reads
                idx.2 > 0
            }
        )
        .map(
            |idx| -> (Vec<u8>, u64, u64)
            {
                // This is just a lookup, so it's fine to do in a loop
                let header = self.header();
                let seq_id = header.tid2name(idx.0 as u32);
                (Vec::from(seq_id), 0, idx.1)
            }
        ).collect();
        Ok(bam_index)
    }
}