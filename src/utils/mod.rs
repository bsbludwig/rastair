use log::warn;
use anyhow::{bail, Ok, Result};
use rust_htslib::bam::{FetchDefinition, IndexedReader, Read, Record};
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

pub trait FetchDefinitionExt {
    fn from_region_string(region: &str) -> Result<FetchDefinition>;
}

impl <'a>FetchDefinitionExt for FetchDefinition<'a> {
    fn from_region_string(region: &str) -> Result<FetchDefinition> {
        if region.len() == 0
        {
            bail!("Empty region string");
        }
        let chr_region_split: Vec<_> = region
            .split(':')
            .collect();

        if chr_region_split.len() == 1
        {
            let chr = chr_region_split[0].as_bytes();
            if chr.len() == 0
            {
                bail!("Unable to parse region string: {}", region);
            }
            return Ok(FetchDefinition::from(chr));
        }
        else if chr_region_split.len() == 2
        {
            let chr = chr_region_split[0].as_bytes();
            if chr.len() == 0
            {
                bail!("Unable to parse region string: {}", region);
            }

            let parts: Vec<_> = chr_region_split[1].split('-').filter(|p| !p.is_empty()).collect();
            if parts.len() != 2
            {
                bail!("Unable to parse region string: {}", region);
            }

            let start = parts[0].parse::<i64>().unwrap_or(0);
            let end = parts[1].parse::<i64>().unwrap_or(0);
            if end < start
            {
                bail!("Unable to parse region string: {}", region);
            }
            return Ok(FetchDefinition::from((chr, start, end)));
        }
        else
        {
            bail!("Unable to parse region string: {}", region);
        }
    }
}
/*====================================================
 = Unit Tests
====================================================*/
#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;
    use rust_htslib::bam::{header::HeaderRecord, Header, HeaderView};

    const READ: &'static str = "A00711:92:HMH3WDSXX:2:1218:15058:1892\t99\t2kb_3_Unmodified\t1\t60\t151M\t=\t364\t514\tCACAGATGTCTGCCTGTTCATCCGCGTCCAGCTCGTTGAGTTTCTCCAGAAGCGTTAATGTCTGGCTTCTGATAAAGCGGGCCATGTTAAGGGCGGTTTTTTCCTGTTTGGTCACTGATGCCTCCGTGTAAGGGGGATTTCTGTTCATGGG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFF\tMC:Z:151M\tMD:Z:151\tPG:Z:MarkDuplicates\tRG:Z:2b683de3-cf0a-4e7d-9c04-2aa4912a4a66\tNM:i:0\tAS:i:151\tXS:i:0";
    fn create_read() -> Result<Record>
    {
        let mut _header = Header::new();
        _header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"2kb_3_Unmodified")
                .push_tag(b"LN", &2018),
        );
        let mut header = HeaderView::from_header(&_header);
        let rec = Record::from_sam(&mut header, READ.as_bytes())?;
        Ok(rec)
    }

    #[test]
    fn orientation_easy() -> Result<()>
    {
        let read = create_read()?;
        assert_eq!(read.read_pair_orientation_lenient(true), F1R2);
        Ok(())
    }

    #[test]
    fn orientation_strict() -> Result<()>
    {
        // make the read R2F1, ie the R2 read starts before the F1
        let mut read = create_read()?;
        read.set_pos(10);
        read.set_insert_size((read.seq_len() -10) as i64);
        read.set_mpos(0);

        assert_eq!(read.read_pair_orientation_lenient(true), R2F1);
        assert_eq!(read.read_pair_orientation_lenient(false), F1R2);
        Ok(())
    }

    #[test]
    fn orientation_ambiguous() -> Result<()>
    {
        // make the read ambiguous, ie F1/R2 start at same place
        let mut read = create_read()?;
        read.set_insert_size(read.seq_len() as i64);
        read.set_mpos(0);

        assert_eq!(read.read_pair_orientation_lenient(true), SequenceReadPairOrientation::None);
        assert_eq!(read.read_pair_orientation_lenient(false), F1R2);
        Ok(())
    }

    #[test]
    fn region_from_string() -> Result<()>
    {
        // make the read ambiguous, ie F1/R2 start at same place
        match FetchDefinition::from_region_string(&"chr1:1-22")?
        {
            FetchDefinition::RegionString(a, b, c) => {
                assert_eq!(std::str::from_utf8(a).unwrap_or_default(), "chr1");
                assert_eq!(b, 1);
                assert_eq!(c, 22);
            },
            _   => assert!(false)
        };

        match FetchDefinition::from_region_string(&"chrX")?
        {
            FetchDefinition::String(a) => {
                assert_eq!(std::str::from_utf8(a).unwrap_or_default(), "chrX");
            },
            _   => assert!(false)
        };

        match FetchDefinition::from_region_string(&"this-is-valid")?
        {
            FetchDefinition::String(a) => {
                assert_eq!(std::str::from_utf8(a).unwrap_or_default(), "this-is-valid");
            },
            _   => assert!(false)
        };
        Ok(())
    }

    #[test]
    fn fail_on_incorrect_string() -> Result<()>
    {
        // make the read ambiguous, ie F1/R2 start at same place
        match FetchDefinition::from_region_string(&"chr1:0")
        {
            Err(_) => assert!(true),
            _      => assert!(false)
        };

        match FetchDefinition::from_region_string(&"chr1:10-6")
        {
            Err(_) => assert!(true),
            _      => assert!(false)
        };

        match FetchDefinition::from_region_string(&"chr1:")
        {
            Err(_) => assert!(true),
            _      => assert!(false)
        };

        match FetchDefinition::from_region_string(&":2-10")
        {
            Err(_) => assert!(true),
            _      => assert!(false)
        };
        Ok(())
    }
}