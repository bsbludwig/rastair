/*
 * Parse each read in a region and return the start coordinate of the read, frag length,
 * flag, number of CpGs, number of modified C/Gs, positions of CpGs, mod status of CpGs
 */


use std::{fmt::{Formatter, Display, Debug}, path::{PathBuf, Path}, fs::File, io::{stdout, Write}};

use log::{debug, info, trace, warn};

use anyhow::{bail, Result};
use bio::bio_types::sequence::SequenceReadPairOrientation::{F1R2, F2R1};
use rust_htslib::bam::{ext::BamRecordExtensions, FetchDefinition, IndexedReader, Read, Record};

use crate::{sequence_segment::{GenomicRegion, SequenceSegment, SequenceSegmentIterator}, utils::{FetchDefinitionExt, IndexedReaderExt}};
use crate::utils::RecordExt;

use super::FLAGS;

 /// Store methylation information for a single read
 pub struct PerReadCount
 {
    /// ID of the sequence
    pub region: GenomicRegion,
    /// Flag of read
    pub flag: u16,
    /// Mapq of read
    pub mapq: u8,
    /// Absolute fragment length (non-directional)
    pub frag_length: u32,
    /// Name of read
    pub read_id: String,
    /// Number of CpGs in a read
    pub cpg_count: u16,
    /// Number of modified CpGs
    pub mod_count: u16,
    /// Positions in read  of modified CpGs
    pub mod_cpgs: Vec<usize>,
    /// Positions in read of unmodified CpGs
    pub unmod_cpgs: Vec<usize>,
 }

 impl Display for PerReadCount
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result
    {
        write!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.region.contig,
            self.region.start,
            self.region.end,
            self.read_id,
            self.mapq,
            if self.flag & 16 == 16 {"-"} else {"+"},
            self.frag_length,
            self.flag,
            self.cpg_count,
            self.mod_count,
            self.mod_cpgs.iter().map(|f| f.to_string()).collect::<Vec<String>>().join(","),
            self.unmod_cpgs.iter().map(|f| f.to_string()).collect::<Vec<String>>().join(","))
    }
}


 /// Iterator over reads in a region, returning an object that
 /// can be easily printed and parsed in e.g. R
 pub struct ReadIterator<'a>
 {
    // A sequence object of that region
    sequence: &'a SequenceSegment,
    // all CpGs in the region
    cpgs: Vec<u64>,
    // lowest index in the cpgs vector to start looking for positions inside the read
    cpg_offset: usize,
    // An indexed bam reader
    bam_reader: &'a mut IndexedReader,
    // config
    config: &'a ReadCounterConfig,
    // Ignore reads if they start beyond this position
    pub right_margin: u64,
 }

 impl <'a> ReadIterator <'a>
 {
    /// Constructor with region and readers
    pub fn with_region_and_reader(segment: &'a SequenceSegment, reader: &'a mut IndexedReader, config: &'a ReadCounterConfig) -> Result<ReadIterator<'a>>
    {
        // Fetch reads in region/set reader the correct subset
        reader.fetch((segment.region.contig.as_bytes(), segment.region.start, segment.region.end))?;
        let all_cpgs: Vec<u64> = segment.find_cpgs()
                                 .unwrap_or_default()
                                 .iter()
                                 .map(|p| p.pos_in_contig())
                                 .collect();
        let tiling_window = config.tiling_window_size;
        let right_margin = if segment.is_last
        {
            segment.region.end
        }
        else if segment.region.end > (tiling_window as u64)
        {
            segment.region.end - (tiling_window as u64)
        }
        else
        {
            0
        };
        if right_margin <= segment.region.start
        {
            bail!("Tiling window larger than segment, no reads will be processed!");
        }
        Ok(Self { sequence: segment, cpgs: all_cpgs, cpg_offset: 0, bam_reader: reader, config, right_margin})
    }

    /// Progress the pointer in the bam file iterator to the next read that matches the criteria set in the config.
    /// Return None if there's no more reads to find in the current segment.
    fn find_next_read(&mut self) -> Option<Record>
    {
        let mut record: Record = Record::new();
        loop {
            match self.bam_reader.read(&mut record)?
            {
                Ok(_)  => (),
                Err(e)  => {
                    info!("Failed to fetch next read: {}", e);
                    return None;
                }
            }

            // if read starts beyond the right margin, terminate
            if (record.pos() as u64) > self.right_margin
            {
                debug!("Skipping read {} that starts beyond right margin: {} <= {}", String::from_utf8(Vec::from(record.qname())).unwrap_or_default(), record.pos(), self.right_margin);
                return None;
            }

            // check if the current read started before the segment start, which means
            // it was processed in a previous iteration. Skip
            if record.pos() > 0 && (record.pos() as u64) <= self.sequence.region.start
            {
                debug!("Skipping read {} that started in a previous segment: {} <= {}", String::from_utf8(Vec::from(record.qname())).unwrap_or_default(), record.pos(), self.sequence.region.start);
                continue;
            }
            // skip reads that lack required flags
            if record.flags() & self.config.required_flags != self.config.required_flags
            {
                debug!("Read {} lacks required flag: {} vs {}", String::from_utf8(Vec::from(record.qname())).unwrap_or_default(), record.flags(), self.config.required_flags);
                continue;
            }
            // skip reads that match an excluded flag
            else if record.flags() & self.config.excluded_flags > 0
            {
                debug!("Read {} has excluded flag: {} vs {}", String::from_utf8(Vec::from(record.qname())).unwrap_or_default(), record.flags(), self.config.excluded_flags);
                continue;
            }

            let read_pair_orientation = record.read_pair_orientation_lenient(self.config.exclude_ambiguous);

            match read_pair_orientation
            {
                F1R2    => (),
                F2R1    => (),
                _       => {
                    debug!("Skipping incorrectly paired read {}. (flag {})", String::from_utf8(Vec::from(record.qname())).unwrap_or_default(), record.flags());
                    continue;
                }
            }

            break;
        };
    Some(record)
    }
 }

 impl Iterator for ReadIterator<'_>
{
    type Item = PerReadCount;

    fn next(&mut self) -> Option<Self::Item>
    {
        let record = self.find_next_read()?;

        /* Check if the read extends beyond the end of the margin
         */
        let alignment = record.cigar();
        if (alignment.end_pos() as u64) >= self.sequence.region.end
        {
            warn!("Read {} extends beyond end of segment {}, tiling window set too short: {}", std::str::from_utf8(record.qname()).unwrap_or_default(), self.sequence, self.config.tiling_window_size);
        }

        let mut next_read_count = PerReadCount {
            region: GenomicRegion{contig: self.sequence.region.contig.clone(), start: record.pos().abs() as u64, end: alignment.end_pos() as u64},
            flag: record.flags(),
            mapq: record.mapq(),
            frag_length: record.insert_size().abs() as u32,
            read_id: String::from_utf8(Vec::from(record.qname())).unwrap_or_default(),
            cpg_count: 0,
            mod_count: 0,
            mod_cpgs: Vec::new(),
            unmod_cpgs: Vec::new()
        };

        if self.cpgs.len() == 0 || self.cpg_offset >= self.cpgs.len()
        {
            debug!("No CpGs in segment or all CpGs already processed, return empty read count");
            return Some(next_read_count);
        }
        // find the first cpg position that could be in this read
        // since reads should arrive here in sorted order, we can
        // be sure that no future read will start to the left of the
        // current position
        loop
        {
            let next_cpg_pos = self.cpgs[self.cpg_offset];
            if next_read_count.region.start <= next_cpg_pos
            {
                if next_read_count.region.end <= next_cpg_pos
                {
                    debug!("Next CpG pos beyond end of read: {} >= {} (#{} in list)", next_cpg_pos, next_read_count.region.end, self.cpg_offset);
                    // we're further than the last element in the alignment
                    // get a new read that might align a bit further right
                    return Some(next_read_count);
                }
                else
                {
                    debug!("Found CpGs within read {} {}:{}-{} ({})", next_read_count.read_id, next_read_count.region.contig, next_read_count.region.start, next_read_count.region.end, next_cpg_pos);
                    // There's something to process
                    break
                }
            }
            self.cpg_offset = self.cpg_offset+1;
            if self.cpg_offset >= self.cpgs.len()
            {
                debug!("No more CpGs in segment {} after {}", self.sequence, next_read_count.read_id);
                return Some(next_read_count);
            }
        }
        let mut offset=0;
        // Find all CpGs in the subset covered by the read
        'block_loop:
        for block in record.aligned_pairs()
        {
            trace!("Processing block {}/{}", block[0], block[1]);
            if (self.cpg_offset + offset) >= self.cpgs.len()
            {
                debug!("Found all CpGs in segment, stop");
                break;
            }

            loop {
                let pos = self.cpgs[self.cpg_offset + offset];
                if pos < block[1] as u64
                {
                    debug!("CpG pos behind current read pos, need to catch up");
                    offset = offset+1;
                    if self.cpg_offset + offset >= self.cpgs.len()
                    {
                        break 'block_loop;
                    }
                }
                else if pos >= next_read_count.region.end
                {
                    debug!("Next CpG ({}) is outside this read ({}-{}), end loop", pos, next_read_count.region.start, next_read_count.region.end);
                    break 'block_loop;
                }
                else
                {
                    break;
                }
            }
            let cpg_pos = self.cpgs[self.cpg_offset + offset];
            if cpg_pos > block[1] as u64
            {
                continue;
            }

            let base = self.sequence.sequence[(cpg_pos - self.sequence.region.start) as usize];
            let read_base = record.seq()[block[0] as usize].to_ascii_uppercase();

            // check if the position in this read is meaningful
            debug!("Processing pos {} in read {}: {} vs {}, flag {}", cpg_pos, next_read_count.read_id, char::from_u32(base as u32).unwrap_or_default(), char::from_u32(read_base as u32).unwrap_or_default(), record.flags());
            let read_pair_orientation = record.read_pair_orientation_lenient(self.config.exclude_ambiguous);
            match read_pair_orientation
            {
                F1R2    =>
                {
                    match base
                    {
                        b'C' =>
                        {
                            next_read_count.cpg_count = next_read_count.cpg_count + 1;
                            if read_base == b'T'
                            {
                                // mod
                                next_read_count.mod_count = next_read_count.mod_count+1;
                                next_read_count.mod_cpgs.push((block[1].abs_diff(next_read_count.region.start as i64)) as usize);
                            }
                            else if read_base == b'C'
                            {
                                //unmod
                                next_read_count.unmod_cpgs.push((block[1].abs_diff(next_read_count.region.start as i64)) as usize);
                            }
                            else
                            {
                                debug!("SNP or sequencing error");
                            }
                        },
                        _       => continue
                    }
                },
                F2R1    =>
                {
                    match base
                    {
                        b'G' =>
                        {
                            next_read_count.cpg_count = next_read_count.cpg_count + 1;
                            if read_base == b'A'
                            {
                                // mod
                                next_read_count.mod_count = next_read_count.mod_count+1;
                                next_read_count.mod_cpgs.push(block[0] as usize);
                            }
                            else if read_base == b'G'
                            {
                                //unmod
                                next_read_count.unmod_cpgs.push(block[0] as usize);
                            }
                            else
                            {
                                debug!("SNP or sequencing error");
                            }
                        },
                        _       => continue
                    }
                },
                _       => continue
            };
        }
        Some(next_read_count)
    }
}

pub struct ReadCounterConfig
{
    pub bam_path: PathBuf,
    min_mapq: u8,
    required_flags: u16,
    excluded_flags: u16,
    htslib_threads: u8,
    all_reads: bool,
    exclude_ambiguous: bool,
    tiling_window_size: usize
}

impl ReadCounterConfig
{
    pub fn with_path(bam_path: impl AsRef<Path> + Debug) -> Result<Self>
    {
        let v = ReadCounterConfig
        {
            bam_path: bam_path.as_ref().to_owned(),
            min_mapq: 1,
            required_flags: FLAGS.is_paired | FLAGS.is_properly_paired,
            excluded_flags: FLAGS.is_failed | FLAGS.is_not_primary | FLAGS.is_unmapped | FLAGS.mate_is_unmapped | FLAGS.is_duplicate | FLAGS.is_supplemental,
            htslib_threads: 0,
            all_reads: false,
            exclude_ambiguous: false,
            tiling_window_size: 200
        };
        Ok(v)
    }
}

pub fn run_caller(
    bam_path: &PathBuf,
    fasta_path: &PathBuf,
    region_option: &Option<String>,
    mapq_option: &Option<u8>,
    chunk_size_option: &Option<usize>,
    read_length_option: &Option<usize>,
    req_flags_option: &Option<u16>,
    excl_flags_option: &Option<u16>,
    all_reads_option: &Option<bool>,
    exclude_ambiguous_option: &Option<bool>,
    read_threads_option: &Option<u8>) -> Result<()>
{
    /* Read fasta index, and open fasta file for tokenising */
    debug!("Reading fasta and index from {}", fasta_path.display());

    let mut config = ReadCounterConfig::with_path(bam_path.clone())?;
    if let Some(min_mapq) = *mapq_option {
        config.min_mapq = min_mapq;
    }

    if let Some(flags) = *req_flags_option {
        config.required_flags = flags;
    }
    if let Some(flags) = *excl_flags_option {
        config.excluded_flags = flags;
    }
    if let Some(threads) = *read_threads_option {
        config.htslib_threads = threads;
    }

    if let Some(all_reads) = *all_reads_option {
        config.all_reads = all_reads;
    }

    if let Some(exclude_ambiguous) = *exclude_ambiguous_option {
        config.exclude_ambiguous = exclude_ambiguous;
    }

    // Create a sequence iterator
    let chunk_size = chunk_size_option.unwrap_or_default();
    let mut iterator: SequenceSegmentIterator<File> =
        if chunk_size == 0
        {
            SequenceSegmentIterator::with_file(fasta_path)?
        }
        else
        {
            SequenceSegmentIterator::with_file_and_stepsize(fasta_path, chunk_size)?
        };

    if let Some(max_read_length) = read_length_option
    {
        config.tiling_window_size = *max_read_length;
    }

    if config.tiling_window_size == 0
    {
        warn!("Setting tiling window to 0 will cause missing data!");
    }
    // We only actually process reads that start up until segment_end-tiling_window.
    // The overhang is needed to get the CpG loci info for reads that extend beyond the margin
    iterator.set_tiling(config.tiling_window_size)?;

    // Create a bam reader
    // TODO this should probably all happen in some constructor somewhere
    let mut bam = IndexedReader::from_path(&config.bam_path)?;
    if config.htslib_threads > 0
    {
        bam.set_threads(config.htslib_threads as usize)?;
    }
    if let Some(region) = region_option
    {
        let new_region =
        // Need to de-parse region and extend by the read-length to include reads that extend beyond the end
        match FetchDefinition::from_region_string(region)?
        {
            FetchDefinition::RegionString(chr, start, end) =>
            {
                format!("{}:{}-{}", std::str::from_utf8(chr).unwrap_or_default(), start, end+(config.tiling_window_size as i64)-1)
            },
            FetchDefinition::String(_) | FetchDefinition::All =>
            {
                region.to_owned()
            },
            _   => bail!("Failed to parse region string: {}", region)
        };
        iterator.subset_to_region(&new_region)?;
    }

    let bam_index = bam.expanded_index()?;
    // only process regions that are actually in the bam file
    iterator.subset_to_intervals(&bam_index)?;

    let mut lock = stdout().lock();
    writeln!(lock, "#chr\tstart\tend\tread_id\tmapq\torientation\tinsert_size\tflag\tnum_cpg\tnum_mod\tmod_cps\tunmod_cpgs")?;

    while let Some(segment) = iterator.next()
    {
        info!("Process reads in {}", segment);
        let mut read_iterator = ReadIterator::with_region_and_reader(&segment, &mut bam, &config)?;
        while let Some(read_info) = read_iterator.next()
        {
            // Skip empty reads unless explicitly requested
            if !config.all_reads && read_info.cpg_count == 0 {
                continue;
            }
            writeln!(lock, "{}", read_info)?;
        }
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
    use anyhow::Result;

    const BAMFILE: &'static str = r"test_data/test.bam";
    const FASTAFILE: &'static str = r"test_data/test.fasta";

    fn create_config() -> Result<ReadCounterConfig>
    {
        let bam_path = PathBuf::from(BAMFILE);
        let rcc = ReadCounterConfig::with_path(bam_path)?;
        Ok(rcc)
    }

    #[test]
    fn can_create_config() -> Result<()>
    {
        let rcc = create_config()?;
        assert_eq!(rcc.bam_path.to_str().unwrap(), BAMFILE);
        Ok(())
    }

    #[test]
    fn can_create_counter() -> Result<()>
    {
        let config = create_config()?;
        let mut bam = IndexedReader::from_path(&config.bam_path)?;
        bam.set_threads(1)?;
        let bam_index = bam.expanded_index()?;

        let mut iterator: SequenceSegmentIterator<File> = SequenceSegmentIterator::with_file(FASTAFILE)?;
        iterator.subset_to_region(&"bacteriophage_lambda_CpG".to_string())?;
        iterator.subset_to_intervals(&bam_index)?;

        if let Some(segment) = iterator.next()
        {
            let counter = ReadIterator::with_region_and_reader(&segment, &mut bam, &config)?;
            assert_eq!(counter.right_margin, segment.region.end-(config.tiling_window_size as u64));
        }
        else
        {
            assert!(false);
        }
        Ok(())
    }
}