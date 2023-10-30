/*
 * Parse each read in a region and return the start coordinate of the read, frag length, 
 * flag, number of CpGs, number of modified C/Gs, positions of CpGs, mod status of CpGs
 */


use std::{fmt::{Formatter, Display, Debug}, path::{PathBuf, Path}, error::Error, fs::File, io::{stdout, Write}};
use log::{debug, info, trace};

use anyhow::Result;
use bio::bio_types::sequence::SequenceReadPairOrientation::{F1R2, F2R1};
use rust_htslib::bam::{IndexedReader, Record, Read, ext::BamRecordExtensions, record::{CigarStringView, CigarString}};

use crate::sequence_segment::{SequenceSegment, SequenceSegmentIterator};

use super::FLAGS;

 /// Store methylation information for a single read
 pub struct PerReadCount
 {
    /// ID of the sequence
    pub contig: String,
    /// start position in sequence coordinate space
    pub start: u64,
    // end position of alignment in seq coordinates
    pub stop: u64,
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
            self.contig, 
            self.start, 
            self.stop, 
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
    // pre-allocated record
    next_read: Record,
    // config
    config: &'a ReadCounterConfig,
    /// Start position of the first read encountered which extends beyond the end of the interval
    pub right_margin: u64
 }

 impl <'a> ReadIterator <'a>
 {
    /// Constructor with region and readers
    pub fn with_region_and_reader(segment: &'a SequenceSegment, reader: &'a mut IndexedReader, config: &'a ReadCounterConfig) -> Result<ReadIterator<'a>>
    {
        // Fetch reads in region/set reader the correct subset
        reader.fetch((&segment.contig[..], segment.start, segment.stop))?;
        let all_cpgs: Vec<u64> = segment.find_cpgs()
                                 .unwrap_or_default()
                                 .iter()
                                 .map(|p| p.pos_in_contig())
                                 .collect();
        Ok(Self { sequence: segment, cpgs: all_cpgs, cpg_offset: 0, bam_reader: reader, next_read: Record::new(), config, right_margin: segment.stop})
    }
 }

 impl Iterator for ReadIterator<'_>
{
    type Item = PerReadCount;

    fn next(&mut self) -> Option<Self::Item>
    {
        if self.cpgs.len() == 0
        {
            info!("No CpGs in segment {}", self.sequence);
            return None;
        }
        else if self.cpg_offset >= self.cpgs.len()
        {
            debug!("All CpGs in segment {} were processed", self.sequence);
            return None;
        }

        // Fetch the next suitable read into the pre-allocated struct
        #[allow(unused_assignments)] // this is just a hack to set the alignment within the loop. There's probably a smarter way to program this...
        let mut alignment = CigarStringView::new(CigarString::try_from(Vec::new()).unwrap(), 0);
        loop {
            match self.bam_reader.read(&mut self.next_read)? 
            {
                Ok(_)  => (),
                Err(e)  => {
                    info!("Failed to fetch next read: {}", e);
                    return None;
                }
            }
            match self.next_read.read_pair_orientation()
            {
                F1R2    => (),
                F2R1    => (),
                _       => {
                    debug!("Skipping incorrectly paired read. (flag {})", self.next_read.flags());
                    continue;
                }
            }
            // skip reads that lack required flags
            if self.next_read.flags() & self.config.required_flags != self.config.required_flags
            {
                debug!("Read lacks required flag: {} vs {}", self.next_read.flags(), self.config.required_flags);
                continue;
            }
            // skip reads that match an excluded flag
            else if self.next_read.flags() & self.config.excluded_flags > 0
            {
                debug!("Read has excluded flag: {} vs {}", self.next_read.flags(), self.config.excluded_flags);
                continue;
            }

            // check if the current read started before the segment start, which means
            // it was processed in a previous iteration. Skip
            if (self.next_read.pos() as u64) < self.sequence.start
            {
                debug!("Skipping read that started in a previous segment");
                continue;
            }
            // this is expensive, so try to do only once
            alignment = self.next_read.cigar();
            if (alignment.end_pos() as u64) > self.sequence.stop
            {
                debug!("Found a read that extends beyond the end of the current segment, will skip");
                if self.right_margin == self.sequence.stop
                {
                    self.right_margin = self.next_read.pos() as u64;
                }
                continue;
            }
            break; 
        }

        let mut next_read_count = PerReadCount {
            contig: self.sequence.contig.clone(), // this is an allocation - can this be saved? prob not, as we can't predict the lifetime of PerReadCount object
            start: self.next_read.pos().abs() as u64,
            stop: alignment.end_pos() as u64,
            flag: self.next_read.flags(),
            mapq: self.next_read.mapq(),
            frag_length: self.next_read.insert_size().abs() as u32,
            read_id: String::from_utf8(Vec::from(self.next_read.qname())).unwrap_or_default(),
            cpg_count: 0,
            mod_count: 0,
            mod_cpgs: Vec::new(),
            unmod_cpgs: Vec::new()
        };

        // find the first cpg position that could be in this read
        // since reads should arrive here in sorted order, we can
        // be sure that no future read will start to the left of the
        // current position
        loop 
        {
            let next_cpg_pos = self.cpgs[self.cpg_offset];
            if next_read_count.start <= next_cpg_pos
            {
                if next_read_count.stop <= next_cpg_pos
                {
                    debug!("Next CpG pos beyond end of read: {} >= {} (#{} in list)", next_cpg_pos, next_read_count.stop, self.cpg_offset);
                    // we're further than the last element in the alignment
                    // get a new read that might align a bit further right
                    return Some(next_read_count);
                }
                else 
                {
                    debug!("Found CpGs within read {} {}:{}-{} ({})", next_read_count.read_id, next_read_count.contig, next_read_count.start, next_read_count.stop, next_cpg_pos);
                    // There's something to process
                    break
                }
            }
            self.cpg_offset = self.cpg_offset+1;
            if self.cpg_offset >= self.cpgs.len()
            {
                return None;
            }
        }
        let mut offset=0;
        // Find all CpGs in the subset covered by the read
        'block_loop:
        for block in self.next_read.aligned_pairs()
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
            
            let base = self.sequence.sequence[(cpg_pos - self.sequence.start) as usize];
            let read_base = self.next_read.seq()[block[0] as usize].to_ascii_uppercase();
            // check if the position in this read is meaningful
            debug!("Processing pos {} in read {}: {} vs {}, flag {}", cpg_pos, next_read_count.read_id, char::from_u32(base as u32).unwrap_or_default(), char::from_u32(read_base as u32).unwrap_or_default(), self.next_read.flags());
            match self.next_read.read_pair_orientation()
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
                                next_read_count.mod_cpgs.push(block[0] as usize);
                            }
                            else if read_base == b'C'
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
    htslib_threads: u8
}

impl ReadCounterConfig
{
    pub fn with_path(bam_path: impl AsRef<Path> + Debug) -> Result<Self>
    {
        let v = ReadCounterConfig
        {
            bam_path: bam_path.as_ref().to_owned(),
            min_mapq: 1,
            required_flags: FLAGS.is_paired & FLAGS.is_properly_paired,
            excluded_flags: FLAGS.is_failed | FLAGS.is_not_primary | FLAGS.is_unmapped | FLAGS.mate_is_unmapped | FLAGS.is_duplicate | FLAGS.is_supplemental,
            htslib_threads: 0
        };
        Ok(v)
    }
}

pub fn run_caller(
    bam_path: &PathBuf,
    fasta_path: &PathBuf,
    mapq_option: &Option<u8>,
    chunk_size_option: &Option<u32>,
    req_flags_option: &Option<u16>,
    excl_flags_option: &Option<u16>,
    read_threads_option: &Option<u8>) -> Result<(), Box<dyn Error>> 
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

    // Create a sequence iterator
    let chunk_size = chunk_size_option.unwrap_or_default();
    let mut iterator: SequenceSegmentIterator<File> = 
        if chunk_size == 0 
        { 
            SequenceSegmentIterator::with_file(fasta_path)?
        } 
        else 
        { 
            SequenceSegmentIterator::with_file_and_stepsize(fasta_path, chunk_size as usize)?
        };
    
    // Create a bam reader
    // TODO this should probably all happen in some constructor somewhere
    let mut bam = IndexedReader::from_path(&config.bam_path)?;
    if config.htslib_threads > 0
    {
        bam.set_threads(config.htslib_threads as usize)?;
    }
    let bam_index: Vec<(Vec<u8>, u64, u64)> = 
        bam
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
                let header = bam.header();
                let seq_id = header.tid2name(idx.0 as u32);
                (Vec::from(seq_id), 0, idx.1)
            }
        ).collect();
    iterator.subset_to_intervals(&bam_index)?;
    let mut lock = stdout().lock();
    writeln!(lock, "#chr\tstart\tend\tread_id\tmapq\torientation\tinsert_size\tflag\tnum_cpg\tnum_mod\tmod_cps\tunmod_cpgs")?;
    while let Some(segment) = iterator.next()
    {
        let mut read_iterator = ReadIterator::with_region_and_reader(&segment, &mut bam, &config)?;
        while let Some(read_info) = read_iterator.next()
        {
            writeln!(lock, "{}", read_info)?;
        }
        if read_iterator.right_margin < segment.stop
        {
            info!("Will change iterator tiling to cover reads overlapping the boundary");
            iterator.set_tiling((segment.stop - read_iterator.right_margin) as usize)?;
        }
    }
    Ok(())
}