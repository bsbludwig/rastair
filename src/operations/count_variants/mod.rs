use rust_htslib::bam::pileup::Alignment;
use rust_htslib::bam::{IndexedReader, Read};

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::{fmt, fs};
use std::error::Error;
use std::io::{stdout, Write};
use log::{debug, warn, error};
use anyhow::{Result, bail};

use crate::sequence_segment::SequenceSegmentIterator;

// Faster hashing than built-in algo
use hashers::fx_hash::FxHasher;
use std::hash::BuildHasherDefault;

use super::MAX_DEPTH;
use super::FLAGS;

pub struct NucleotideCount
{
    pub a: i32,
    pub c: i32,
    pub g: i32,
    pub t: i32,
    pub n: i32,
}

impl fmt::Display for NucleotideCount
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result
    {
        write!(f, "A: {} C: {} G: {} T: {} N: {}", self.a, self.c, self.g, self.t, self.n)
    }
}
pub struct VariantCount
{
    pub contig: String,
    pub pos: u64,
    pub ref_base: u8,
    pub top: NucleotideCount,
    pub bottom: NucleotideCount,
}

impl VariantCount
{
    pub fn new() -> Self
    {
        let fw = NucleotideCount{a: 0, c: 0, g: 0, t: 0, n: 0};
        let rv = NucleotideCount{a: 0, c: 0, g: 0, t: 0, n: 0};
        VariantCount { contig: String::new(), pos: 0, ref_base: 0, top: fw, bottom: rv }
    }
}

impl fmt::Display for VariantCount
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> fmt::Result
    {
        let char = char::from_u32(self.ref_base as u32).unwrap_or_default();
        write!(f, "{}:{} ({})\nFW\t{}\nRV\t{}", self.contig, self.pos, char, self.top, self.bottom)
    }
}

pub struct VariantCounterConfig<P>
{
    pub fasta_path: P,
    pub bam_path: P,
    pub min_mapq: u8,
    pub min_baseq: u8,
    pub max_depth: u32,
    pub chunk_size: usize,
    pub required_flags: u16,
    pub excluded_flags: u16,
    pub keep_overlaps: bool
}

impl <P: AsRef<Path> + std::fmt::Debug> VariantCounterConfig<P>
{
    pub fn with_paths(fasta_path: P, bam_path: P) -> Result<Self>
    {
        let v = VariantCounterConfig
        {
            fasta_path,
            bam_path,
            min_mapq: 1,
            min_baseq: 10,
            max_depth: MAX_DEPTH,
            chunk_size: 100000,
            required_flags: FLAGS.is_paired | FLAGS.is_properly_paired,
            excluded_flags: FLAGS.is_failed | FLAGS.is_not_primary | FLAGS.is_unmapped | FLAGS.mate_is_unmapped | FLAGS.is_duplicate | FLAGS.is_supplemental,
            keep_overlaps: false,
        };
        Ok(v)
    }
}

pub struct VariantCounter<P: AsRef<Path> + std::fmt::Debug>
{
    config: VariantCounterConfig<P>,
    bam:	IndexedReader,
    fasta: SequenceSegmentIterator<fs::File>,
}

impl <P: AsRef<Path> + std::fmt::Debug> VariantCounter<P>
{
    pub fn with_config(config: VariantCounterConfig<P>) -> Result<Self>
    {
        let mut bam = IndexedReader::from_path(&config.bam_path)?;
        let mut fasta = SequenceSegmentIterator::with_file_and_stepsize(&config.fasta_path, config.chunk_size)?;

        // intersect the index files
        let bam_index: Vec<(Vec<u8>, u64, u64)> = bam
                        .index_stats()
                        .unwrap()
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
                                let header = bam.header();
                                let seq_id = header.tid2name(idx.0 as u32);
                                (Vec::from(seq_id), 0, idx.1)
                            }
                        ).collect();
        match fasta.subset_to_intervals(&bam_index)
        {
            Some(_) => Ok(
                VariantCounter
                {
                    config,
                    bam,
                    fasta
                }
            ),
            None    => bail!("No sequences intersect between fasta and bam")
        }
    }
}

impl <P: AsRef<Path> + std::fmt::Debug> Iterator for VariantCounter<P>
{
    type Item = Vec<VariantCount>;

    fn next(&mut self) -> Option<Self::Item>
    {
        let segment = self.fasta.next()?;

        debug!("Process {}", &segment);
        //TODO this needs changing to make it more generic, ie allow different types of subsets
        // Search the string segment for all CpG positions
        let cpg_positions = segment.find_cpgs().unwrap_or(Vec::new());

        if cpg_positions.len() == 0
        {
            return Some(Vec::new());
        }

        // Allocate enough space for output
        let mut output: Vec<VariantCount> = Vec::with_capacity(cpg_positions.len());

        /* Fetch the pileup for the region from the bam file, and go
        * through all CpG positions, performing whatever calculation
        * needs to be performed. Stream the output to a writer that
        * writes the results to STDOUT or some file.
        */
        match self.bam.fetch((&segment.contig[..], segment.start, segment.stop))
        {
            Ok(_) => (),
            Err(e) =>
            {
                warn!("Error fetching sequence for {}: {}", &segment, e);
                return None;
            }
        }

        let mut pileup_iterator = self.bam.pileup();
        pileup_iterator.set_max_depth(self.config.max_depth);

        // Pre-create a hash to keep; pre-allocate memory to hold as many reads as the max read depth allows
        let mut read_hash:HashMap<Vec<u8>, (u8, u8), BuildHasherDefault<FxHasher>> = HashMap::with_capacity_and_hasher(self.config.max_depth as usize, BuildHasherDefault::<FxHasher>::default());

        // Find the next CpG position that is greater than or equal to the current pos
        // If none exists, we've reached the end
        let mut cpg_index = 0;
        
        'pileup_loop:
        for pileups in pileup_iterator
        {
            let pileup = match pileups
            {
                Ok(p) => p,
                Err(e)  =>
                {
                    error!("Error reading pileup at {} {}", &segment, e);
                    continue;
                }
            };
            let pileup_pos = pileup.pos() as u64;
            let this_position = &cpg_positions[cpg_index];

            if pileup_pos == this_position.pos_in_contig()
            {
                debug!("Found a CpG site at {}", this_position);
                cpg_index = cpg_index + 1;

                // Count conversion vs non-conversion
                let mut var_count = VariantCount::new();
                var_count.contig = segment.contig.clone();
                var_count.pos = this_position.pos_in_contig();
                var_count.ref_base = this_position.base();

                let filter_closure = |alignment: &Alignment| -> bool
                {
                    let record = alignment.record();
                    let mut filter = !(alignment.is_del() || alignment.is_refskip());
                    filter &= record.mapq() >= self.config.min_mapq;
                    // Require all "required" flags (properly paired/aligned)
                    filter &= (record.flags() & self.config.required_flags) == self.config.required_flags;
                    // exclude reads matching _any_ excluded flag
                    filter &= (record.flags() & self.config.excluded_flags) == 0;
                    filter
                };

                'alignment_loop:
                for alignment in
                                                pileup
                                                .alignments()
                                                .filter(filter_closure)
                {
                    let pos = alignment.qpos()?;
                    // This copies memory, so don't repeat
                    let record = alignment.record();
                    let seq = record.seq();

                    if seq.len() == 0
                    {
                        continue 'alignment_loop;
                    }

                    let qual = record.qual()[pos];
                    if qual < self.config.min_baseq
                    {
                        continue 'alignment_loop;
                    }

                    let base = seq[pos];

                    let nuc_counts =
                    {
                        if (record.flags() & (FLAGS.is_first_in_pair | FLAGS.mate_is_reverse_strand) == 0)
                        || (record.flags() & (FLAGS.is_second_in_pair | FLAGS.is_reverse_strand) == 0)
                        {
                            &mut var_count.top
                        }
                        else
                        {
                            &mut var_count.bottom
                        }
                    };

                    // Check for overlapping reads if requested
                    if !self.config.keep_overlaps
                    {
                        if let Some(pos_tuple) = read_hash.get(record.qname())
                        {
                            debug!("Found overlapping read pair {} at pos {} with bases {} vs {}", std::str::from_utf8(record.qname()).unwrap_or_default(), pos, char::from_u32(base as u32).unwrap_or_default(), char::from_u32(pos_tuple.0 as u32).unwrap_or_default());
                            //debug!("Found an overlapping read at pos {} in fragment {} ({} vs {})", this_position, String::from_utf8(Vec::from(record.qname())).unwrap_or_default(), pos, pos_tuple.0);
                            if pos_tuple.0 == base
                            {
                                // same sequence in each pair, do not double-count but keep previous
                                continue 'alignment_loop;
                            }
                            else 
                            {
                                // Mismatch! Remove previously counted base and ignore this one
                                // TODO: decide whether to follow the example of Methyldackel and
                                // count the higher-quality base - however, unlike Methyldackel
                                // I check for overlaps only _after_ baseq cutoff, so some overlaps
                                // will have already been treated in that way. If there's two high-quality
                                // disagreeing calls, I feel it's better to ignore the whole fragment
                                increment_counter_by(nuc_counts, pos_tuple.0, -1);
                                continue 'alignment_loop;
                            }
                        } 
                        else 
                        {
                            // TODO I'm storing the qual here _in case_ I want to implement quality-based
                            // decision on which read disagreeing base to keep in the future
                            read_hash.insert(Vec::from(record.qname()), (base, qual));
                        }
                    }

                    match increment_counter_by(nuc_counts, base, 1)
                    {
                        Some(()) => (),
                        None => 
                        {
                            let char = char::from_u32(base as u32).unwrap_or_default();
                            warn!("Encountered unknown char {} at {}", char, this_position);
                        }
                    }
                }
                // Empty read hash
                read_hash.clear();
                debug!("{}", var_count);
                output.push(var_count);
                // Check if we're done
                if cpg_index >= cpg_positions.len()
                {
                    break 'pileup_loop;
                }
            }
        }

        Some(output)
    }
}

fn increment_counter_by(nuc_counts: &mut NucleotideCount, base: u8, amount: i32) -> Option<()>
{
    match base
    {
        b'a' => nuc_counts.a += amount,
        b'c' => nuc_counts.c += amount,
        b'g' => nuc_counts.g += amount,
        b't' => nuc_counts.t += amount,
        b'n' => nuc_counts.n += amount,
        b'A' => nuc_counts.a += amount,
        b'C' => nuc_counts.c += amount,
        b'G' => nuc_counts.g += amount,
        b'T' => nuc_counts.t += amount,
        b'N' => nuc_counts.n += amount,
        _   =>
        {
            return None;
        }
    };
    Some(())
}


pub fn run_caller(
    bam_path: &PathBuf,
    fasta_path: &PathBuf,
    mapq_option: &Option<u8>,
    baseq_option: &Option<u8>,
    max_depth_option: &Option<u32>,
    chunk_size_option: &Option<usize>,
    req_flags_option: &Option<u16>,
    excl_flags_option: &Option<u16>) -> Result<(), Box<dyn Error>> 
{
    /* Read fasta index, and open fasta file for tokenising */
    debug!("Reading fasta and index from {}", fasta_path.display());
    
    let mut config = VariantCounterConfig::with_paths(fasta_path, bam_path).unwrap();
    if let Some(min_mapq) = mapq_option {
        config.min_mapq = *min_mapq;
    }
    if let Some(min_baseq) = baseq_option {
        config.min_baseq = *min_baseq;
    }
    if let Some(max_depth) = max_depth_option {
        config.max_depth = *max_depth;
    }
    if let Some(cs) = chunk_size_option {
        config.chunk_size = *cs;
    }
    if let Some(flags) = req_flags_option {
        config.required_flags = *flags;
    }
    if let Some(flags) = excl_flags_option {
        config.excluded_flags = *flags;
    }
    let counter = VariantCounter::with_config(config).unwrap();

    let mut lock = stdout().lock();
    for cpgs in counter
    {
        for cpg in cpgs
        {
            if cpg.ref_base == b'C'
            { // C
                writeln!(lock, "{}\t{}\t{}\t.\t.\t+\t{}\t{}\t{}\t{}\t{}", cpg.contig, cpg.pos, cpg.pos+1, cpg.top.c, cpg.top.t, cpg.bottom.c, cpg.bottom.t, cpg.top.a + cpg.top.c + cpg.top.g + cpg.top.t + cpg.bottom.a + cpg.bottom.c + cpg.bottom.g + cpg.bottom.t).unwrap();
            }
            else
            { // G
                writeln!(lock, "{}\t{}\t{}\t.\t.\t-\t{}\t{}\t{}\t{}\t{}", cpg.contig, cpg.pos, cpg.pos+1, cpg.bottom.g, cpg.bottom.a, cpg.top.g, cpg.top.a, cpg.top.a + cpg.top.c + cpg.top.g + cpg.top.t + cpg.bottom.a + cpg.bottom.c + cpg.bottom.g + cpg.bottom.t).unwrap();
            }
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
}