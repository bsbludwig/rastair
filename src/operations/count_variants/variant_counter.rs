use rust_htslib::bam::pileup::Alignment;
use rust_htslib::bam::{IndexedReader, Read};

use std::collections::HashMap;
use std::path::Path;
use std::fs;
use log::{debug, warn, error};
use anyhow::Result;

use crate::sequence_segment::{SequenceSegmentIterator, SequenceSegment};

// Faster hashing than built-in algo
use hashers::fx_hash::FxHasher;
use std::hash::BuildHasherDefault;

use super::{MAX_DEPTH, FLAGS, ReadMaskSetting, ReadMask, VariantCount, NucleotideCount};

/// Configuration of a variant counter
pub struct VariantCounterConfig<P>
{
    /// Path to alignment file
    pub bam_path: P,
    /// Min mapping quality
    pub min_mapq: u8,
    /// Min base quality
    pub min_baseq: u8,
    /// Max depth per alignment position
    pub max_depth: u32,
    /// Only reads that match all these flags will be considered
    pub required_flags: u16,
    /// Any read that matches these flags will be ignored
    pub excluded_flags: u16,
    /// Do not reconcile overlapping read pairs
    pub keep_overlaps: bool,
    /// Mask out a certain number of bases from the left and right for OT reads
    pub ot_mask: ReadMaskSetting,
    /// Mask out a certain number of bases from the left and right for OB reads
    pub ob_mask: ReadMaskSetting,
    /// set the number of threads to use in htslib internally
    pub htslib_threads: usize
}

impl <P: AsRef<Path> + std::fmt::Debug> VariantCounterConfig<P>
{
    pub fn with_path(bam_path: P) -> Result<Self>
    {
        let v = VariantCounterConfig
        {
            bam_path,
            min_mapq: 1,
            min_baseq: 10,
            max_depth: MAX_DEPTH,
            required_flags: FLAGS.is_paired | FLAGS.is_properly_paired,
            excluded_flags: FLAGS.is_failed | FLAGS.is_not_primary | FLAGS.is_unmapped | FLAGS.mate_is_unmapped | FLAGS.is_duplicate | FLAGS.is_supplemental,
            keep_overlaps: false,
            ot_mask: ReadMaskSetting { r1: ReadMask(0, 0), r2: ReadMask(0, 0) },
            ob_mask: ReadMaskSetting { r1: ReadMask(0, 0), r2: ReadMask(0, 0) },
            htslib_threads: 0
        };
        Ok(v)
    }
}

/// Count variants (ie modifications) in sequence chunks
pub struct VariantCounter<P: AsRef<Path> + std::fmt::Debug>
{
    config: VariantCounterConfig<P>,
    bam:	IndexedReader,
    bam_index: Vec<(Vec<u8>, u64, u64)>
}

impl <P: AsRef<Path> + std::fmt::Debug> VariantCounter<P>
{
    /// Initiate a new reader from a configuration object
    pub fn with_config(config: VariantCounterConfig<P>) -> Result<Self>
    {
        let mut bam = IndexedReader::from_path(&config.bam_path)?;
        if config.htslib_threads > 0
        {
            bam.set_threads(config.htslib_threads)?;
        }
        
        // intersect the index files
        let bam_index: Vec<(Vec<u8>, u64, u64)> = bam
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
        Ok(VariantCounter {
            config,
            bam,
            bam_index
        })
    }

    pub fn count_from_file(&mut self, fasta_path: P) -> Result<VariantCounterIterator<P>>
    {
        let iterator = VariantCounterIterator::with_file_and_counter(fasta_path, self)?;
        Ok(iterator)
    }

    pub fn count_from_file_with_step_size(&mut self, fasta_path: P, step_size: usize) -> Result<VariantCounterIterator<P>>
    {
        let iterator = VariantCounterIterator::with_file_and_counter_and_size(fasta_path, self, step_size)?;
        Ok(iterator)
    }
    
    fn index(&self) -> &Vec<(Vec<u8>, u64, u64)>
    {
        &self.bam_index
    }

    /// Generate a closure that can be used to filter alignments, given the configuration settings
    /// This is implemented as a class and not a member function to avoid mutable/immutable ref issues
    fn generate_alignemnt_filter<'a>(config: &'a VariantCounterConfig<P>) -> impl Fn(&Alignment<'a>) -> bool
    {
        let filter_closure = |alignment: &Alignment| -> bool
        {
            if alignment.is_del() || alignment.is_refskip() {
                return false;
            }
            let qpos = alignment.qpos().unwrap(); // safe cause we checked for deletions before
            let record = alignment.record();
            let seq_len = record.seq_len();
            
            if seq_len == 0
            {
                return false;
            }

            let mut filter = record.mapq() >= config.min_mapq;
            // Require all "required" flags (properly paired/aligned)
            filter &= (record.flags() & config.required_flags) == config.required_flags;
            // exclude reads matching _any_ excluded flag
            filter &= (record.flags() & config.excluded_flags) == 0;
            if filter == false
            {
                return false;
            }
            
            let qual = record.qual()[qpos];
            
            println!("qual: {}", qual);
            if qual < config.min_baseq
            {
                return false;
            }

            // first in pair
            if record.flags() & FLAGS.is_first_in_pair > 0
            {
                // F1R2
                if record.flags() & FLAGS.mate_is_reverse_strand > 0 {
                    // Ensure that there's at least one base left after soft-trimming
                    if seq_len < config.ot_mask.r1.0+config.ot_mask.r1.1 + 1 {
                        return false;
                    }

                    if qpos < config.ot_mask.r1.0 || qpos > seq_len-config.ot_mask.r1.1-1
                    {
                        return false;
                    }
                }
                else // F2R1
                {
                    if seq_len < config.ob_mask.r1.0+config.ob_mask.r1.1 + 1 {
                        return false;
                    }

                    // I'm flipping the start/end here, because the R1 of the OB is reversed but 
                    // samtools reports it in ref direction, so if I want to remove 5 bases from the start
                    // of the read, that's actually the "end" in the coordinate system that htslib provides
                    if qpos < config.ob_mask.r1.1 || qpos > seq_len-config.ob_mask.r1.0-1
                    {
                        return false;
                    }
                }
            }
            else 
            {
                // F1R2
                if record.flags() & FLAGS.is_reverse_strand > 0 {
                    if seq_len < config.ot_mask.r2.0+config.ot_mask.r2.1 + 1 {
                        return false;
                    }
                    // also flipped the end/start mask, cause the read is mapped in reverse
                    if qpos < config.ot_mask.r2.1 || qpos > seq_len-config.ot_mask.r2.0-1
                    {
                        return false;
                    }
                }
                else // F2R1
                {
                    if seq_len < config.ob_mask.r2.0+config.ob_mask.r2.1 + 1 {
                        return false;
                    }
                    if qpos < config.ob_mask.r2.0 || qpos > seq_len-config.ob_mask.r2.1-1
                    {
                        return false;
                    }
                }
            }
            
            true
        };
        filter_closure
    }

    /// For a given sequence segment, return the number of observed nucleotide counts at all covered
    /// positions
    fn count_variants_in_segment(&mut self, segment: SequenceSegment) -> Option<Vec<VariantCount>>
    {
        //TODO this needs changing to make it more generic, ie allow different types of subsets
        // Search the string segment for all CpG positions
        let cpg_positions = segment.find_cpgs().unwrap_or(Vec::new());

        if cpg_positions.len() == 0
        {
            return Some(Vec::new());
        }
        
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

        // Allocate enough space for output
        let mut output: Vec<VariantCount> = Vec::with_capacity(cpg_positions.len());

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

            let filter_closure = VariantCounter::generate_alignemnt_filter(&self.config);

            if pileup_pos == this_position.pos_in_contig()
            {
                debug!("Found a CpG site at {}", this_position);
                cpg_index = cpg_index + 1;

                // Count conversion vs non-conversion
                let mut var_count = VariantCount::new();
                var_count.contig = segment.contig.clone();
                var_count.pos = this_position.pos_in_contig();
                var_count.ref_base = this_position.base();

                'alignment_loop:
                for alignment in pileup
                                                .alignments()
                                                .filter(filter_closure)
                {
                    let pos = alignment.qpos()?;
                    // This copies memory, so don't repeat
                    let record = alignment.record();
                    let seq = record.seq();

                    let base = seq[pos];
                    let qual = record.qual()[pos];

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

pub struct VariantCounterIterator<'a, P: AsRef<Path> + std::fmt::Debug> 
{
    counter: &'a mut VariantCounter<P>,
    fasta: SequenceSegmentIterator<fs::File>,
}

impl <'a, P: AsRef<Path> + std::fmt::Debug> VariantCounterIterator<'a, P>
{
    pub fn with_file_and_counter(fasta_path: P, counter:&'a mut VariantCounter<P>) -> Result<Self>
    {
        let mut fasta = SequenceSegmentIterator::with_file(&fasta_path)?;
        fasta.subset_to_intervals(counter.index())?;
        
        Ok(VariantCounterIterator {
            counter,
            fasta
        })
    }

    pub fn with_file_and_counter_and_size(fasta_path: P, counter: &'a mut VariantCounter<P>, chunk_size: usize) -> Result<Self>
    {
        let mut fasta = SequenceSegmentIterator::with_file_and_stepsize(&fasta_path, chunk_size)?;
        fasta.subset_to_intervals(counter.index())?;
        
        Ok(VariantCounterIterator {
            counter,
            fasta
        })
    }
}
impl <'a, P: AsRef<Path> + std::fmt::Debug> Iterator for VariantCounterIterator<'a, P>
{
    type Item = Vec<VariantCount>;

    fn next(&mut self) -> Option<Self::Item>
    {
        let segment = self.fasta.next()?;

        debug!("Process {}", &segment);

        self.counter.count_variants_in_segment(segment)
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

/*====================================================
 = Unit Tests
====================================================*/
#[cfg(test)]
mod tests {
    
    // For testing
    use super::*;
    use std::str::FromStr;
    use std::path::PathBuf;

    fn create_test_config() -> Result<VariantCounterConfig<PathBuf>>
    {
        let fasta_path = PathBuf::from(r"test_data/test.fasta");
        let bam_path = PathBuf::from(r"test_data/test.bam");
        let new_config = VariantCounterConfig::with_path(bam_path)?;
        Ok(new_config)
    }

    #[test]
    fn can_create_config() -> Result<()>
    {
        let config = create_test_config()?;

        assert_eq!(config.max_depth, MAX_DEPTH);
        assert_eq!(config.required_flags, 3);
        assert_eq!(config.excluded_flags, 3852);
        Ok(())
    }

    #[test]
    fn can_filter_alignments() -> Result<()>
    {
        let config = create_test_config()?;
        let filter = VariantCounter::generate_alignemnt_filter(&config);
        let mut bam = IndexedReader::from_path(&config.bam_path)?;

        bam.fetch(("bacteriophage_lambda_CpG", 0, 100))?;
        if let Ok(pileup) = bam.pileup().nth(0).unwrap() {
            let alignments: Vec<Alignment> = pileup.alignments().collect();
            assert_eq!(alignments.len(), 7);
            assert_eq!(alignments.into_iter().filter(filter).count(), 3);
        };
        
        Ok(())
    }

    #[test]
    fn can_filter_with_masking_1() -> Result<()>
    {
        let mut config = create_test_config()?;
        config.ot_mask = ReadMaskSetting::from_str(&"1,0,0,1").unwrap();
        config.ob_mask = ReadMaskSetting::from_str(&"0,1,1,0").unwrap();
        let filter = VariantCounter::generate_alignemnt_filter(&config);
        let mut bam = IndexedReader::from_path(&config.bam_path)?;

        bam.fetch(("bacteriophage_lambda_CpG", 0, 100))?;
        if let Ok(pileup) = bam.pileup().nth(0).unwrap() {
            let alignments: Vec<Alignment> = pileup.alignments().collect();
            assert_eq!(alignments.len(), 7);
            assert_eq!(alignments.into_iter().filter(filter).count(), 0);
        };
        
        Ok(())
    }

    #[test]
    fn can_filter_with_masking_2() -> Result<()>
    {
        let mut config = create_test_config()?;
        config.ot_mask = ReadMaskSetting::from_str(&"1,0,0,1").unwrap();
        config.ob_mask = ReadMaskSetting::from_str(&"1,0,0,1").unwrap();
        let filter = VariantCounter::generate_alignemnt_filter(&config);
        let mut bam = IndexedReader::from_path(&config.bam_path)?;

        bam.fetch(("bacteriophage_lambda_CpG", 0, 100))?;
        if let Ok(pileup) = bam.pileup().nth(0).unwrap() {
            let alignments: Vec<Alignment> = pileup.alignments().collect();
            assert_eq!(alignments.len(), 7);
            assert_eq!(alignments.into_iter().filter(filter).count(), 1);
        };
        
        Ok(())
    }
}