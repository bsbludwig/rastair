use rust_htslib::bam::pileup::Alignment;
use rust_htslib::bam::{IndexedReader, Read};
use bio::bio_types::sequence::SequenceReadPairOrientation::{F1R2, F2R1};

use std::collections::HashMap;
use std::path::{Path, PathBuf};
//use std::fs;
use std::fmt::Debug;
use log::{trace, debug, info, warn, error};
use anyhow::Result;

use crate::sequence_segment::SequenceSegment;
use crate::utils::extensions::{IndexedReaderExt, RecordExt};

// Faster hashing than built-in algo
use fxhash::FxBuildHasher;

use super::{ErrorRate, ReadMask, ReadMaskSetting, VariantCount, ERRORRATES, FLAGS, MAX_DEPTH};

#[derive(Clone, Debug)]
/// Configuration of a variant counter
pub struct VariantCounterConfig
{
    /// Path to alignment file
    pub bam_path: PathBuf,
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
    /// Ignore read pairs that start/end at the same position
    pub exclude_ambiguous: bool,
    /// Mask out a certain number of bases from the left and right for OT reads
    pub ot_mask: ReadMaskSetting,
    /// Mask out a certain number of bases from the left and right for OB reads
    pub ob_mask: ReadMaskSetting,
    /// set the number of threads to use in htslib internally
    pub htslib_threads: usize,
    /// Optional region to fetch
    pub region: Option<String>,
    /// Errormodel to use for genotyping
    pub error_model: ErrorRate
}

impl VariantCounterConfig
{
    pub fn with_path(bam_path: impl AsRef<Path> + Debug) -> Result<Self>
    {
        let v = VariantCounterConfig
        {
            bam_path: bam_path.as_ref().to_owned(),
            min_mapq: 1,
            min_baseq: 10,
            max_depth: MAX_DEPTH,
            required_flags: FLAGS.is_paired | FLAGS.is_properly_paired,
            excluded_flags: FLAGS.is_failed | FLAGS.is_not_primary | FLAGS.is_unmapped | FLAGS.mate_is_unmapped | FLAGS.is_duplicate | FLAGS.is_supplemental,
            keep_overlaps: false,
            exclude_ambiguous: false,
            ot_mask: ReadMaskSetting { r1: ReadMask(0, 0), r2: ReadMask(0, 0) },
            ob_mask: ReadMaskSetting { r1: ReadMask(0, 0), r2: ReadMask(0, 0) },
            htslib_threads: 0,
            region: None,
            error_model: ERRORRATES.novaseq_6000
        };
        Ok(v)
    }
}

/// Count variants (ie modifications) in sequence chunks
pub struct VariantCounter
{
    config: VariantCounterConfig,
    bam:	IndexedReader,
    bam_index: Vec<(Vec<u8>, u64, u64)>
}

impl VariantCounter
{
    /// Initiate a new reader from a configuration object
    pub fn with_config(config: VariantCounterConfig) -> Result<Self>
    {
        let mut bam = IndexedReader::from_path(&config.bam_path)?;
        if config.htslib_threads > 0
        {
            bam.set_threads(config.htslib_threads)?;
        }

        // cache the expanded index
        let bam_index: Vec<(Vec<u8>, u64, u64)> = bam.expanded_index()?;
        Ok(VariantCounter {
            config,
            bam,
            bam_index
        })
    }

    pub fn index(&self) -> &Vec<(Vec<u8>, u64, u64)>
    {
        &self.bam_index
    }

    /// Generate a closure that can be used to filter alignments, given the configuration settings
    /// This is implemented as a class and not a member function to avoid mutable/immutable ref issues
    fn generate_alignment_filter<'a>(config: &'a VariantCounterConfig) -> impl Fn(&Alignment<'a>) -> bool
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

            if qual < config.min_baseq
            {
                return false;
            }

            let read_pair_orientation = record.read_pair_orientation_lenient(config.exclude_ambiguous);

            match read_pair_orientation
            {
                F1R2 =>
                {
                    if record.is_first_in_template()
                    {
                        // Ensure that there's at least one base left after soft-trimming
                        if seq_len < config.ot_mask.r1.0+config.ot_mask.r1.1 + 1 {
                            return false;
                        }

                        if qpos < config.ot_mask.r1.0 || qpos > seq_len-config.ot_mask.r1.1-1
                        {
                            return false;
                        }
                    }
                    else
                    {
                        if seq_len < config.ot_mask.r2.0+config.ot_mask.r2.1 + 1 {
                            return false;
                        }
                        // also flipped the end/start mask, cause the read is mapped in reverse
                        if qpos < config.ot_mask.r2.1 || qpos > seq_len-config.ot_mask.r2.0-1
                        {
                            return false;
                        }
                    }
                },
                F2R1 =>
                {
                    if record.is_first_in_template()
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
                    else
                    {
                        if seq_len < config.ob_mask.r2.0+config.ob_mask.r2.1 + 1 {
                            return false;
                        }
                        if qpos < config.ob_mask.r2.0 || qpos > seq_len-config.ob_mask.r2.1-1
                        {
                            return false;
                        }
                    }
                },
                _   =>
                {
                    warn!("Unexpected read orientation or mates on different chromosomes for record {}", String::from_utf8(Vec::from(record.qname())).unwrap_or_default());
                    return false;
                },
            };

            true
        };
        filter_closure
    }

    /// For a given sequence segment, return the number of observed nucleotide counts at all covered
    /// positions
    pub fn count_variants_in_segment(&mut self, segment: SequenceSegment) -> Option<Vec<VariantCount>>
    {
        debug!("Processing segment {}", &segment);
        //TODO this needs changing to make it more generic, ie allow different types of subsets
        // Search the string segment for all CpG positions
        let Some(cpg_positions) = segment.find_cpgs() else {
            info!("No CpGs in {}", &segment);
            return None;
        };

        /* Fetch the pileup for the region from the bam file, and go
        * through all CpG positions, performing whatever calculation
        * needs to be performed. Stream the output to a writer that
        * writes the results to STDOUT or some file.
        */
        match self.bam.fetch((segment.region.contig.as_bytes(), segment.region.start, segment.region.end))
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
        let mut read_hash:HashMap<Vec<u8>, (u8, u8), FxBuildHasher> = HashMap::with_capacity_and_hasher(self.config.max_depth as usize, FxBuildHasher::default());

        // Find the next CpG position that is greater than or equal to the current pos
        // If none exists, we've reached the end
        let mut cpg_index = 0;
        let mut last_pos = 0;
        let mut repeat_err_count: u8 = 0;
        'pileup_loop:
        for pileups in pileup_iterator
        {
            let pileup = match pileups
            {
                Ok(p) => p,
                Err(e)  =>
                {
                    // check if we're in a death-loop
                    if repeat_err_count < 3
                    {
                        repeat_err_count += 1;
                        warn!("Error reading pileup at {} {}", &segment, e);
                        continue;
                    }
                    else
                    {
                        error!("Error reading pileup at {} {}", &segment, e);
                        return None;
                    }

                }
            };
            // reset error counter
            repeat_err_count = 0;

            let pileup_pos = pileup.pos() as u64;
            last_pos = pileup_pos;
            trace!("start: {} pileup_pos: {}", segment.region.start, pileup.pos());
            let mut this_position = &cpg_positions[cpg_index];

            if pileup_pos < this_position.pos_in_contig()
            {
                // Alignment behind
                continue;
            }
            else if pileup_pos > this_position.pos_in_contig() {
                // Skipped some positions because there were no reads,
                // ie cpg_pos behind
                loop {
                    cpg_index = cpg_index+1;
                    if cpg_index >= cpg_positions.len()
                    {
                        break 'pileup_loop;
                    }
                    this_position = &cpg_positions[cpg_index];
                    if this_position.pos_in_contig() == pileup_pos
                    {
                        break;
                    }
                    else if this_position.pos_in_contig() > pileup_pos
                    {
                        continue 'pileup_loop;
                    }
                }
            }

            debug!("Found CpG site {}, col: {}", this_position, pileup_pos);
            cpg_index = cpg_index + 1;

            // Count conversion vs non-conversion
            let mut var_count = VariantCount::new();
            var_count.contig = segment.region.contig.clone();
            var_count.pos = pileup_pos;
            var_count.ref_base = this_position.base();
            // Need to be done here so that the lifetime of pileup exceeds the lifetime of the closure
            let filter_closure = VariantCounter::generate_alignment_filter(&self.config);
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
                debug!("Processing pos {} in read {} with base {} and qual {} at col {}", pos, std::str::from_utf8(record.qname()).unwrap_or_default(), char::from_u32(base as u32).unwrap_or_default(), qual, pileup_pos);
                let read_pair_orientation = record.read_pair_orientation_lenient(self.config.exclude_ambiguous);

                let nuc_counts =
                {
                    match read_pair_orientation {
                        F1R2 => &mut var_count.top,
                        F2R1 => &mut var_count.bottom,
                        _   => {
                            error!("Cannot process ambiguous read-pair, should have been filtered earlier!");
                            continue 'alignment_loop;
                        }
                    }
                };

                // Check for overlapping reads if requested
                if !self.config.keep_overlaps
                {
                    if let Some(pos_tuple) = read_hash.get(record.qname())
                    {
                        trace!("Found overlapping read pair {} at pos {} with bases {} vs {}", std::str::from_utf8(record.qname()).unwrap_or_default(), pos, char::from_u32(base as u32).unwrap_or_default(), char::from_u32(pos_tuple.0 as u32).unwrap_or_default());
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
                            nuc_counts.increment_counter_by(pos_tuple.0, -1);
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

                match nuc_counts.increment_counter_by(base, 1)
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
            // only report positions with at least _some_ coverage
            if var_count.total_count() > 0
            {
                output.push(var_count);
            }
            // Check if we're done
            if cpg_index >= cpg_positions.len()
            {
                break 'pileup_loop;
            }
        }
        if cpg_index < cpg_positions.len()
        {
            warn!("Not all CpG positions were processed! Should have fetched {}, but stopped at {}. {} CpGs skipped.", &segment, last_pos, cpg_positions.len()-cpg_index);
            debug!("Next CpG in the list: {}", cpg_positions[cpg_index]);
        }
        Some(output)
    }
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

    fn create_test_config() -> Result<VariantCounterConfig>
    {
        //let fasta_path = PathBuf::from(r"test_data/test.fasta");
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
        let mut config = create_test_config()?;
        config.max_depth = 1000;
        let mut bam = IndexedReader::from_path(&config.bam_path)?;

        bam.fetch(("bacteriophage_lambda_CpG", 0, 100))?;
        let mut count_all = 0;
        let mut count_filt = 0;
        for p in bam.pileup()
        {
            if let Ok(pileup) = p {
                let filter = VariantCounter::generate_alignment_filter(&config);
                let alignments: Vec<Alignment> = pileup.alignments().collect();
                count_all += alignments.iter().count();
                count_filt += alignments.into_iter().filter(&filter).count();
            }
        }
        assert_eq!(count_all, 1227);
        assert_eq!(count_filt, 1224);

        Ok(())
    }

    #[test]
    fn can_filter_with_masking_1() -> Result<()>
    {
        let mut config = create_test_config()?;
        config.ot_mask = ReadMaskSetting::from_str(&"1,0,0,1").unwrap();
        config.ob_mask = ReadMaskSetting::from_str(&"0,1,1,0").unwrap();
        let filter = VariantCounter::generate_alignment_filter(&config);
        let mut bam = IndexedReader::from_path(&config.bam_path)?;

        bam.fetch(("bacteriophage_lambda_CpG", 0, 100))?;
        if let Ok(pileup) = bam.pileup().nth(0).unwrap() {
            let alignments: Vec<Alignment> = pileup.alignments().collect();
            assert_eq!(alignments.len(), 6);
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
        let filter = VariantCounter::generate_alignment_filter(&config);
        let mut bam = IndexedReader::from_path(&config.bam_path)?;

        bam.fetch(("bacteriophage_lambda_CpG", 1000, 1200))?;
        if let Ok(pileup) = bam.pileup().nth(100).unwrap() {
            let alignments: Vec<Alignment> = pileup.alignments().collect();
            assert_eq!(alignments.len(), 10);
            assert_eq!(alignments.into_iter().filter(filter).count(), 9);
        };

        Ok(())
    }
}