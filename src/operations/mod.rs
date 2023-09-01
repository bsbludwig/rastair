use rust_htslib::bam::pileup::Alignment;
use rust_htslib::bam::{IndexedReader, Read};
use std::{fmt, fs, path::Path};
use log::{trace, debug, warn, error};
use anyhow::Result;

use crate::sequence_segment::SequenceSegmentIterator;
const MAX_DEPTH: u32 = 500;
struct Flags
{
    is_paired: u16,
    is_properly_paired: u16,
    is_unmapped: u16,
    mate_is_unmapped: u16,
    is_reverse_strand: u16,
    mate_is_reverse_strand: u16,
    is_first_in_pair: u16,
    is_second_in_pair: u16,
    is_not_primary: u16,
    is_failed: u16,
    is_duplicate: u16,
    is_supplemental: u16,
}
const FLAGS: Flags = Flags
{
    is_paired: 0x1,
    is_properly_paired: 0x2,
    is_unmapped: 0x4,
    mate_is_unmapped: 0x8,
    is_reverse_strand: 0x10,
    mate_is_reverse_strand: 0x20,
    is_first_in_pair: 0x40,
    is_second_in_pair: 0x80,
    is_not_primary: 0x100,
    is_failed: 0x200,
    is_duplicate: 0x400,
    is_supplemental: 0x800
};

pub struct NucleotideCount
{
    pub a: u32,
    pub c: u32,
    pub g: u32,
    pub t: u32,
    pub n: u32,
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
    pub chunk_size: u64,
    pub required_flags: u16,
    pub excluded_flags: u16
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
            chunk_size: 10000,
            required_flags: FLAGS.is_paired | FLAGS.is_properly_paired,
            excluded_flags: FLAGS.is_failed | FLAGS.is_not_primary | FLAGS.is_unmapped | FLAGS.mate_is_unmapped | FLAGS.is_duplicate | FLAGS.is_supplemental,
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
        let bam = IndexedReader::from_path(&config.bam_path)?;
        let fasta = SequenceSegmentIterator::with_file_and_stepsize(&config.fasta_path, config.chunk_size)?;

        Ok(
            VariantCounter
            {
                config,
                bam,
                fasta
            }
        )
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
            return None;
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

        'pileuploop:
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
            // Find the next CpG position that is greater than or equal to the current pos
            // If none exists, we've reached the end
            let this_position = match cpg_positions
            .iter()
            .find(|pos| pos.pos_in_contig() >= pileup_pos)
            {
                Some(pos)   => pos,
                None =>
                {
                    debug!("Found all CpGs, next segment");
                    break 'pileuploop
                }
            };
            if pileup_pos == this_position.pos_in_contig()
            {
                debug!("Found a CpG site at {}", this_position);

                // Count conversion vs non-conversion
                let mut var_count = VariantCount::new();
                var_count.contig = segment.contig.clone();
                var_count.pos = this_position.pos_in_contig();
                var_count.ref_base = this_position.base();

                let filter_closure = |alignment: &Alignment| -> bool
                {
                    let mut filter = !(alignment.is_del() || alignment.is_refskip());
                    filter &= alignment.record().mapq() >= self.config.min_mapq;
                    // Require all "required" flags (properly paired/aligned)
                    filter &= (alignment.record().flags() & self.config.required_flags) == self.config.required_flags;
                    // exclude reads matching _any_ excluded flag
                    filter &= (alignment.record().flags() & self.config.excluded_flags) == 0;
                    filter
                };
                for alignment in
                                                pileup
                                                .alignments()
                                                .filter(filter_closure)
                {
                    let qual = match qual_at_position(&alignment)
                    {
                        Some(b) => b,
                        None => continue
                    };

                    if qual < self.config.min_baseq
                    {
                        continue;
                    }

                    let base = match base_at_position(&alignment)
                    {
                        Some(b) => b,
                        None => continue
                    };
                    let nuc_counts =
                    {
                        if (alignment.record().flags() & (FLAGS.is_first_in_pair | FLAGS.mate_is_reverse_strand) == 0)
                        || (alignment.record().flags() & (FLAGS.is_second_in_pair | FLAGS.is_reverse_strand) == 0)
                        {
                            &mut var_count.top
                        }
                        else
                        {
                            &mut var_count.bottom
                        }
                    };

                    match base
                    {
                        b'a' => nuc_counts.a += 1,
                        b'c' => nuc_counts.c += 1,
                        b'g' => nuc_counts.g += 1,
                        b't' => nuc_counts.t += 1,
                        b'n' => nuc_counts.n += 1,
                        b'A' => nuc_counts.a += 1,
                        b'C' => nuc_counts.c += 1,
                        b'G' => nuc_counts.g += 1,
                        b'T' => nuc_counts.t += 1,
                        b'N' => nuc_counts.n += 1,
                        _   =>
                        {
                            let char = char::from_u32(base as u32).unwrap_or_default();
                            warn!("Encountered unknown char {} at {}", char, this_position);
                        }
                    }
                }
                debug!("{}", var_count);
                output.push(var_count);
            }
        }

        Some(output)
    }
}

fn qual_at_position(alignment: &Alignment) -> Option<u8>
{
    let pos = alignment.qpos()?;
    if alignment.is_del()
    {
        return None;
    }
    let record = alignment.record();
    Some(record.qual()[pos])
}

fn base_at_position(alignment: &Alignment) -> Option<u8>
{
    let pos = alignment.qpos()?;
    if alignment.is_del()
    {
        return None;
    }
    let record = alignment.record();
    let seq = record.seq();
    if seq.len() == 0
    {
        return None;
    }
    Some(seq[pos])
}
