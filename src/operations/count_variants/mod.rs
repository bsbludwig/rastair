pub mod variant_counter;

use std::str::FromStr;
use std::fmt::{Debug, Display, Formatter};
use std::error::Error;
use std::io::{stdout, Write};
use std::fs::File;
use std::path::{Path, PathBuf};
use log::{debug, info, warn, error};

use thiserror::Error;
use anyhow::Result;
use r2d2::{ManageConnection, Pool};
use pariter::IteratorExt as _;

use crate::sequence_segment::{SequenceSegmentIterator, SequenceSegment};

pub use super::{ReadMaskSetting, ReadMask};
use super::{FLAGS, MAX_DEPTH};
use variant_counter::{VariantCounter, VariantCounterConfig};

/// A simple struct to represent counts of nucleotides
pub struct NucleotideCount
{
    pub a: i32,
    pub c: i32,
    pub g: i32,
    pub t: i32,
    pub n: i32,
}

impl Display for NucleotideCount
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result
    {
        write!(f, "A: {} C: {} G: {} T: {} N: {}", self.a, self.c, self.g, self.t, self.n)
    }
}

/// A representation of a variant position
pub struct VariantCount
{
    /// ID of the sequence
    pub contig: String,
    /// position in sequence coordinate space
    pub pos: u64,
    /// Nucleotide in the reference sequence
    pub ref_base: u8,
    /// Counts of nucleotides observed on the OB
    pub top: NucleotideCount,
    /// Counts of nucleotides observed on the OT
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

impl Display for VariantCount
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result
    {
        let char = char::from_u32(self.ref_base as u32).unwrap_or_default();
        write!(f, "{}:{} ({})\nFW\t{}\nRV\t{}", self.contig, self.pos, char, self.top, self.bottom)
    }
}

struct VariantCounterConnectionManager<P>
{
    config: VariantCounterConfig<P>
}

impl <P: AsRef<Path> + Clone + Debug> VariantCounterConnectionManager<P>
{
    fn with_config(config: VariantCounterConfig<P>) -> Result<Self>
    {
        Ok(VariantCounterConnectionManager{
            config
        })
    }
}

#[derive(Error, Debug)]
pub enum VariantCounterConnectionError {
    #[error("Error connecting to the bam file")]
    ConnectionError( #[from] anyhow::Error )
}

impl <P: AsRef<Path> + Clone + Debug + std::marker::Send + std::marker::Sync + 'static> ManageConnection for VariantCounterConnectionManager<P>
{
    type Connection = VariantCounter<P>;
    // TODO create a proper custom error type
    type Error = VariantCounterConnectionError;

    fn connect(&self) -> Result<Self::Connection, Self::Error> {
        match VariantCounter::with_config(self.config.clone())
        {
            Ok(counter) => Ok(counter),
            Err(e)  => Err(VariantCounterConnectionError::ConnectionError(e))
        }
    }

    fn is_valid(&self, _conn: &mut Self::Connection) -> Result<(), Self::Error> {
        //TODO better check for valid connection?
        Ok(())
    }

    fn has_broken(&self, _conn: &mut Self::Connection) -> bool {
        false
    }
}

#[allow(non_snake_case)]
pub fn run_caller(
    bam_path: &PathBuf,
    fasta_path: &PathBuf,
    mapq_option: &Option<u8>,
    baseq_option: &Option<u8>,
    max_depth_option: &Option<u32>,
    chunk_size_option: &Option<u32>,
    req_flags_option: &Option<u16>,
    excl_flags_option: &Option<u16>,
    nOT_option: &Option<String>,
    nOB_option: &Option<String>,
    read_threads_option: &Option<u8>,
    threads_option: &Option<u8>) -> Result<(), Box<dyn Error>> 
{
    /* Read fasta index, and open fasta file for tokenising */
    debug!("Reading fasta and index from {}", fasta_path.display());

    let chunk_size = if let Some(cs) = chunk_size_option {
        *cs as usize
    } 
    else 
    {
        0
    };

    let threads = if let Some(t) = threads_option {
        *t as usize
    } 
    else 
    {
        1
    };

    let mut config = VariantCounterConfig::with_path(bam_path.clone())?;
    if let Some(min_mapq) = mapq_option {
        config.min_mapq = *min_mapq;
    }
    if let Some(min_baseq) = baseq_option {
        config.min_baseq = *min_baseq;
    }
    if let Some(max_depth) = max_depth_option {
        config.max_depth = *max_depth;
    }
    
    if let Some(flags) = req_flags_option {
        config.required_flags = *flags;
    }
    if let Some(flags) = excl_flags_option {
        config.excluded_flags = *flags;
    }
    if let Some(threads) = read_threads_option {
        config.htslib_threads = *threads as usize;
    }
    
    #[allow(non_snake_case)]
    if let Some(nOT_s) = nOT_option {
        if let Ok(ot_mask) = ReadMaskSetting::from_str(nOT_s) {
            config.ot_mask = ot_mask;
        }
    }
    #[allow(non_snake_case)]
    if let Some(nOB_s) = nOB_option {
        if let Ok(ob_mask) = ReadMaskSetting::from_str(nOB_s) {
            config.ob_mask = ob_mask;
        }
    }

    let manager = VariantCounterConnectionManager::with_config(config)?;
    let pool = r2d2::Pool::builder()
        .max_size(threads as u32)
        .build(manager)?;
    //let mut counter = VariantCounter::with_config(config)?;
    
    /* TODO this needs refactoring to allow multi-threading:
     * 1a. Loop over genomic segments, and then inject the segment into the VariantCounter [x]
     * 1b. Ensure GenomicSegments are cloneable, so I can distribute them? Or will a thread scope suffice?
     * 2. Change the VariantCounter to not be an iterator, but create a separate VariantCounterIterator
     *    that calls a generic method to extract the nucleotide counts [x]
     * 3. Create a pool of VariantCounters, each with a fixed open bam file, and use one of them per thread,
     *    using e.g. [R2D2](https://docs.rs/r2d2/latest/r2d2/trait.ManageConnection.html) [x]
     * 4. Use e.g. [pariter](https://lib.rs/crates/pariter) and fetch a VariantCounter (with attached bam handle)
     *    for each closure invokation
    */
    let mut iterator: SequenceSegmentIterator<File> = 
        if chunk_size == 0 
        { 
            SequenceSegmentIterator::with_file(fasta_path)?
        } 
        else 
        { 
            SequenceSegmentIterator::with_file_and_stepsize(fasta_path, chunk_size)?
        };
    let counter = pool.get()?;
    iterator.subset_to_intervals(counter.index())?;
    drop(counter);

    // Get a write lock on STDOUT
    let mut lock = stdout().lock();

    iterator
    .map(move |segment| {
        SegmentCounterPair {
            segment,
            pool: pool.clone()
        }
    })
    .parallel_map_custom(|b| b.threads(threads), |sp|
    {
        let segment = sp.segment;
        debug!("Will try to count variants in {}", segment);
        let mut counter = sp.pool.get().unwrap();
        debug!("Counter address: {:p}", &counter);
        counter.count_variants_in_segment(segment).unwrap()
    })
    .flatten()
    .for_each(|cpg| 
    {
        if cpg.ref_base == b'C'
        { // C
            writeln!(lock, "{}\t{}\t{}\t.\t.\t+\t{}\t{}\t{}\t{}\t{}", cpg.contig, cpg.pos, cpg.pos+1, cpg.top.c, cpg.top.t, cpg.bottom.c, cpg.bottom.t, cpg.top.a + cpg.top.c + cpg.top.g + cpg.top.t + cpg.bottom.a + cpg.bottom.c + cpg.bottom.g + cpg.bottom.t).unwrap();
        }
        else
        { // G
            writeln!(lock, "{}\t{}\t{}\t.\t.\t-\t{}\t{}\t{}\t{}\t{}", cpg.contig, cpg.pos, cpg.pos+1, cpg.bottom.g, cpg.bottom.a, cpg.top.g, cpg.top.a, cpg.top.a + cpg.top.c + cpg.top.g + cpg.top.t + cpg.bottom.a + cpg.bottom.c + cpg.bottom.g + cpg.bottom.t).unwrap();
        }
    });

    Ok(())
}

struct SegmentCounterPair {
    segment: SequenceSegment,
    pool: Pool<VariantCounterConnectionManager<PathBuf>>
}