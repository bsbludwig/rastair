pub mod variant_counter;

use std::str::FromStr;
use std::fmt::{Debug, Display, Formatter};
use std::error::Error;
use std::io::{stdout, Write};
use std::path::PathBuf;
use clap::ValueEnum;
use log::{debug, error};

use num_cpus;
use thiserror::Error;
use anyhow::{bail, Result};
use r2d2::ManageConnection;
use pariter::IteratorExt as _;
use probability::prelude::*;

use crate::sequence_segment::SequenceSegmentIterator;
use crate::utils::file_helpers::open_file;
use super::{ReadMaskSetting, ReadMask};
use crate::utils::constants::*;
use variant_counter::{VariantCounter, VariantCounterConfig};

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum ErrorModel {
    Miseq,
    Miniseq,
    Nextseq500,
    Nextseq550,
    Hiseq2500,
    Novaseq6000,
    HiseqXTen
}

/// A simple struct to represent counts of nucleotides
pub struct NucleotideCount
{
    pub a: i32,
    pub c: i32,
    pub g: i32,
    pub t: i32,
    pub n: i32,
}

impl NucleotideCount {
    pub fn total(&self) -> i32
    {
        return self.a + self.c + self.g + self.t + self.n;
    }

    fn increment_counter_by(&mut self, base: u8, amount: i32) -> Option<()>
{
    match base
    {
        b'a' => self.a += amount,
        b'c' => self.c += amount,
        b'g' => self.g += amount,
        b't' => self.t += amount,
        b'n' => self.n += amount,
        b'A' => self.a += amount,
        b'C' => self.c += amount,
        b'G' => self.g += amount,
        b'T' => self.t += amount,
        b'N' => self.n += amount,
        _   =>
        {
            return None;
        }
    };
    Some(())
}
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

    pub fn total_count(&self) -> i32
    {
        return self.top.total() + self.bottom.total();
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

struct VariantCounterConnectionManager
{
    config: VariantCounterConfig
}

impl VariantCounterConnectionManager
{
    fn with_config(config: VariantCounterConfig) -> Result<Self>
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

impl ManageConnection for VariantCounterConnectionManager
{
    type Connection = VariantCounter;
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
    region_option: &Option<String>,
    mapq_option: &Option<u8>,
    baseq_option: &Option<u8>,
    max_depth_option: &Option<u32>,
    chunk_size_option: &Option<u32>,
    req_flags_option: &Option<u16>,
    excl_flags_option: &Option<u16>,
    error_model_option: &Option<ErrorModel>,
    exclude_ambiguous_option: &Option<bool>,
    nOT_option: &Option<String>,
    nOB_option: &Option<String>,
    read_threads_option: &Option<u8>,
    threads_option: &Option<u8>) -> Result<(), Box<dyn Error>>
{
    /* Read fasta index, and open fasta file for tokenising */
    debug!("Reading fasta and index from {}", fasta_path.display());

    let chunk_size = chunk_size_option.unwrap_or_default();

    let max_threads = num_cpus::get();
    let threads = std::cmp::min(threads_option.unwrap_or(1) as usize, max_threads);

    let mut config = VariantCounterConfig::with_path(bam_path.clone())?;
    if let Some(min_mapq) = *mapq_option {
        config.min_mapq = min_mapq;
    }
    if let Some(min_baseq) = *baseq_option {
        config.min_baseq = min_baseq;
    }
    if let Some(max_depth) = *max_depth_option {
        config.max_depth = max_depth;
    }
    if let Some(exclude_ambiguous) = *exclude_ambiguous_option {
        config.exclude_ambiguous = exclude_ambiguous;
    }
    if let Some(flags) = *req_flags_option {
        config.required_flags = flags;
    }
    if let Some(flags) = *excl_flags_option {
        config.excluded_flags = flags;
    }
    if let Some(threads) = *read_threads_option {
        config.htslib_threads = threads as usize;
    }

    // pass this right into the config object, will de-parse when needed
    config.region = region_option.to_owned();

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
    if let Some(error_model) = error_model_option {
        match error_model {
            ErrorModel::Miseq => config.error_model = ERRORRATES.miseq,
            ErrorModel::Miniseq => config.error_model = ERRORRATES.miniseq,
            ErrorModel::Nextseq500 => config.error_model = ERRORRATES.nextseq_500,
            ErrorModel::Nextseq550 => config.error_model = ERRORRATES.nextseq_550,
            ErrorModel::Hiseq2500 => config.error_model = ERRORRATES.hiseq_2500,
            ErrorModel::Novaseq6000 => config.error_model = ERRORRATES.novaseq_6000,
            ErrorModel::HiseqXTen => config.error_model = ERRORRATES.hiseq_x_ten,
        }
    }
    // neet to do here, before move into manager
    let error_model = config.error_model;

    let manager = VariantCounterConnectionManager::with_config(config)?;
    let pool = r2d2::Pool::builder()
        .max_size(threads as u32)
        .build(manager)?;
    //let mut counter = VariantCounter::with_config(config)?;

    // Load fasta file as an indexedreader
    // need to do this manually to enable bgzip-compressed input
    let fasta_file = open_file(&fasta_path)?;
    let index_path = PathBuf::from(format!("{}.fai", fasta_path.to_str().unwrap_or_default()).as_str());
    let fasta_index = bio::io::fasta::Index::from_file(&index_path)?;
    let indexed_reader = bio::io::fasta::IndexedReader::with_index(fasta_file, fasta_index);

    let mut iterator =
        if chunk_size == 0
        {
            SequenceSegmentIterator::with_reader(indexed_reader)?
        }
        else
        {
            SequenceSegmentIterator::with_reader_and_stepsize(indexed_reader, chunk_size as usize)?
        };

    // Optionally restrict to region of interest
    if let Some(region) = region_option
    {
        iterator.subset_to_region(region)?;
    }
    // Subset to those sequences that are actually in the bam/fasta file
    let counter = pool.get()?;
    iterator.subset_to_intervals(counter.index())?;
    drop(counter); // Ugly, but I need to free that counter up for later


    // Get a write lock on STDOUT
    let mut lock = stdout().lock();
    writeln!(lock, "#chr\tstart\tend\tname\tbeta_est\tstrand\tunmod\tmod\tno_snp\tsnp\tcoverage\tgenotype\tgt_p_score\tgt_conf_score")?;
    iterator
    .map(move |segment| {
        (segment, pool.clone())
    })
    .parallel_map_custom(|b| b.threads(threads), |sp|
    {
        let segment = sp.0;
        let pool = sp.1;
        debug!("Will try to count variants in {}", segment);
        let Ok(mut counter) = pool.get() else
        {
            panic!("Failed to get counter from pool, too many threads?");
        };

        match counter.count_variants_in_segment(segment)
        {
            Some(res) => res,
            None => Vec::new()
        }
    })
    .flatten()
    .for_each(|cpg|
    {
        let (unmod_c, mod_c, nosnp, snp, strand) =
        if cpg.ref_base == b'C'
        {
            (cpg.top.c, cpg.top.t, cpg.bottom.c, cpg.bottom.t, "+")
        }
        else
        {
            (cpg.bottom.g, cpg.bottom.a, cpg.top.g, cpg.top.a, "-")
        };
        let gt = EstimatedGenotype::calculate(nosnp, snp, error_model).unwrap_or_default();
        let beta: f32 = match gt.genotype {
            Genotype::CC => (mod_c as f32)/(mod_c + unmod_c) as f32,
            Genotype::CT => ((mod_c as f32)/2.0)/((mod_c as f32)/2.0 + unmod_c as f32),
            Genotype::TT => 0.0,
        };
        let gt_string = if cpg.ref_base == b'C'
        { // C
            match gt.genotype {
                Genotype::CC  => "C/C",
                Genotype::CT  => "C/T",
                Genotype::TT  => "T/T",
            }
        }
        else
        { // G
            match gt.genotype {
                Genotype::CC  => "G/G",
                Genotype::CT  => "G/A",
                Genotype::TT  => "A/A",
            }
        };
        writeln!(lock, "{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", cpg.contig, cpg.pos, cpg.pos+1, ".", beta, strand, unmod_c, mod_c, nosnp, snp, cpg.top.a + cpg.top.c + cpg.top.g + cpg.top.t + cpg.bottom.a + cpg.bottom.c + cpg.bottom.g + cpg.bottom.t, gt_string, prob_to_phred(1.0-gt.likelihood), prob_to_phred(1.0-gt.confidence)).unwrap();
    });

    Ok(())
}

fn prob_to_phred (prob: f64) -> u8
{
    let phred = -10.0 * prob.log10();
    if phred >= 99.0
    {
        return 99;
    }
    if phred <= f64::MIN_POSITIVE
    {
        return 0;
    }
    else
    {
        return phred.round() as u8;
    }
}

/// Public only because it's exposed in ErrorRates
pub type ErrorRate = f64;
/// Empirically derived error rates as published here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8002175/
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct ErrorRates {
    miseq: ErrorRate,
    miniseq: ErrorRate,
    nextseq_500: ErrorRate,
    nextseq_550: ErrorRate,
    hiseq_2500: ErrorRate,
    novaseq_6000: ErrorRate,
    hiseq_x_ten: ErrorRate
}

const ERRORRATES: ErrorRates = ErrorRates {
    miseq: 0.00473,
    miniseq: 0.00613,
    nextseq_500: 0.00429,
    nextseq_550: 0.00593,
    hiseq_2500: 0.00112,
    novaseq_6000: 0.00109,
    hiseq_x_ten: 0.00087
};

enum Genotype {
    CC,
    CT,
    TT
}
struct EstimatedGenotype {
    genotype: Genotype,
    likelihood: f64,
    confidence: f64
}
impl Default for EstimatedGenotype {
    fn default() -> Self {
        Self { genotype: Genotype::CC, likelihood: 0.0, confidence: 0.0 }
    }
}

impl EstimatedGenotype {
    fn calculate(ref_count: i32, alt_count: i32, error_rate: ErrorRate) -> Result<Self>
{
    if ref_count == alt_count && alt_count == 0
    {
        bail!("No ref or alt read counts, cannot compute likelihood");
    }
    if error_rate <= f64::MIN
    {
        bail!("Error rate too small, cannot calculate likelihood");
    }
    // This is a simple estimate of genotype, based on the following consideration:
    // A site is either het or hom, where hom could be CC or TT.
    // If alt_count > ref_count, the latter is more likely, otherwise the former.

    // First, I calculate the likelihood to observe this many alt_reads
    // under the assumption that ref and alt are equally likely, ie this is a het position.
    // TODO This assumes a simple diploid sample with no purity issues. For
    // cancer samples, we could make this a setting to allow for different cancer fraction?

    let mut binom = Binomial::new((ref_count + alt_count) as usize, 0.5); // 0.5 because a het position
    let p_het = binom.mass(alt_count as usize);
    let p_het_max = binom.mass(((alt_count+ref_count) as f32 / 2.0).round() as usize);

    // Then, I calculate the probability that this many or more alt_count/ref_count reads
    // are observed by error, assuming independence of reads and errors.
    binom = Binomial::new((ref_count + alt_count) as usize, error_rate);

    if ref_count >= alt_count {
        let p_hom = binom.mass(alt_count as usize) + (1.0 - binom.distribution(alt_count as f64));

        if p_het < p_hom
        {
            debug!("Assuming CC: ({} vs {}) -> ({:.5} < {:.5})", ref_count, alt_count, p_het, p_hom);
            return Ok(EstimatedGenotype {
                genotype: Genotype::CC,
                likelihood: p_hom,
                confidence: (p_hom - p_het)/p_hom
            });
        }
        else
        {
            debug!("Assuming CT: ({} vs {}) -> ({:.5} >= {:.5})", ref_count, alt_count, p_het, p_hom);
            return Ok(EstimatedGenotype {
                genotype: Genotype::CT,
                likelihood: p_het / p_het_max,
                confidence: (p_het - p_hom)/p_het
            });
        }
    }
    else
    {
        let p_hom = binom.mass(ref_count as usize) + (1.0-binom.distribution(ref_count as f64));
        if p_het < p_hom
        {
            debug!("Assuming TT: ({} vs {}) -> ({:.5} < {:.5})", ref_count, alt_count, p_het, p_hom);
            return Ok(EstimatedGenotype {
                genotype: Genotype::TT,
                likelihood: p_hom,
                confidence: (p_hom - p_het)/p_hom
            });
        }
        else
        {
            debug!("Assuming TC: ({} vs {}) -> ({:.5} >= {:.5})", ref_count, alt_count, p_het, p_hom);
            return Ok(EstimatedGenotype {
                genotype: Genotype::CT,
                likelihood: p_het / p_het_max,
                confidence: (p_het - p_hom)/p_het
            });
        }
    }
}
}