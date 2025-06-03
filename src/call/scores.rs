// Calculate the following scores:
// - binomial test probability of "false positive", given background error likelihoods
// - Root Mean Square (RMS) of mapping quality (mapQ) of reads evidence for reference allele/alternative allele
// - RMS of base quality (baseQ) for reference allele/alternative allele
// - Chi-square test of position-in-reads of alternative alleles being uniformly distributed (e.g. bin 150bp into 10 regions for Chi-square test)
// - @RMS of variant positions per read (excluding likely @methylation), for reads supporting @RefAllele and @AltAllele
// - Entropy of sequence in region (i.e. whether the @variant is in a @repeatedSequence)
// - Strand bias between @OT and @OB for @RefAllele and @AltAllele
// - Tri- or penta-nucleotide context#note[consider 2 or 4 bases immediately flanking the variant (one on each side)] of variant for @SNV:pl
// - Repeat/homopolymer length of region for @indel:pl
// - Realignment score/hamming-distance difference for reads covering @variant#mark[assumes realignment]

use crate::utils::RootMeanSquare;
use color_eyre::eyre::{Context, Result};
use probability::{distribution::Binomial, prelude::Discrete};
use smallvec::SmallVec;
use std::fmt;

use super::variants::VariantCandidatePileup;

impl VariantCandidatePileup {
    pub fn metrics(&self) -> Result<VariantCandidatePileupMetrics> {
        let reference_bases = self.bases.iter().filter(|b| b.base == self.reference_base);
        let reference_count = reference_bases.clone().count();
        let alt_bases = self.bases.iter().filter(|b| b.base != self.reference_base);
        let alt_count = alt_bases.clone().count();

        let vaf = VariantAlleleFrequency {
            reference_count: reference_count as u64,
            alt_count: alt_count as u64,
        };

        let binomial = BinomialTest {
            reference_count: reference_count as u64,
            alt_count: alt_count as u64,
            error_rate: 0.01, // fixme: use real error rates
        };

        let mapq = MappingQuality {
            reference_mapq: SmallVec::from_iter(reference_bases.clone().map(|b| b.mapq)),
            alt_mapq: SmallVec::from_iter(alt_bases.clone().map(|b| b.mapq)),
        };

        let baseq = BaseQuality {
            reference_baseq: SmallVec::from_iter(reference_bases.clone().map(|b| b.qual)),
            alt_baseq: SmallVec::from_iter(alt_bases.clone().map(|b| b.qual)),
        };

        Ok(VariantCandidatePileupMetrics {
            reference_count,
            alt_count,
            vaf: vaf.calculate().wrap_err("vaf")?,
            binomial: binomial.calculate().wrap_err("binomial")?,
            mapq: mapq.calculate().wrap_err("mapq")?,
            baseq: baseq.calculate().wrap_err("baseq")?,
            strand_bias: StrandBias {
                reference_ot: reference_bases.clone().filter(|b| !b.reverse).count() as u64,
                reference_ob: reference_bases.clone().filter(|b| b.reverse).count() as u64,
                alt_ot: alt_bases.clone().filter(|b| !b.reverse).count() as u64,
                alt_ob: alt_bases.clone().filter(|b| b.reverse).count() as u64,
            }
            .calculate()
            .wrap_err("strand bias")?,
        })
    }
}

/// Metrics for a variant candidate based on its pileup
#[derive(Debug)]
pub struct VariantCandidatePileupMetrics {
    pub reference_count: usize,
    pub alt_count: usize,
    /// Variant Allele Frequency, ratio of alternative allele to total alleles
    pub vaf: f64,
    /// Probability of "false positive" given a specified read error rate.
    /// Values close to 1 indicate a high probability of false positive.
    /// Typically seen values are between 0.0001 and 0.2.
    pub binomial: f64,
    /// RMS of mapping quality for (reference allele, alternative allele)
    pub mapq: RefVsAlt<RootMeanSquare>,
    /// RMS of base quality (baseQ) for (reference allele, alternative allele)
    pub baseq: RefVsAlt<RootMeanSquare>,
    /// Strand bias between OT and OB for (reference allele, alternative-allele)
    pub strand_bias: RefVsAlt<f64>,
}

/// Simple wrapper for metrics for reference and alternative alleles
pub struct RefVsAlt<T> {
    pub reference: T,
    pub alt: T,
}

impl<T> RefVsAlt<T> {
    pub fn new(reference: T, alt: T) -> Self {
        RefVsAlt { reference, alt }
    }
}

#[cfg(not(tarpaulin_include))]
impl<T: fmt::Debug> fmt::Debug for RefVsAlt<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "(ref: {:?}, alt: {:?})", self.reference, self.alt)
    }
}

pub trait Calc {
    type Output;
    fn calculate(&self) -> Result<Self::Output>;
}

pub struct VariantAlleleFrequency {
    pub reference_count: u64,
    pub alt_count: u64,
}

impl Calc for VariantAlleleFrequency {
    type Output = f64;

    fn calculate(&self) -> Result<f64> {
        let total = self.reference_count + self.alt_count;
        if total == 0 {
            return Ok(0.0);
        }
        Ok(self.alt_count as f64 / total as f64)
    }
}

pub struct BinomialTest {
    pub reference_count: u64,
    pub alt_count: u64,
    pub error_rate: f64,
}

impl Calc for BinomialTest {
    type Output = f64;

    fn calculate(&self) -> Result<f64> {
        let total = self.reference_count + self.alt_count;
        if total == 0 {
            return Ok(0.0);
        }
        let binomial =
            Binomial::new(usize::try_from(total).wrap_err("total > usize")?, self.error_rate);

        if self.reference_count <= self.alt_count {
            Ok(binomial
                .mass(usize::try_from(self.reference_count).wrap_err("reference_count > usize")?))
        } else {
            Ok(binomial.mass(usize::try_from(self.alt_count).wrap_err("alt_count > usize")?))
        }
    }
}

// Root Mean Square (RMS) of mapping quality (mapQ) of reads evidence for reference allele/alternative allele
pub struct MappingQuality {
    pub reference_mapq: SmallVec<u8, 16>,
    pub alt_mapq: SmallVec<u8, 16>,
}

impl Calc for MappingQuality {
    type Output = RefVsAlt<RootMeanSquare>;

    fn calculate(&self) -> Result<Self::Output> {
        let reference_rms = RootMeanSquare::new(&self.reference_mapq);
        let alt_rms = RootMeanSquare::new(&self.alt_mapq);
        Ok(RefVsAlt::new(reference_rms, alt_rms))
    }
}

// RMS of base quality (baseQ) for reference allele/alternative allele
pub struct BaseQuality {
    pub reference_baseq: SmallVec<u8, 16>,
    pub alt_baseq: SmallVec<u8, 16>,
}

impl Calc for BaseQuality {
    type Output = RefVsAlt<RootMeanSquare>;

    fn calculate(&self) -> Result<Self::Output> {
        let reference_rms = RootMeanSquare::new(&self.reference_baseq);
        let alt_rms = RootMeanSquare::new(&self.alt_baseq);
        Ok(RefVsAlt::new(reference_rms, alt_rms))
    }
}

// Strand bias between OT and OB for reference allele and alternative-allele
pub struct StrandBias {
    /// Counts of bases matching the reference base on the forward strand
    pub reference_ot: u64,
    /// Counts of bases matching the reference base on the reverse strand
    pub reference_ob: u64,
    /// Counts of bases matching the alternative base on the forward strand
    pub alt_ot: u64,
    /// Counts of bases matching the alternative base on the reverse strand
    pub alt_ob: u64,
}

impl Calc for StrandBias {
    type Output = RefVsAlt<f64>;

    fn calculate(&self) -> Result<Self::Output> {
        let reference_total = self.reference_ot + self.reference_ob;
        let alt_total = self.alt_ot + self.alt_ob;

        let reference_bias = if reference_total == 0 {
            0.0
        } else {
            self.reference_ot as f64 / reference_total as f64
        };

        let alt_bias = if alt_total == 0 { 0.0 } else { self.alt_ot as f64 / alt_total as f64 };

        Ok(RefVsAlt::new(reference_bias, alt_bias))
    }
}
