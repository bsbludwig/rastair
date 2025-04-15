// Calculate the following scores:
// - Variant Allele Frequency
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

use std::fmt;

use probability::{distribution::Binomial, prelude::Discrete};
use smallvec::SmallVec;

use crate::utils::RootMeanSquare;

pub struct RefVsAlt<T> {
    reference: T,
    alt: T,
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
    fn calculate(&self) -> Self::Output;
}

pub struct VariantAlleleFrequency {
    pub reference_count: u64,
    pub alt_count: u64,
}

impl Calc for VariantAlleleFrequency {
    type Output = f64;

    fn calculate(&self) -> f64 {
        let total = self.reference_count + self.alt_count;
        if total == 0 {
            return 0.0;
        }
        self.alt_count as f64 / total as f64
    }
}

pub struct BinomialTest {
    pub reference_count: u64,
    pub alt_count: u64,
    pub error_rate: f64,
}

impl Calc for BinomialTest {
    type Output = f64;

    fn calculate(&self) -> f64 {
        let total = self.reference_count + self.alt_count;
        if total == 0 {
            return 0.0;
        }
        let binomial = Binomial::new(total as usize, self.error_rate);

        if self.reference_count <= self.alt_count {
            binomial.mass(self.reference_count as usize)
        } else {
            binomial.mass(self.alt_count as usize)
        }
    }
}

// - Root Mean Square (RMS) of mapping quality (mapQ) of reads evidence for reference allele/alternative allele
pub struct MappingQuality {
    pub reference_mapq: SmallVec<u8, 16>,
    pub alt_mapq: SmallVec<u8, 16>,
}

impl Calc for MappingQuality {
    type Output = RefVsAlt<RootMeanSquare>;

    fn calculate(&self) -> Self::Output {
        let reference_rms = RootMeanSquare::new(&self.reference_mapq);
        let alt_rms = RootMeanSquare::new(&self.alt_mapq);
        RefVsAlt::new(reference_rms, alt_rms)
    }
}

// RMS of base quality (baseQ) for reference allele/alternative allele
pub struct BaseQuality {
    pub reference_baseq: SmallVec<u8, 16>,
    pub alt_baseq: SmallVec<u8, 16>,
}

impl Calc for BaseQuality {
    type Output = RefVsAlt<RootMeanSquare>;

    fn calculate(&self) -> Self::Output {
        let reference_rms = RootMeanSquare::new(&self.reference_baseq);
        let alt_rms = RootMeanSquare::new(&self.alt_baseq);
        RefVsAlt::new(reference_rms, alt_rms)
    }
}

// Strand bias between OT and OB for reference allele and alternative-allele
pub struct StrandBias {
    pub reference_ot: u64,
    pub reference_ob: u64,
    pub alt_ot: u64,
    pub alt_ob: u64,
}

impl Calc for StrandBias {
    type Output = RefVsAlt<f64>;

    fn calculate(&self) -> Self::Output {
        let reference_total = self.reference_ot + self.reference_ob;
        let alt_total = self.alt_ot + self.alt_ob;

        let reference_bias = if reference_total == 0 {
            0.0
        } else {
            self.reference_ot as f64 / reference_total as f64
        };

        let alt_bias = if alt_total == 0 { 0.0 } else { self.alt_ot as f64 / alt_total as f64 };

        RefVsAlt::new(reference_bias, alt_bias)
    }
}
