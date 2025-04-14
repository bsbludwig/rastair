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

use probability::{distribution::Binomial, prelude::Discrete};
use smallvec::SmallVec;

pub struct VariantAlleleFrequency {
    pub reference_count: u64,
    pub alt_count: u64,
}

impl VariantAlleleFrequency {
    pub fn calculate(&self) -> f64 {
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

impl BinomialTest {
    pub fn calculate(&self) -> f64 {
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

impl MappingQuality {
    pub fn calculate(&self) -> (f64, f64) {
        let reference_rms = rms(&self.reference_mapq);
        let alt_rms = rms(&self.alt_mapq);
        (reference_rms, alt_rms)
    }
}

fn rms(mapq: &[u8]) -> f64 {
    if mapq.is_empty() {
        return 0.0;
    }
    let sum_of_squares: f64 = mapq.iter().map(|&x| (x as f64).powi(2)).sum();
    (sum_of_squares / mapq.len() as f64).sqrt()
}
