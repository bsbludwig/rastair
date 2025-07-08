use crate::{utils::Base, vcf};
use std::fmt;

/// Helper struct to hold values for top and bottom strands
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ByStrand<T> {
    /// Base of the allele
    pub base: Base,
    /// Value for the top strand
    pub ot: T,
    /// Value for the bottom strand
    pub ob: T,
}

impl<T: Copy> Copy for ByStrand<T> {}

impl<T: Default> Default for ByStrand<T> {
    fn default() -> Self {
        ByStrand { base: Base::Unknown, ot: T::default(), ob: T::default() }
    }
}

/// Helper methods for more concise queries
impl vcf::Record {
    /// Returns true if the record has the given base as an alternative allele
    pub fn has_alt(&self, base: Base) -> bool {
        self.main.alt.iter().any(|alt| alt == &base)
    }

    /// Returns true if the record has any alternative alleles other than the given base
    pub fn has_alts_other_than(&self, base: Base) -> bool {
        self.main.alt.iter().any(|alt| alt != &base)
    }

    /// Returns the allele frequency for the given base
    pub fn strand_count(&self, base: Base) -> Result<vcf::ByStrand<u32>, NoStrandBiasForBaseError> {
        self.info
            .allele_specific_strand_bias
            .iter()
            .find(|x| x.base == base)
            .copied()
            .ok_or(NoStrandBiasForBaseError { base })
    }
}

#[derive(Debug, thiserror::Error)]
pub struct NoStrandBiasForBaseError {
    base: Base,
}

impl fmt::Display for NoStrandBiasForBaseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "No strand bias information for base {}", self.base)
    }
}

pub trait NoStrandBiasForBaseErrorExt {
    /// Returns the strand bias counts for the base, or a default value if the error is encountered.
    fn or_empty(&self) -> vcf::ByStrand<u32>;
}

impl NoStrandBiasForBaseErrorExt for Result<vcf::ByStrand<u32>, NoStrandBiasForBaseError> {
    fn or_empty(&self) -> vcf::ByStrand<u32> {
        match self {
            Ok(counts) => *counts,
            Err(NoStrandBiasForBaseError { base }) => vcf::ByStrand { base: *base, ot: 0, ob: 0 },
        }
    }
}
