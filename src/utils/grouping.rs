use crate::utils::Base;

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

/// Helper struct to group values by allele
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ByAllele<T> {
    /// Base of the allele
    pub base: Base,
    /// Value for the allele
    pub value: T,
}

impl<T> std::ops::Deref for ByAllele<T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        &self.value
    }
}

impl<T: Copy> Copy for ByAllele<T> {}

impl<T: Default> Default for ByAllele<T> {
    fn default() -> Self {
        ByAllele { base: Base::Unknown, value: T::default() }
    }
}
