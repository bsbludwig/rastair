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
