use crate::utils::Base;
use std::ops::Deref;

/// De-novo CPG candidate: Could the alt alleles create a new CpG site?
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub enum DeNovoCpGCandidate {
    /// No, the variant cannot create a new CpG site
    #[default]
    NotCandidate,
    /// Yes, the variant could create a new CpG site
    Candidate {
        /// The reference base at the variant position
        ref_base: Base,
        /// The alternative base that creates a new CpG site
        alt_base: Base,
    },
    /// Yes, the variant is adjacent to a de-novo CpG site
    Adjecent {
        /// The reference base at this adjacent position
        ///
        /// If the other position has a C variant, this is G; and if the other
        /// position has a G variant, this is C.
        ref_base: Base,
    },
}

impl Deref for DeNovoCpGCandidate {
    type Target = bool;

    fn deref(&self) -> &Self::Target {
        match self {
            DeNovoCpGCandidate::NotCandidate => &false,
            DeNovoCpGCandidate::Candidate { .. } => &true,
            DeNovoCpGCandidate::Adjecent { .. } => &true,
        }
    }
}
