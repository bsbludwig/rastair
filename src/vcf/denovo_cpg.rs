use crate::utils::Base;
use color_eyre::{Result, eyre::WrapErr};
use rastair_vcf::{InfoFieldNumber, VcfField as _};
use rust_htslib::bcf::Record;
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
        /// The index of the alternative allele that creates a new CpG site
        alt_index: usize,
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

impl DeNovoCpGCandidate {
    pub fn alt_base(&self) -> Option<Base> {
        match self {
            DeNovoCpGCandidate::NotCandidate => None,
            DeNovoCpGCandidate::Candidate { alt_base, .. } => Some(*alt_base),
            DeNovoCpGCandidate::Adjecent { ref_base } => Some(*ref_base),
        }
    }
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

impl rastair_vcf::VcfField for DeNovoCpGCandidate {
    const ID: &'static cstr8::CStr8 = cstr8::cstr8!("CPGnovo");
}

impl rastair_vcf::HeaderField for DeNovoCpGCandidate {
    const DESCRIPTION: &'static str =
        "De-novo CPG candidate: Could the alt alleles create a new CpG site?";
}

impl rastair_vcf::InfoField for DeNovoCpGCandidate {
    type Type = bool;
    const NUMBER: InfoFieldNumber = InfoFieldNumber::Flag;

    fn write(&self, record: &mut Record) -> Result<()> {
        let is_candidate = matches!(self, DeNovoCpGCandidate::Candidate { .. });
        <bool as rastair_vcf::InfoFieldValue>::write(record, Self::ID, &[is_candidate])
            .wrap_err("Failed to write info flag CPGnovo")
    }
}
