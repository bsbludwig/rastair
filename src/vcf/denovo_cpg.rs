use crate::utils::Base;
use color_eyre::{Result, eyre::WrapErr};
use rastair2_vcf::{InfoFieldNumber, VcfField as _};
use rust_htslib::bcf::Record;
use std::ops::Deref;

/// De-novo CPG candidate: Could the alt alleles create a new CpG site?
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub enum DeNovoCpGCandidate {
    NotCandidate,
    Candidate { ref_base: Base, alt_base: Base, alt_index: usize },
}

impl DeNovoCpGCandidate {
    pub fn alt_base(&self) -> Option<Base> {
        match self {
            DeNovoCpGCandidate::NotCandidate => None,
            DeNovoCpGCandidate::Candidate { alt_base, .. } => Some(*alt_base),
        }
    }
}

impl Deref for DeNovoCpGCandidate {
    type Target = bool;

    fn deref(&self) -> &Self::Target {
        match self {
            DeNovoCpGCandidate::NotCandidate => &false,
            DeNovoCpGCandidate::Candidate { .. } => &true,
        }
    }
}

impl rastair2_vcf::VcfField for DeNovoCpGCandidate {
    const ID: &'static cstr8::CStr8 = cstr8::cstr8!("CPGnovo");
}

impl rastair2_vcf::HeaderField for DeNovoCpGCandidate {
    const DESCRIPTION: &'static str =
        "De-novo CPG candidate: Could the alt alleles create a new CpG site?";
}

impl rastair2_vcf::InfoField for DeNovoCpGCandidate {
    type Type = bool;
    const NUMBER: InfoFieldNumber = InfoFieldNumber::Flag;

    fn write(&self, record: &mut Record) -> Result<()> {
        let is_candidate = matches!(self, DeNovoCpGCandidate::Candidate { .. });
        <bool as rastair2_vcf::InfoFieldValue>::write(record, Self::ID, &[is_candidate])
            .wrap_err("Failed to write info flag CPGnovo")
    }
}
