use crate::{
    call::variants::VariantCandidatePileup,
    utils::{Base, Base::*},
};
use color_eyre::Result;
use color_eyre::eyre::WrapErr;
use rastair_vcf::VcfField as _;
use rust_htslib::bcf::Record;
use std::fmt;
use std::ops::Deref;

/// Is this a CpG site?
#[derive(Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum InCpG {
    /// No, this is not a CpG site.
    No,
    /// Yes, first position in CpG.
    C,
    /// Yes, second position in CpG.
    G,
}

impl InCpG {
    pub fn new(base: Base, before: Option<Base>, after: Option<Base>) -> Self {
        if base == C && after == G {
            InCpG::C
        } else if base == G && before == C {
            InCpG::G
        } else {
            InCpG::No
        }
    }

    pub fn alt_base(&self) -> Option<Base> {
        match self {
            InCpG::C => Some(T),
            InCpG::G => Some(A),
            InCpG::No => None,
        }
    }
}

impl fmt::Debug for InCpG {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InCpG::No => write!(f, "NoCpg"),
            InCpG::C => write!(f, "CpG::C"),
            InCpG::G => write!(f, "CpG::G"),
        }
    }
}

impl Deref for InCpG {
    type Target = bool;

    fn deref(&self) -> &Self::Target {
        match self {
            InCpG::No => &false,
            _ => &true,
        }
    }
}

impl From<&VariantCandidatePileup> for InCpG {
    fn from(pileup: &VariantCandidatePileup) -> Self {
        InCpG::new(pileup.reference_base, pileup.ref_before(), pileup.ref_after())
    }
}

impl rastair_vcf::VcfField for InCpG {
    const ID: &'static cstr8::CStr8 = cstr8::cstr8!("CPG");
}

impl rastair_vcf::HeaderField for InCpG {
    const DESCRIPTION: &'static str = "Is this a CpG site?";
}

impl rastair_vcf::InfoField for InCpG {
    type Type = bool;
    const NUMBER: rastair_vcf::InfoFieldNumber = rastair_vcf::InfoFieldNumber::Flag;

    fn write(&self, record: &mut Record) -> Result<()> {
        <bool as rastair_vcf::InfoFieldValue>::write(record, Self::ID, &[*self != InCpG::No])
            .wrap_err("Failed to write info flag CPG")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_in_cpg() {
        assert_eq!(InCpG::new(A, None, None), InCpG::No);

        assert_eq!(InCpG::new(C, None, Some(G)), InCpG::C);
        assert_eq!(InCpG::new(C, Some(A), Some(G)), InCpG::C);

        assert_eq!(InCpG::new(G, Some(C), None), InCpG::G);
        assert_eq!(InCpG::new(G, Some(C), Some(T)), InCpG::G);

        assert_eq!(InCpG::new(T, Some(C), Some(G)), InCpG::No);
        assert_eq!(InCpG::new(C, Some(T), Some(A)), InCpG::No);
    }
}
