use crate::{
    utils::Base,
    vcf::{self, DeNovoCpGCandidate, InCpG},
};

impl vcf::Record {
    /// Returns the base of the methylation evidence
    pub fn m_base(&self) -> Option<Base> {
        if self.info.in_cp_g == InCpG::C {
            Some(Base::T)
        } else if self.info.in_cp_g == InCpG::G {
            Some(Base::C)
        } else if let DeNovoCpGCandidate::Candidate { alt_base, .. } =
            self.info.de_novo_cp_g_candidate
        {
            Some(alt_base)
        } else {
            // We're not looking at a CpG site or de-novo CpG
            None
        }
    }
}
