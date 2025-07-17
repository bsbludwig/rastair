use crate::{
    utils::Base,
    vcf::{self},
};

impl vcf::Record {
    /// Returns the base of the methylation evidence
    pub fn m_base(&self) -> Option<Base> {
        self.info.in_cp_g.alt_base().or_else(|| self.info.de_novo_cp_g_candidate.alt_base())
    }
}
