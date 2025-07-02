use crate::{
    call::variants::VariantCandidatePileup,
    utils::Base::*,
    vcf::{self},
};
use color_eyre::Result;

#[derive(Debug, Clone, clap::Args)]
pub struct DenovoParams {
    /// Minimum reads needed in support of de-novo CpG
    #[clap(long, default_value_t = 2)]
    pub cpg_novo_min_depth: usize,
    /// Minimum base quality for de-novo CpGs
    #[clap(long, default_value_t = 15.)]
    pub cpg_novo_min_baseq: f64,
    /// Minimum mapping quality for de-novo CpGs
    #[clap(long, default_value_t = 50.)]
    pub cpg_novo_min_mapq: f64,
    /// Minimum variant allele frequency for de-novo CpGs
    #[clap(long, default_value_t = 0.2)]
    pub cpg_novo_min_vaf: f64,
}

impl From<&VariantCandidatePileup> for vcf::DeNovoCpGCandidate {
    fn from(pileup: &VariantCandidatePileup) -> Self {
        let ref_base = pileup.reference_base;
        let (alt_base, alt_index) = if let Some(pos) = pileup.alts().iter().position(|x| *x == C)
            && pileup.ref_after() == G
        {
            (C, pos)
        } else if let Some(pos) = pileup.alts().iter().position(|x| *x == G)
            && pileup.ref_before() == C
        {
            (G, pos)
        } else {
            return vcf::DeNovoCpGCandidate::NotCandidate;
        };

        vcf::DeNovoCpGCandidate::Candidate { ref_base, alt_base, alt_index }
    }
}

impl DenovoParams {
    pub fn filter(&self, record: &mut vcf::Record) -> Result<()> {
        let vcf::DeNovoCpGCandidate::Candidate { alt_index: idx, .. } =
            record.info.de_novo_cp_g_candidate
        else {
            return Ok(());
        };

        if record.info.allele_read_depth.get(idx) <= Some(&self.cpg_novo_min_depth) {
            record.filters.add(vcf::dnCpG_lowDp);
        }

        if record.info.allele_base_quality.get(idx) <= Some(&self.cpg_novo_min_baseq) {
            record.filters.add(vcf::dnCpG_bq);
        }

        if record.info.allele_map_quality.get(idx) <= Some(&self.cpg_novo_min_mapq) {
            record.filters.add(vcf::dnCpG_mapq);
        }

        if record.info.allele_frequency.get(idx) <= Some(&self.cpg_novo_min_vaf) {
            record.filters.add(vcf::dnCpG_vaf);
        }

        Ok(())
    }
}
