use crate::{
    call::variants::VariantCandidatePileup,
    utils::Base,
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
        let (alt_base, alt_index) = if let Some(pos) =
            pileup.alts().iter().position(|x| *x == Base::C)
            && pileup.ref_after() == Some(Base::G)
        {
            (Base::C, pos)
        } else if let Some(pos) = pileup.alts().iter().position(|x| *x == Base::G)
            && pileup.ref_before() == Some(Base::C)
        {
            (Base::G, pos)
        } else {
            return vcf::DeNovoCpGCandidate::NotCandidate;
        };

        vcf::DeNovoCpGCandidate::Candidate { ref_base, alt_base, alt_index }
    }
}

impl DenovoParams {
    pub fn filter(
        &self,
        record: &mut vcf::Record,
        _before: Option<&vcf::Record>,
        _after: Option<&vcf::Record>,
    ) -> Result<()> {
        let vcf::DeNovoCpGCandidate::Candidate { alt_index, .. } =
            record.info.de_novo_cp_g_candidate
        else {
            return Ok(());
        };

        let critera: [(bool, Box<dyn Fn(&mut vcf::Record)>); 4] = [
            (
                record.info.allele_read_depth.get(alt_index).copied().unwrap_or_default()
                    >= self.cpg_novo_min_depth,
                Box::new(|record| record.filters.add(vcf::dnCpG_lowDp)),
            ),
            (
                record.info.allele_base_quality.get(alt_index).copied().unwrap_or_default()
                    >= self.cpg_novo_min_baseq,
                Box::new(|record| record.filters.add(vcf::dnCpG_bq)),
            ),
            (
                record.info.allele_map_quality.get(alt_index).copied().unwrap_or_default()
                    >= self.cpg_novo_min_mapq,
                Box::new(|record| record.filters.add(vcf::dnCpG_mapq)),
            ),
            (
                record.info.allele_frequency.get(alt_index).copied().unwrap_or_default()
                    >= self.cpg_novo_min_vaf,
                Box::new(|record| record.filters.add(vcf::dnCpG_vaf)),
            ),
        ];

        let critera_len = critera.len();
        let mut met_criteria = 0;
        for (criteria_met, filter_fn) in critera {
            if criteria_met {
                met_criteria += 1;
            } else {
                filter_fn(record);
            }
        }

        record.samples[0].de_novo_cpg =
            vcf::DeNovoCpg(Some(f64::from(met_criteria) / critera_len as f64));
        Ok(())
    }
}
