use crate::{utils::Base, vcf};
use color_eyre::Result;
use tracing::warn;

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

impl DenovoParams {
    pub fn filter(
        &self,
        record: &mut vcf::Record,
        _before: Option<&vcf::Record>,
        _after: Option<&vcf::Record>,
    ) -> Result<()> {
        if !record.info.de_novo_cp_g_candidate.0 {
            return Ok(());
        }

        let alt_idx = {
            if let Some(pos) = record.fixed_fields.alt.iter().position(|alt| alt == "C")
                && record.info.sequence_context.after_1 == Some(Base::G)
            {
                pos
            } else if let Some(pos) = record.fixed_fields.alt.iter().position(|alt| alt == "G")
                && record.info.sequence_context.before_1 == Some(Base::C)
            {
                pos
            } else {
                warn!(
                    "Position previously tagged as de-novo CpG candidate does not have C or G as alternate allele. This is likely a bug in the variant calling pipeline."
                );
                return Ok(());
            }
        };

        let critera: [(bool, Box<dyn Fn(&mut vcf::Record)>); 4] = [
            (
                record.info.read_depth_per_allel.get(alt_idx).copied().unwrap_or_default()
                    >= self.cpg_novo_min_depth,
                Box::new(|record| record.filters.add(vcf::dnCpG_lowDp)),
            ),
            (
                record.info.allel_base_quality.get(alt_idx).copied().unwrap_or_default()
                    >= self.cpg_novo_min_baseq,
                Box::new(|record| record.filters.add(vcf::dnCpG_bq)),
            ),
            (
                record.info.allel_map_quality.get(alt_idx).copied().unwrap_or_default()
                    >= self.cpg_novo_min_mapq,
                Box::new(|record| record.filters.add(vcf::dnCpG_mapq)),
            ),
            (
                record.info.allel_frequency.get(alt_idx).copied().unwrap_or_default()
                    >= self.cpg_novo_min_vaf,
                Box::new(|record| record.filters.add(vcf::dnCpG_af)),
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
