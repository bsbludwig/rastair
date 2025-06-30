use crate::vcf;
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

impl DenovoParams {
    pub fn filter(
        &self,
        record: &mut vcf::Record,
        before: Option<&vcf::Record>,
        after: Option<&vcf::Record>,
    ) -> Result<()> {
        if !record.info.de_novo_cp_g_candidate.0 {
            return Ok(());
        }

        let alt_idx = {
            if let Some(pos) = record.fixed_fields.alt.iter().position(|alt| alt == "C")
                && let Some(after) = after
                && after.fixed_fields.r#ref == "G"
            {
                pos
            } else if let Some(pos) = record.fixed_fields.alt.iter().position(|alt| alt == "G")
                && let Some(before) = before
                && before.fixed_fields.r#ref == "C"
            {
                pos
            } else {
                return Ok(());
            }
        };

        let critera = [
            record.info.read_depth_per_allel.get(alt_idx).copied().unwrap_or_default()
                >= self.cpg_novo_min_depth,
            record.info.allel_base_quality.get(alt_idx).copied().unwrap_or_default()
                >= self.cpg_novo_min_baseq,
            record.info.allel_map_quality.get(alt_idx).copied().unwrap_or_default()
                >= self.cpg_novo_min_mapq,
            record.info.allel_frequency.get(alt_idx).copied().unwrap_or_default()
                >= self.cpg_novo_min_vaf,
        ];
        let criteria_met: u32 = critera.iter().map(|&c| u32::from(c)).sum();

        // TODO: Add filtering logic based on the parameters
        record.samples[0].de_novo_cpg =
            vcf::DeNovoCpg(Some(f64::from(criteria_met) / critera.len() as f64));
        Ok(())
    }
}
