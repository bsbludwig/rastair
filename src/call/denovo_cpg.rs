use crate::vcf;
use color_eyre::Result;

#[derive(Debug, Clone, clap::Args)]
pub struct DenovoParams {
    /// Minimum reads needed in support of de-novo CpG
    #[clap(long, default_value_t = 2)]
    pub cpg_novo_min_depth: u32,
    /// Minimum base quality for de-novo CpGs
    #[clap(long, default_value_t = 15)]
    pub cpg_novo_min_baseq: u32,
    /// Minimum mapping quality for de-novo CpGs
    #[clap(long, default_value_t = 50)]
    pub cpg_novo_min_mapq: u32,
    /// Minimum variant allele frequency for de-novo CpGs
    #[clap(long, default_value_t = 0.2)]
    pub cpg_novo_min_vaf: f64,
}

impl DenovoParams {
    pub fn filter(&self, record: &mut vcf::Record) -> Result<()> {
        if !record.info.de_novo_cp_g_candidate.0 {
            return Ok(());
        }

        // TODO: Add filtering logic based on the parameters
        record.samples[0].de_novo_cpg = vcf::DeNovoCpg(Some(0.));
        Ok(())
    }
}
