//! De-novo CpG calling and filtering
//!
//! A "de-novo CpG" is a pair of adjacent nucleotides where in the reference
//! genome they are not a CpG (i.e. not "C" followed by "G"), but in the sample
//! they are changed to a CpG by a SNP (i.e. a reference "C" followed by a base
//! that changed from "A" to "G").

use crate::{
    call::pileup::Pileup,
    utils::{Base::*, cli},
    vcf::{self},
};
use better_default::Default;
use color_eyre::Result;

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct DenovoParams {
    /// Minimum reads needed in support of de-novo CpG
    #[arg(long, default_value_t = 2)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(2)]
    pub cpg_novo_min_depth: usize,
    /// Minimum base quality for de-novo CpGs
    #[arg(long, default_value_t = 15.)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(15.)]
    pub cpg_novo_min_baseq: f64,
    /// Minimum mapping quality for de-novo CpGs
    #[arg(long, default_value_t = 50.)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(50.)]
    pub cpg_novo_min_mapq: f64,
    /// Minimum variant allele frequency for de-novo CpGs
    #[arg(long, default_value_t = 0.2)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(0.2)]
    pub cpg_novo_min_vaf: f64,
}

impl From<&Pileup> for vcf::DeNovoCpGCandidate {
    fn from(pileup: &Pileup) -> Self {
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
    fn filter_vcf(&self, record: &mut vcf::Record) -> Result<()> {
        let vcf::DeNovoCpGCandidate::Candidate { alt_base, alt_index: idx, .. } =
            record.info.de_novo_cp_g_candidate
        else {
            return Ok(());
        };

        if record.info.allele_read_depth.get(idx) <= Some(&self.cpg_novo_min_depth) {
            record.filters.add_per_allele(alt_base, vcf::dnCpG_lowDp);
        }

        if record.info.allele_base_quality.get(idx) <= Some(&self.cpg_novo_min_baseq) {
            record.filters.add_per_allele(alt_base, vcf::dnCpG_bq);
        }

        if record.info.allele_map_quality.get(idx) <= Some(&self.cpg_novo_min_mapq) {
            record.filters.add_per_allele(alt_base, vcf::dnCpG_mapq);
        }

        if record.info.allele_frequency.get(idx) <= Some(&self.cpg_novo_min_vaf) {
            record.filters.add_per_allele(alt_base, vcf::dnCpG_vaf);
        }

        Ok(())
    }
}

// #[cfg(test)]
// mod tests {
//     use crate::call::{test_helpers::variant_pileup, variant_calling::VariantCallingParams};

//     use super::*;

//     #[test]
//     #[ignore = "wip"]
//     fn denovo_position() -> Result<()> {
//         let now_c = variant_pileup("bacteriophage_lambda_CpG", 2199)?
//             .variant_metrics(&VariantCallingParams::default())?;
//         let ref_g = variant_pileup("bacteriophage_lambda_CpG", 2200)?
//             .variant_metrics(&VariantCallingParams::default())?;

//         DenovoParams::default().add_if_adjecent(&mut now_c.clone(), None, Some(&ref_g));
//         assert!(matches!(
//             now_c.info.de_novo_cp_g_candidate,
//             vcf::DeNovoCpGCandidate::Adjecent { ref_base: C }
//         ));

//         Ok(())
//     }
// }
