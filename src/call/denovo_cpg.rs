//! De-novo CpG calling and filtering
//!
//! A "de-novo CpG" is a pair of adjacent nucleotides where in the reference
//! genome they are not a CpG (i.e. not "C" followed by "G"), but in the sample
//! they are changed to a CpG by a SNP (i.e. a reference "C" followed by a base
//! that changed from "A" to "G").

use crate::{
    call::{methylation::ThresholdParams, pileup::Pileup},
    metrics::{Filters, MetricsForAlt},
    utils::{Base::*, IntoF64 as _, cli},
    vcf,
};
use better_default::Default;
use color_eyre::Result;

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct DenovoParams {
    /// Minimum reads needed in support of de-novo CpG
    #[arg(long, default_value_t = 2)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(2)]
    pub cpg_novo_min_depth: u32,
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

pub fn add_filters(config: &DenovoParams, current: &MetricsForAlt) -> Result<Filters> {
    let mut filters = Filters::default();

    let alt = &current.alt;

    if !*alt.denovo {
        // Not a de-novo CpG site, skipping filters
        return Ok(filters);
    }

    filters.add(vcf::dnCpG_lowDp, || alt.depth < config.cpg_novo_min_depth);
    filters.add(vcf::dnCpG_bq, || alt.baseq.f() < config.cpg_novo_min_baseq);
    filters.add(vcf::dnCpG_mapq, || alt.mapq.f() < config.cpg_novo_min_mapq);
    filters.add(vcf::dnCpG_vaf, || alt.allele_frequency.f() < config.cpg_novo_min_vaf);

    Ok(filters)
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
