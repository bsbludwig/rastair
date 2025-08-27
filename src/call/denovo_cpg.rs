//! De-novo CpG calling and filtering
//!
//! A "de-novo CpG" is a pair of adjacent nucleotides where in the reference
//! genome they are not a CpG (i.e. not "C" followed by "G"), but in the sample
//! they are changed to a CpG by a SNP (i.e. a reference "C" followed by a base
//! that changed from "A" to "G").

use crate::{
    call::variants::VariantCandidatePileup,
    utils::Base::*,
    vcf::{self},
};
use better_default::Default;
use color_eyre::Result;

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct DenovoParams {
    /// Minimum reads needed in support of de-novo CpG
    #[clap(long, default_value_t = 2)]
    #[default(2)]
    pub cpg_novo_min_depth: usize,
    /// Minimum base quality for de-novo CpGs
    #[clap(long, default_value_t = 15.)]
    #[default(15.)]
    pub cpg_novo_min_baseq: f64,
    /// Minimum mapping quality for de-novo CpGs
    #[clap(long, default_value_t = 50.)]
    #[default(50.)]
    pub cpg_novo_min_mapq: f64,
    /// Minimum variant allele frequency for de-novo CpGs
    #[clap(long, default_value_t = 0.2)]
    #[default(0.2)]
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

    pub fn add_if_adjecent(
        &self,
        current: &mut vcf::Record,
        before: Option<&vcf::Record>,
        after: Option<&vcf::Record>,
    ) {
        if *current.info.de_novo_cp_g_candidate {
            // already a candidate
        } else if let Some(before) = before
            && let vcf::DeNovoCpGCandidate::Candidate { alt_base, .. } =
                before.info.de_novo_cp_g_candidate
            && alt_base == C
            && current.main.r#ref == G
        {
            // previous position has a C alt, so this G is now part of a CpG
            current.info.de_novo_cp_g_candidate = vcf::DeNovoCpGCandidate::Adjecent { ref_base: G };
        } else if let Some(after) = after
            && let vcf::DeNovoCpGCandidate::Candidate { alt_base, .. } =
                after.info.de_novo_cp_g_candidate
            && alt_base == G
            && current.main.r#ref == C
        {
            // next position has a G alt, so this C is now part of a CpG
            current.info.de_novo_cp_g_candidate = vcf::DeNovoCpGCandidate::Adjecent { ref_base: C };
        }
        // todo: do both AdjecentRef and AdjecentAlt?
        // ("alt" means both positions create the cpg with alts… but then who created it and who is adjecent?)
    }
}

#[cfg(test)]
mod tests {
    use crate::call::{test_helpers::variant_pileup, variant_calling::VariantCallingParams};

    use super::*;

    #[test]
    #[ignore = "wip"]
    fn denovo_position() -> Result<()> {
        let now_c = variant_pileup("bacteriophage_lambda_CpG", 2199)?
            .variant_metrics(&VariantCallingParams::default())?;
        let ref_g = variant_pileup("bacteriophage_lambda_CpG", 2200)?
            .variant_metrics(&VariantCallingParams::default())?;

        DenovoParams::default().add_if_adjecent(&mut now_c.clone(), None, Some(&ref_g));
        assert!(matches!(
            now_c.info.de_novo_cp_g_candidate,
            vcf::DeNovoCpGCandidate::Adjecent { ref_base: C }
        ));

        Ok(())
    }
}
