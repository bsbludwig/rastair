use crate::{
    bed::rastair1::{BedRecordsConvertParams, Rastair1BedFormat},
    utils::{Base::*, logging::ThisIsABug as _},
    vcf::{Record as Rastair2Record, utils::NoStrandBiasForBaseErrorExt as _},
};
use color_eyre::{Result, eyre::Context as _};
use rastair_types::Probability;
use tracing::trace;

impl Rastair1BedFormat {
    #[allow(clippy::cast_possible_truncation)]
    pub fn from_record(
        record: &Rastair2Record,
        params: &BedRecordsConvertParams,
    ) -> Result<Option<Self>> {
        // Skip positions without evidence
        if !params.filters.include_empty && *record.info.read_depth == 0 {
            return Ok(None);
        }

        // If a position is covered by both a ref CpG site and a de-novo CpG
        // site, the ref case should take precedence.
        let de_novo = !*record.info.in_cp_g && *record.info.de_novo_cp_g_candidate;
        if de_novo && !record.filters.pass() {
            // Only report de novo candidates that we are confident about
            trace!(
                pos=%record.main,
                ml=?record.samples[0].machine_learning_prediction,
                "de novo candidate with low score"
            );
            return Ok(None);
        }

        let r#ref = record.main.r#ref.clone();

        let (unmod, r#mod, no_snp, snp) = if r#ref == "C" {
            (
                record.strand_count(C).or_empty().ot,
                record.strand_count(T).or_empty().ot,
                record.strand_count(C).or_empty().ob,
                record.strand_count(T).or_empty().ob,
            )
        } else if r#ref == "G" {
            (
                record.strand_count(G).or_empty().ob,
                record.strand_count(A).or_empty().ob,
                record.strand_count(G).or_empty().ot,
                record.strand_count(A).or_empty().ot,
            )
        } else {
            (0, 0, 0, 0)
        };

        // If this looks like SNP based on the genotype calling, set beta to 0
        // TODO: make ML threshold configurable
        let beta = if let Some(alt) = record.info.in_cp_g.alt_base()
            && let Some(alt_idx) = record.main.alt.iter().position(|b| b == &alt)
            && let Some(ml_score) = record.samples[0].machine_learning_prediction.get(alt_idx)
            && *ml_score > *params.ml_threshold
        {
            // ML score indicates this is a SNP
            trace!(%record.main, %ml_score, "SNP detected at CpG site based on ML score");
            Some(0.0)
        } else if *record.info.in_cp_g && record.samples[0].genotype.homozygous_not_ref() {
            Some(0.0)
        } else {
            record.samples[0].methylated.beta()
        };
        let beta = if let Some(beta) = beta {
            Some(Probability::new(beta).wrap_err("Beta value out of range").this_is_a_bug()?)
        } else {
            None
        };

        Ok(Some(super::Rastair1BedFormat {
            contig: record.main.chrom.clone(),
            pos: record.main.pos as usize,
            r#ref,
            beta,
            unmod,
            r#mod,
            no_snp,
            snp,
            coverage: *record.info.read_depth,
            genotype: record.samples[0].genotype.clone(),
            genotype_likelihood: record.samples[0].genotype_likelihood.clone(),
            genotype_confidence: record.samples[0].genotype_confidence.clone(),
            de_novo,
        }))
    }
}
