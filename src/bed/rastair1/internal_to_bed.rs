use crate::{
    bed::rastair1::Rastair1BedFormat,
    utils::Base::*,
    vcf::{Record as Rastair2Record, utils::NoStrandBiasForBaseErrorExt as _},
};
use color_eyre::Result;
use tracing::trace;

impl Rastair1BedFormat {
    #[allow(clippy::cast_possible_truncation)]
    pub fn from_record(record: &Rastair2Record) -> Result<Option<Self>> {
        if *record.info.de_novo_cp_g_candidate && !record.filters.pass() {
            // Only report de novo candidates that we are confident about
            trace!(
                chr=%record.main.chrom, pos=record.main.pos,
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
        let beta = if *record.info.in_cp_g && record.samples[0].genotype.homozygous_not_ref() {
            Some(0.0)
        } else {
            record.samples[0].methylated.beta()
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
            de_novo: !*record.info.in_cp_g && *record.info.de_novo_cp_g_candidate,
        }))
    }
}
