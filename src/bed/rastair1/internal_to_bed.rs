use crate::{
    bed::rastair1::Rastair1BedFormat,
    utils::Base::*,
    vcf::{Record as Rastair2Record, utils::NoStrandBiasForBaseErrorExt as _},
};
use color_eyre::eyre::Report;

impl TryFrom<&Rastair2Record> for Rastair1BedFormat {
    type Error = Report;

    #[allow(clippy::cast_possible_truncation)]
    fn try_from(record: &Rastair2Record) -> Result<Self, Self::Error> {
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

        Ok(super::Rastair1BedFormat {
            contig: record.main.chrom.clone(),
            pos: record.main.pos as usize,
            r#ref,
            beta: record.samples[0].methylated.beta().unwrap_or_default() as f32,
            unmod,
            r#mod,
            no_snp,
            snp,
            coverage: *record.info.read_depth,
            genotype: record.samples[0].genotype.clone(),
            genotype_likelihood: record.samples[0].genotype_likelihood.clone(),
            genotype_confidence: record.samples[0].genotype_confidence.clone(),
            de_novo: !*record.info.in_cp_g && *record.info.de_novo_cp_g_candidate,
        })
    }
}
