use crate::{
    bed::rastair1::{BedRecordsConvertParams, Rastair1BedFormat},
    call::variant_calling::GenotypeTag,
    metrics::PileupMetrics,
    utils::{Base::*, ByStrand, logging::ThisIsABug as _},
    vcf::{GenotypeConfidence, GenotypeLikelihood},
};
use color_eyre::{Result, eyre::Context as _};
use rastair_types::{Phred, Probability};
use smallvec::smallvec_inline;
use tracing::{instrument, trace, warn};

impl Rastair1BedFormat {
    #[instrument(level = "trace", skip_all, fields(pos=%pileup.contig_pos()))]
    pub fn from_metrics(
        pileup: &PileupMetrics,
        params: &BedRecordsConvertParams,
    ) -> Result<Option<Self>> {
        if !params.filters.include_empty && pileup.pos_metrics.depth == 0 {
            trace!("no coverage, skipping");
            return Ok(None);
        }
        if !(*pileup.pos_metrics.cpg || pileup.forms_denovo()) {
            trace!("in neither ref CpG nor de-novo CpG, skipping");
            return Ok(None);
        }

        // If a position is covered by both a ref CpG site and a de-novo CpG
        // site, the ref case should take precedence.
        let de_novo = !*pileup.pos_metrics.cpg && pileup.forms_denovo();

        let Some(gt) = pileup.pos_metrics.genotype else {
            warn!("No genotype for record");
            return Ok(None);
        };

        #[derive(Debug, Default)]
        struct Counts {
            unmod: u32,
            r#mod: u32,
            no_snp: u32,
            snp: u32,
        }

        impl Counts {
            fn new(r: ByStrand<u32>, alt: ByStrand<u32>) -> Self {
                Counts { unmod: r.ot, r#mod: alt.ot, no_snp: r.ob, snp: alt.ob }
            }

            fn total(&self) -> u32 {
                self.unmod + self.r#mod + self.no_snp + self.snp
            }
        }

        let counts = if pileup.ref_base() == C {
            let r = pileup.ref_metrics.strand_count;
            let alt = pileup.alt(T).map(|a| a.strand_count).unwrap_or_default();

            Counts::new(r, alt)
        } else if pileup.ref_base() == G {
            let r = pileup.ref_metrics.strand_count;
            let alt = pileup.alt(A).map(|a| a.strand_count).unwrap_or_default();

            Counts::new(r, alt)
        } else {
            // Counts::default()
            todo!("Writing BED but ref is neither C nor G, so this is a de-novo candidate?")
        };

        // If this looks like SNP, set beta to 0
        // - Does ML say this is a true variant?
        let beta = if let Some(alt_base) = pileup.pos_metrics.cpg.alt_base()
            && let Some(alt) = pileup.alt_filters(alt_base)
            && let Some(score) = alt.ml
            && let threshold = params.ml_threshold
            && score >= threshold
        {
            Some(Probability::new_panicky(0.0))
        } else if *pileup.pos_metrics.cpg
            && let Some(gt) = pileup.pos_metrics.genotype
            && gt.genotype == GenotypeTag::CT
        {
            Some(Probability::new_panicky(0.0))
        } else if let Some(beta) = pileup.pos_metrics.methylated.beta() {
            Some(Probability::new(beta).wrap_err("Beta value out of range").this_is_a_bug()?)
        } else {
            warn!(in_cpg=?pileup.pos_metrics.cpg, genotype=?pileup.pos_metrics.genotype, "why no beta?");
            None
        };

        let bed = super::Rastair1BedFormat {
            contig: pileup.contig(),
            pos: pileup.pos() as usize,
            r#ref: pileup.ref_base().into(),
            beta,
            unmod: counts.unmod,
            r#mod: counts.r#mod,
            no_snp: counts.no_snp,
            snp: counts.snp,
            coverage: counts.total() as usize,
            genotype: gt.genotype.into(),
            genotype_likelihood: GenotypeLikelihood(smallvec_inline![Some(
                Phred::from_probability(gt.likelihood)
                    .wrap_err("Genotype likelihood out of range")
                    .this_is_a_bug()?
            )]),
            genotype_confidence: GenotypeConfidence(smallvec_inline![Some(
                Phred::from_probability(gt.confidence)
                    .wrap_err("Genotype likelihood out of range")
                    .this_is_a_bug()?
            )]),
            de_novo,
        };

        if cfg!(debug_assertions) {
            bed.sanity_check().wrap_err("BED record failed sanity check")?;
        }

        Ok(Some(bed))
    }
}
