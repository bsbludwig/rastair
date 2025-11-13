use crate::{
    bed::rastair1::{BedRecordsConvertParams, Rastair1BedFormat},
    call::variant_calling::{EstimatedGenotype, GenotypeTag},
    metrics::PileupMetrics,
    utils::{Base::*, logging::ThisIsABug as _},
};
use color_eyre::{
    Result, Section as _, SectionExt as _,
    eyre::{Context as _, eyre},
};
use rastair_types::{Phred, Probability};
use tracing::{debug, instrument, trace, warn};

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
        let ml_threshold = params.ml_threshold;

        let ref_base = pileup.ref_base();

        // If a position is covered by both a ref CpG site and a de-novo CpG
        // site, the ref case should take precedence.
        let cpg = *pileup.pos_metrics.cpg;
        let de_novo = !cpg && pileup.forms_denovo();
        if de_novo
            && !*pileup.pos_metrics.denovo_adj
            && let Some(alt) = pileup.alts.iter().find(|a| *a.metrics.denovo)
            && !alt.filters.pass(Some(ml_threshold))
        {
            // Only report de novo candidates that we are confident about
            trace!(%alt.base, "de novo candidate with low score");
            return Ok(None);
        }

        let gt = if let Some(gt) = pileup.pos_metrics.genotype {
            gt
        } else {
            debug!("No genotype for record");
            EstimatedGenotype {
                genotype: GenotypeTag::CC,
                likelihood: Probability::ONE,
                confidence: Probability::ONE,
            }
        };

        let counts = if cpg && ref_base == C {
            let r = pileup.ref_metrics.strand_count;
            let alt = pileup.alt(T).map(|a| a.strand_count).unwrap_or_default();
            Counts { unmod: r.ot, r#mod: alt.ot, no_snp: r.ob, snp: alt.ob }
        } else if cpg && ref_base == G {
            let r = pileup.ref_metrics.strand_count;
            let alt = pileup.alt(A).map(|a| a.strand_count).unwrap_or_default();
            Counts { unmod: r.ob, r#mod: alt.ob, no_snp: r.ot, snp: alt.ot }
        } else {
            // TOOD: Writing BED but ref is neither C nor G, so this is a de-novo candidate? Handle it!
            Counts::default()
        };

        // If this looks like SNP, set beta to 0
        let beta = if let Some(alt_base) = pileup.pos_metrics.cpg.alt_base()
            && let Some(alt) = pileup.alt_filters(alt_base)
            && let Some(score) = alt.ml
        {
            if score >= ml_threshold {
                // - Does ML say this is a true variant?
                Some(Probability::ZERO)
            } else if let Some(beta) = pileup.pos_metrics.methylated.beta() {
                // - ML says not a variant, use beta if available
                Some(Probability::new(beta).wrap_err("Beta value out of range").this_is_a_bug()?)
            } else {
                // - ML says not a variant, but no beta available
                debug!(%score, %ml_threshold, "ML says not a variant, but no beta available");
                Some(Probability::ZERO)
            }
        } else if *pileup.pos_metrics.cpg
            && let Some(gt) = pileup.pos_metrics.genotype
            && gt.genotype == GenotypeTag::CT
        {
            // - Is it in a CpG and called as heterozygous?
            Some(Probability::ZERO)
        } else if let Some(beta) = pileup.pos_metrics.methylated.beta() {
            // - Just use beta if available
            Some(Probability::new(beta).wrap_err("Beta value out of range").this_is_a_bug()?)
        } else {
            // - No beta available
            warn!(in_cpg=?pileup.pos_metrics.cpg, ?de_novo, genotype=?pileup.pos_metrics.genotype, "why no beta?");
            None
        };

        let bed = super::Rastair1BedFormat {
            contig: pileup.contig(),
            pos: pileup.pos() as usize,
            r#ref: ref_base.into(),
            beta,
            unmod: counts.unmod,
            r#mod: counts.r#mod,
            no_snp: counts.no_snp,
            snp: counts.snp,
            // coverage: counts.total() as usize, // TODO: is coverage meant to include other alts?
            coverage: pileup.pileup.reads.len(),
            genotype: gt.genotype.into(),
            genotype_likelihood: Phred::from(gt.likelihood),
            genotype_confidence: Phred::from(gt.confidence),
            de_novo,
        };

        if cfg!(debug_assertions)
            && let Some(err) = bed.sanity_check()
        {
            Err(eyre!("invalid bed record")).section(err.header("BED errors")).this_is_a_bug()?;
        }

        Ok(Some(bed))
    }
}

#[derive(Debug, Default)]
struct Counts {
    unmod: u32,
    r#mod: u32,
    no_snp: u32,
    snp: u32,
}

impl Counts {
    fn total(&self) -> u32 {
        self.unmod + self.r#mod + self.no_snp + self.snp
    }
}
