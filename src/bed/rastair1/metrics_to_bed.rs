use crate::{
    bed::rastair1::{BedRecordsConvertParams, Rastair1BedFormat},
    call::variant_calling::{EstimatedGenotype, GenotypeTag},
    metrics::{DenovoAdjecent, FormsDenovo, PileupMetrics},
    utils::{Base::*, logging::ThisIsABug as _},
    vcf::InCpG,
};
use color_eyre::{Result, Section as _, SectionExt as _, eyre::eyre};
use rastair_types::{Phred, Probability, Strand};
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
        if de_novo {
            if *pileup.pos_metrics.denovo_adj && !pileup.pos_filters.other_pos_in_denovo_passes {
                // This is an adjacent de-novo CpG site, but the other position did not pass filters
                return Ok(None);
            } else if let Some(alt) = pileup.alts.iter().find(|a| *a.metrics.denovo)
                && !alt.filters.pass(Some(ml_threshold))
            {
                // Only report de novo candidates that we are confident about
                trace!(%alt.base, "de novo candidate with low score");
                return Ok(None);
            }
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
            Counts::from_c(pileup)
        } else if cpg && ref_base == G {
            Counts::from_g(pileup)
        } else if let Some(denovo) = pileup.alts.iter().find(|a| *a.metrics.denovo) {
            if denovo.metrics.denovo == FormsDenovo::ThisBecomesC {
                Counts::from_c(pileup)
            } else if denovo.metrics.denovo == FormsDenovo::ThisBecomesG {
                Counts::from_g(pileup)
            } else {
                // Should be never happen since we filtered above
                Counts::default()
            }
        } else if pileup.pos_metrics.denovo_adj == DenovoAdjecent::ThisIsTheMatchingC {
            Counts::from_c(pileup)
        } else if pileup.pos_metrics.denovo_adj == DenovoAdjecent::ThisIsTheMatchingG {
            Counts::from_g(pileup)
        } else {
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
                Some(beta)
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
            Some(beta)
        } else {
            // - No beta available
            warn!(in_cpg=?pileup.pos_metrics.cpg, ?de_novo, genotype=?pileup.pos_metrics.genotype, "why no beta?");
            None
        };

        let strand = guess_strand_from_pileup(pileup);

        let bed = super::Rastair1BedFormat {
            contig: pileup.contig(),
            pos: pileup.pos() as usize,
            r#ref: ref_base.into(),
            beta,
            strand,
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
    fn from_c(pileup: &PileupMetrics) -> Self {
        let c = pileup.allele(C).map(|a| a.strand_count).unwrap_or_default();
        let t = pileup.allele(T).map(|a| a.strand_count).unwrap_or_default();

        Counts { unmod: c.ot, r#mod: t.ot, no_snp: c.ob, snp: t.ob }
    }

    fn from_g(pileup: &PileupMetrics) -> Self {
        let g = pileup.allele(G).map(|a| a.strand_count).unwrap_or_default();
        let a = pileup.allele(A).map(|a| a.strand_count).unwrap_or_default();

        Counts { unmod: g.ob, r#mod: a.ob, no_snp: g.ot, snp: a.ot }
    }

    //     todo: choose whether to use this
    //     fn total(&self) -> u32 {
    //         self.unmod + self.r#mod + self.no_snp + self.snp
    //     }
}

fn guess_strand_from_pileup(pileup: &PileupMetrics) -> Strand {
    if pileup.pos_metrics.cpg == InCpG::C {
        Strand::OT
    } else if pileup.pos_metrics.cpg == InCpG::G {
        Strand::OB
    } else if pileup.pos_metrics.denovo_adj == DenovoAdjecent::ThisIsTheMatchingC {
        Strand::OT
    } else if pileup.pos_metrics.denovo_adj == DenovoAdjecent::ThisIsTheMatchingG {
        Strand::OB
    } else if let Some(denovo) = pileup.alts.iter().filter_map(|a| a.metrics.denovo.some()).next() {
        if denovo == FormsDenovo::ThisBecomesC {
            Strand::OT
        } else if denovo == FormsDenovo::ThisBecomesG {
            Strand::OB
        } else {
            // should never happen since .some above filters these out
            Strand::Unknown
        }
    } else {
        Strand::Unknown
    }
}
