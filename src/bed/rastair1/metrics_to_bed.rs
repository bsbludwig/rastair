use crate::{
    bed::rastair1::{Rastair1BedFormat, format::GenotypeString},
    call::variant_calling::GenotypeTag,
    metrics::{MethylationEvidenceStrandInfo, PileupMetrics, methylation},
    utils::logging::ThisIsABug as _,
};
use color_eyre::{Result, Section as _, SectionExt as _, eyre::eyre};
use tracing::{debug, instrument, trace, warn};

impl Rastair1BedFormat {
    #[instrument(level = "trace", skip_all, fields(pos=%pileup.contig_pos()))]
    pub fn from_metrics(pileup: &PileupMetrics) -> Result<Option<Self>> {
        if !pileup.cpg_record() {
            trace!("in neither ref CpG nor de-novo CpG, skipping");
            return Ok(None);
        }

        let ref_base = pileup.ref_base();

        // If a position is covered by both a ref CpG site and a de-novo CpG
        // site, the ref case should take precedence.
        let de_novo = pileup.reference_cpg().is_none() && pileup.in_denovo_cpg();

        let (genotype, quality) = if let Some(gt) = &pileup.pos_metrics.genotype {
            (gt.genotype, Some(gt.quality))
        } else {
            debug!("No genotype for record");
            (GenotypeTag::HomRef, None)
        };
        let quality = super::format::gq_column(quality);

        let counts = MethylationEvidenceStrandInfo::from_pileup(pileup);

        // β, strand and counts are all the evidence CpG's.
        let evidence = methylation::evidence(pileup);
        let beta = evidence
            .and_then(|(_, origin)| pileup.pos_metrics.methylated.get(origin))
            .map(|cpg| cpg.beta);
        let strand = evidence.map(|(side, _)| side.strand());

        let alt_bases: Vec<_> = pileup.alts.iter().map(|a| a.base).collect();

        let bed = super::Rastair1BedFormat {
            contig: pileup.region.contig.clone(),
            pos: pileup.pos().as_usize(),
            r#ref: ref_base.into(),
            beta,
            strand,
            unmod: counts.unmod,
            r#mod: counts.modified,
            no_snp: counts.no_snp,
            snp: counts.snp,
            // Coverage is total read depth at this position, including all bases/alts
            coverage: pileup.pos_metrics.depth as usize,
            genotype: GenotypeString::from_genotype_tag(genotype, ref_base, &alt_bases),
            // The VCF has one genotype number, `GQ`; the legacy column repeats it.
            genotype_likelihood: quality,
            genotype_confidence: quality,
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
