use crate::{
    bed::rastair1::{BedRecordsConvertParams, Rastair1BedFormat},
    utils::{Base, ByStrand, logging::ThisIsABug},
    vcf::{
        AlleleSpecificStrandBias, DeNovoCpGCandidate, GenotypeConfidence, GenotypeLikelihood,
        InCpG, Methylated,
    },
};
use color_eyre::{
    Result, Section as _, SectionExt as _,
    eyre::{Context as _, ContextCompat as _, Report, ensure, eyre},
};
use rastair_types::{Phred, Probability};
use rastair_vcf::{
    StrandSpecificInfoField as _, VcfField as _,
    standard_fields::{Genotype, PASS, ReadDepth},
};
use rust_htslib::bcf::Record as HtslibRecord;
use smallvec::SmallVec;
use smol_str::SmolStr;
use tracing::{instrument, trace};

impl Rastair1BedFormat {
    #[allow(clippy::cast_possible_truncation)]
    #[instrument(level = "trace", skip_all, fields(pos=%r.pos()))]
    pub fn from_vcf(r: &HtslibRecord, params: &BedRecordsConvertParams) -> Result<Option<Self>> {
        let contig = r
            .rid()
            .wrap_err("Record has no ID")
            .and_then(|id| r.header().rid2name(id).wrap_err("Header does not contain ID"))
            .and_then(|name| str::from_utf8(name).wrap_err("Contig name is not valid UTF-8"))
            .map(SmolStr::new)
            .wrap_err("Could not fetch contig name")?;

        if tracing::enabled!(tracing::Level::TRACE) {
            tracing::Span::current().record("contig", tracing::field::display(&contig));
        }

        let alleles = r
            .alleles()
            .iter()
            .map(|a| str::from_utf8(a).map(SmolStr::new))
            .collect::<Result<SmallVec<_, 4>, _>>()
            .wrap_err("Failed to parse alleles")?;
        let r#ref = alleles.first().wrap_err("Record has no reference allele")?.clone();

        let beta = if let Ok(buffer) = r.format(Methylated::ID.as_bytes()).float()
            && let Some(betas) = buffer.first()
            && let Some(beta) = betas.first()
        {
            Some(f64::from(*beta))
        } else {
            None
        };

        let read_depth = if let Some(buffer) = r
            .info(ReadDepth::ID.as_bytes())
            .integer()
            .wrap_err("Could not fetch read depth from record")?
            && let Some(depth) = buffer.first()
        {
            *depth
        } else {
            0
        };

        // Skip positions without evidence
        if !params.filters.include_empty && read_depth == 0 {
            return Ok(None);
        }

        let assb = AlleleSpecificStrandBias::from_vcf(r)
            .wrap_err("Failed to read Allele-specific Strand Bias field")?;

        let count =
            StrandCount::from_assb(&assb).wrap_err("Failed to parse strand counts from record")?;
        let (unmod, r#mod, no_snp, snp) = if r#ref == "C" {
            (count.c.ot, count.t.ot, count.c.ob, count.t.ob)
        } else if r#ref == "G" {
            (count.g.ob, count.a.ob, count.g.ot, count.a.ot)
        } else {
            trace!(
                %r#ref,
                "Writing BED but ref is neither C nor G, so this is a de-novo candidate?"
            );
            (0, 0, 0, 0)
        };

        let genotype = if let Ok(gs) = r.genotypes() {
            gs.get(0).iter().map(|x| (*x).into()).collect()
        } else {
            SmallVec::new()
        };
        let genotype = Genotype(genotype);
        let genotype_likelihood = if let Ok(buffer) =
            r.format(GenotypeLikelihood::ID.as_bytes()).integer()
            && let Some(first) = buffer.first()
            && let Some(val) = first.first()
        {
            Phred::from_phred(*val)
        } else {
            trace!(?genotype, "No genotype likelihood field found");
            Phred::from_phred(0)
        };
        let genotype_confidence = if let Ok(buffer) =
            r.format(GenotypeConfidence::ID.as_bytes()).integer()
            && let Some(first) = buffer.first()
            && let Some(val) = first.first()
        {
            Phred::from_phred(*val)
        } else {
            trace!(?genotype, "No genotype confidence field found");
            Phred::from_phred(0)
        };
        let in_cpg = r.info(InCpG::ID.as_bytes()).flag().unwrap_or(false);
        let de_novo = r.info(DeNovoCpGCandidate::ID.as_bytes()).flag().unwrap_or(false);

        if de_novo && !r.has_filter(&PASS) {
            return Ok(None);
        }

        // TODO: Use ML
        let beta = if in_cpg && genotype.homozygous_not_ref() { Some(0.0) } else { beta };
        let beta = if let Some(beta) = beta {
            Some(Probability::new(beta).wrap_err("Beta value out of range").this_is_a_bug()?)
        } else {
            trace!(pos=%contig, pos=r.pos(), ?in_cpg, ?genotype, "why no beta?");
            None
        };

        let bed = Rastair1BedFormat {
            contig,
            pos: r.pos() as usize,
            r#ref,
            beta,
            unmod,
            r#mod,
            no_snp,
            snp,
            coverage: read_depth as usize,
            genotype,
            genotype_likelihood,
            genotype_confidence,
            de_novo: !in_cpg && de_novo,
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
struct StrandCount {
    a: ByStrand<u32>,
    c: ByStrand<u32>,
    g: ByStrand<u32>,
    t: ByStrand<u32>,
}

impl StrandCount {
    fn from_assb(assb: &AlleleSpecificStrandBias) -> Result<StrandCount, Report> {
        // AS_SB is encoded with two integers per allele, so we need to parse it accordingly
        let mut counts = StrandCount::default();
        for count in assb.iter() {
            match count.base {
                Base::A => counts.a = *count,
                Base::C => counts.c = *count,
                Base::G => counts.g = *count,
                Base::T => counts.t = *count,
                _ => {}
            }
        }

        Ok(counts)
    }
}

impl AlleleSpecificStrandBias {
    fn from_vcf(r: &HtslibRecord) -> Result<Self> {
        let alleles = r
            .alleles()
            .iter()
            .map(|a| str::from_utf8(a).map(SmolStr::new))
            .collect::<Result<SmallVec<_, 4>, _>>()
            .wrap_err("Failed to parse alleles")?;

        // we have two fields: AS_SB_OT and AS_SB_OB
        let assb_ot = r
            .info(AlleleSpecificStrandBias::ID_OT.as_bytes())
            .integer()
            .wrap_err("Failed to fetch AS_SB_OT field")?
            .wrap_err("AS_SB_OT field not set")?;
        let assb_ob = r
            .info(AlleleSpecificStrandBias::ID_OB.as_bytes())
            .integer()
            .wrap_err("Failed to fetch AS_SB_OB field")?
            .wrap_err("AS_SB_OB field not set")?;

        ensure!(
            assb_ot.len() == assb_ob.len(),
            "AS_SB_OT and AS_SB_OB fields have different lengths"
        );

        Ok(AlleleSpecificStrandBias(
            alleles
                .iter()
                .enumerate()
                .map(|(i, allele)| ByStrand {
                    base: Base::from(allele),
                    ot: assb_ot.get(i).copied().unwrap_or(0) as u32,
                    ob: assb_ob.get(i).copied().unwrap_or(0) as u32,
                })
                .collect(),
        ))
    }
}
