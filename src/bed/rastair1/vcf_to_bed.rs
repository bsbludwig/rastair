use crate::{
    bed::rastair1::Rastair1BedFormat,
    utils::Base,
    vcf::{
        AlleleSpecificStrandBias, ByStrand, DeNovoCpGCandidate, GenotypeConfidence,
        GenotypeLikelihood, InCpG, Methylated,
    },
};
use color_eyre::{
    Result,
    eyre::{Context as _, ContextCompat as _, Report, ensure},
};
use rastair2_types::Phred;
use rastair2_vcf::{
    StrandSpecificInfoField as _, VcfField as _,
    standard_fields::{Genotype, ReadDepth},
};
use rust_htslib::bcf::Record as HtslibRecord;
use smallvec::SmallVec;
use smol_str::SmolStr;

impl TryFrom<&HtslibRecord> for Rastair1BedFormat {
    type Error = Report;

    #[allow(clippy::cast_possible_truncation)]
    fn try_from(r: &HtslibRecord) -> Result<Self, Self::Error> {
        let contig = r
            .rid()
            .wrap_err("Record has no ID")
            .and_then(|id| r.header().rid2name(id).wrap_err("Header does not contain ID"))
            .and_then(|name| str::from_utf8(name).wrap_err("Contig name is not valid UTF-8"))
            .map(SmolStr::new)
            .wrap_err("Could not fetch contig name")?;

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
            *beta
        } else {
            0.0
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

        let assb = AlleleSpecificStrandBias::from_vcf(r)
            .wrap_err("Failed to read Allele-specific Strand Bias field")?;

        let count =
            StrandCount::from_assb(&assb).wrap_err("Failed to parse strand counts from record")?;
        let (unmod, r#mod, no_snp, snp) = if r#ref == "C" {
            (count.c.ot, count.t.ot, count.c.ob, count.t.ob)
        } else if r#ref == "G" {
            (count.g.ob, count.a.ob, count.g.ot, count.a.ot)
        } else {
            (0, 0, 0, 0)
        };

        let genotype = if let Ok(gs) = r.genotypes() {
            gs.get(0).iter().map(|x| (*x).into()).collect()
        } else {
            SmallVec::new()
        };
        let genotype_likelihood = GenotypeLikelihood(SmallVec::from_buf([{
            if let Ok(buffer) = r.format(GenotypeLikelihood::ID.as_bytes()).integer()
                && let Some(first) = buffer.first()
                && let Some(val) = first.first()
            {
                Some(Phred::from_phred(*val))
            } else {
                None
            }
        }]));
        let genotype_confidence = GenotypeConfidence(SmallVec::from_buf([{
            if let Ok(buffer) = r.format(GenotypeConfidence::ID.as_bytes()).integer()
                && let Some(first) = buffer.first()
                && let Some(val) = first.first()
            {
                Some(Phred::from_phred(*val))
            } else {
                None
            }
        }]));
        let in_cpg = r.info(InCpG::ID.as_bytes()).flag().unwrap_or(false);
        let de_novo = r.info(DeNovoCpGCandidate::ID.as_bytes()).flag().unwrap_or(false);

        Ok(Rastair1BedFormat {
            contig,
            pos: r.pos() as usize,
            r#ref,
            beta,
            unmod,
            r#mod,
            no_snp,
            snp,
            coverage: read_depth as usize,
            genotype: Genotype(genotype),
            genotype_likelihood,
            genotype_confidence,
            de_novo: !in_cpg && de_novo,
        })
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
