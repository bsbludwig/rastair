use crate::{
    bed::rastair1::{BedRecordsFilterParams, Rastair1BedFormat, format::GenotypeString},
    call::variant_calling::GenotypeTag,
    metrics::MethylationEvidenceStrandInfo,
    utils::logging::ThisIsABug,
};
use color_eyre::{
    Result, Section as _, SectionExt as _,
    eyre::{Context as _, ContextCompat as _, eyre},
};
use seqair::vcf::{
    FilterStatus, FormatFloatOpts, FormatGt, FormatInt, Info, InfoFlag, InfoInt, InfoInts, Keys,
    Locus, OnlySample, Site,
};
use seqair_types::{Base, SmallVec, SmolStr};
use seqair_types::{Phred, Probability};
use tracing::{instrument, trace};

/// What `convert` reads of a VCF record to make its BED row.
pub type VcfRecord = (Locus, Site, Info, OnlySample);

/// The fields a BED row is made from, where the header declares them.
#[derive(Debug)]
pub struct VcfKeys {
    cpg: Option<InfoFlag>,
    cpg_novo: Option<InfoFlag>,
    depth: Option<InfoInt>,
    strands: Option<InfoInts>,
    genotype: Option<FormatGt>,
    genotype_quality: Option<FormatInt>,
    beta: Option<FormatFloatOpts>,
}

impl VcfKeys {
    pub fn of(keys: &Keys<'_>) -> Result<Self> {
        Ok(Self {
            cpg: keys.info_key_opt("CPG")?,
            cpg_novo: keys.info_key_opt("CPGnovo")?,
            depth: keys.info_key_opt("DP")?,
            strands: keys.info_key_opt("M5mC_Strands")?,
            genotype: keys.format_key_opt("GT")?,
            genotype_quality: keys.format_key_opt("GQ")?,
            beta: keys.format_key_opt("M5mC")?,
        })
    }
}

impl Rastair1BedFormat {
    #[instrument(level = "trace", skip_all, fields(at = %locus))]
    pub fn from_vcf(
        keys: &VcfKeys,
        (locus, site, info, sample): &VcfRecord,
        params: &BedRecordsFilterParams,
    ) -> Result<Option<Self>> {
        let flag = |key: &Option<InfoFlag>| -> Result<bool> {
            Ok(match key {
                Some(key) => info.get(key)?,
                None => false,
            })
        };
        let in_cpg = flag(&keys.cpg)?;
        let de_novo = flag(&keys.cpg_novo)?;
        let is_pass = matches!(site.filters().status(), FilterStatus::Pass);

        let relevant = in_cpg || (de_novo && is_pass);
        if !relevant {
            return Ok(None);
        }

        let contig = locus.contig().name().clone();
        let pos = locus.start().pos0().wrap_err("Record before the first base")?.as_usize();

        let alleles = site
            .alleles()
            .iter()
            .map(|allele| str::from_utf8(allele.as_bytes()).map(SmolStr::new))
            .collect::<Result<SmallVec<_, 4>, _>>()
            .wrap_err("Failed to parse alleles")?;
        let r#ref = alleles.first().wrap_err("Record has no reference allele")?.clone();

        // `call --bed`'s β: the reference CpG's, else the de-novo one's,
        // whichever a read shows first.
        let beta = match &keys.beta {
            Some(key) => sample.get(key)?.and_then(|betas| betas.flatten().next()).map(f64::from),
            None => None,
        };

        let read_depth = match &keys.depth {
            Some(key) => info.get(key).wrap_err("Could not fetch read depth from record")?,
            None => None,
        }
        .unwrap_or(0);

        // Skip positions without evidence
        if !params.include_empty && read_depth == 0 {
            return Ok(None);
        }

        let count = MethylationEvidenceStrandInfo::from_vcf(keys, info)
            .wrap_err("Failed to read methylation evidence strand info")?;

        // Taken as written: a record never mixes '.' with real alleles.
        let genotype_tag = match &keys.genotype {
            Some(key) => sample.get(key)?.and_then(|genotype| {
                let mut numbers = genotype.alleles().map(|allele| allele.allele);
                match (numbers.next(), numbers.next(), numbers.next()) {
                    (Some(Some(first)), Some(Some(second)), None) => {
                        GenotypeTag::from_allele_numbers([first, second]).ok()
                    }
                    _ => None,
                }
            }),
            None => None,
        };
        let alt_bases: Vec<Base> = alleles.iter().skip(1).map(Base::from).collect();
        let genotype = if let Some(gt) = genotype_tag {
            GenotypeString::from_genotype_tag(gt, Base::from(&r#ref), &alt_bases)
        } else {
            // Fallback for invalid genotype
            GenotypeString(Base::from(&r#ref), Base::from(&r#ref))
        };
        let genotype_confidence = match &keys.genotype_quality {
            Some(key) => sample.get(key)?,
            None => None,
        }
        .map(|quality| Phred::from_phred(u8::try_from(quality.clamp(0, 255)).unwrap_or(0)));
        if genotype_confidence.is_none() {
            trace!(?genotype, "No genotype confidence field found");
        }
        let genotype_confidence = super::format::gq_column(genotype_confidence);

        // `M5mC` is the call's β: 0 where both chromosomes lose the CpG, the
        // remaining allele's at a heterozygous SNV, missing where no read shows it.
        let beta = beta
            .map(|beta| Probability::new(beta).wrap_err("Beta value out of range").this_is_a_bug())
            .transpose()?;

        let strand = if in_cpg {
            if r#ref == "C" {
                Some(seqair_types::Strand::OT)
            } else {
                Some(seqair_types::Strand::OB)
            }
        } else if de_novo {
            // Infer strand from de-novo CpG information:
            // De-novo CpGs are created when:
            // - Any base → C followed by G creates CG (C is methylation site, OT strand)
            // - C followed by any base → G creates CG (G is on OB strand)
            //
            // For a position with CPGnovo flag:
            // - If alt contains C → this position becomes C (OT strand)
            // - If alt contains G → this position becomes G (OB strand)
            // - If ref=G and no alt → adjacent G to a C variant (OB strand)
            // - If ref=C and no alt → adjacent C to a G variant (OT strand)
            let has_alt_c = alleles.iter().skip(1).any(|a| a.as_str() == "C");
            let has_alt_g = alleles.iter().skip(1).any(|a| a.as_str() == "G");

            if has_alt_c {
                // This position has a variant creating a C
                Some(seqair_types::Strand::OT)
            } else if has_alt_g {
                // This position has a variant creating a G
                Some(seqair_types::Strand::OB)
            } else if r#ref == "G" {
                // Adjacent G position (partner to a C variant)
                Some(seqair_types::Strand::OB)
            } else if r#ref == "C" {
                // Adjacent C position (partner to a G variant)
                Some(seqair_types::Strand::OT)
            } else {
                None
            }
        } else {
            None
        };

        let bed = Rastair1BedFormat {
            contig: contig.clone(),
            pos,
            r#ref,
            beta,
            strand,
            unmod: count.unmod,
            r#mod: count.modified,
            no_snp: count.no_snp,
            snp: count.snp,
            coverage: read_depth as usize,
            genotype,
            genotype_likelihood: genotype_confidence,
            genotype_confidence,
            de_novo: !in_cpg && de_novo,
        };

        if cfg!(debug_assertions)
            && let Some(err) = bed.sanity_check()
        {
            Err(eyre!("invalid bed record"))
                .section(err.header("BED errors"))
                .with_note(|| format!("Position {locus}"))
                .with_note(|| format!("CPG={in_cpg}, CPGnovo={de_novo}"))
                .this_is_a_bug()?;
        }

        Ok(Some(bed))
    }
}

impl MethylationEvidenceStrandInfo {
    #[instrument(level = "trace", skip_all)]
    fn from_vcf(keys: &VcfKeys, info: &Info) -> Result<Self> {
        let key = keys.strands.as_ref().wrap_err("the header does not declare `M5mC_Strands`")?;
        let nums = info.get(key).wrap_err("Failed to fetch field")?.wrap_err("field not set")?;
        let count = |value: i32| u32::try_from(value).wrap_err("negative read count");
        let [unmod, modified, no_snp, snp] =
            nums.to_array::<4>().wrap_err("field has invalid length")?;
        Ok(MethylationEvidenceStrandInfo {
            unmod: count(unmod)?,
            modified: count(modified)?,
            no_snp: count(no_snp)?,
            snp: count(snp)?,
        })
    }
}
