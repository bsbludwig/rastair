use crate::{
    bed::rastair1::{BedRecordsConvertParams, Rastair1BedFormat, format::GenotypeString},
    call::variant_calling::GenotypeTag,
    metrics::MethylationEvidenceStrandInfo,
    utils::logging::ThisIsABug,
    vcf::{DeNovoCpGCandidate, GenotypeConfidence, GenotypeLikelihood, InCpG, Methylated},
};
use color_eyre::{
    Result, Section as _, SectionExt as _,
    eyre::{Context as _, ContextCompat as _, ensure, eyre},
};
use rastair_types::SmolStr;
use rastair_types::{Base, SmallVec};
use rastair_types::{Phred, Probability};
use rastair_vcf::{
    VcfField as _,
    standard_fields::{PASS, ReadDepth},
};
use rust_htslib::bcf::Record as HtslibRecord;
use tracing::{instrument, trace};

impl Rastair1BedFormat {
    #[allow(clippy::cast_possible_truncation, reason = "htslib likes i64")]
    #[instrument(level = "trace", skip_all, fields(pos=%r.pos()))]
    pub fn from_vcf(r: &HtslibRecord, params: &BedRecordsConvertParams) -> Result<Option<Self>> {
        let in_cpg = r.info(InCpG::ID.as_bytes()).flag().unwrap_or(false);
        let de_novo = r.info(DeNovoCpGCandidate::ID.as_bytes()).flag().unwrap_or(false);
        let is_pass = r.has_filter(&PASS);

        let relevant = in_cpg || (de_novo && is_pass);
        if !relevant {
            return Ok(None);
        }

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

        let count = MethylationEvidenceStrandInfo::from_vcf(r)
            .wrap_err("Failed to read methylation evidence strand info")?;

        let genotype_alleles: SmallVec<_, 2> = if let Ok(gs) = r.genotypes() {
            // No more genotype remapping needed since VCF no longer mixes '.' with real variants
            gs.get(0).iter().map(|x| (*x).into()).collect()
        } else {
            SmallVec::new()
        };
        let genotype_tag = GenotypeTag::try_from(&genotype_alleles[..]).ok();
        let alt_bases: Vec<Base> = alleles.iter().skip(1).map(Base::from).collect();
        let genotype = if let Some(gt) = genotype_tag {
            GenotypeString::from_genotype_tag(gt, Base::from(&r#ref), &alt_bases)
        } else {
            // Fallback for invalid genotype
            GenotypeString(Base::from(&r#ref), Base::from(&r#ref))
        };
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

        // If this is a CpG position with a called CpG-relevant variant, set beta to 0.
        // We check both:
        // 1. The CpG-relevant alt allele (T for C, A for G) exists with PASS filter
        // 2. The genotype actually calls the variant (not hom ref 0/0)
        // Without the genotype check, methylation evidence (which shows as T/A reads
        // at CpG sites) would incorrectly set beta=0 when using --no-ml.
        let is_called_variant = genotype_tag.is_some_and(|gt| !gt.is_hom_ref());
        let has_cpg_snp = in_cpg && is_pass && is_called_variant && {
            // Determine which alt base would be the CpG-relevant SNP (T for C, A for G)
            let cpg_snp_base = if r#ref == "C" { "T" } else { "A" };
            // Check if this alt exists in the alleles list
            alleles.iter().skip(1).any(|a| a.as_str() == cpg_snp_base)
        };

        // Set beta to 0 if there's a called CpG-relevant SNP, otherwise use beta from VCF
        let beta = if has_cpg_snp { Some(0.0) } else { beta };
        let beta = if let Some(beta) = beta {
            Some(Probability::new(beta).wrap_err("Beta value out of range").this_is_a_bug()?)
        } else {
            trace!(pos=%contig, pos=r.pos(), ?in_cpg, ?genotype, "why no beta?");
            Some(Probability::ZERO)
        };

        let strand = if in_cpg {
            if r#ref == "C" { rastair_types::Strand::OT } else { rastair_types::Strand::OB }
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
                rastair_types::Strand::OT
            } else if has_alt_g {
                // This position has a variant creating a G
                rastair_types::Strand::OB
            } else if r#ref == "G" {
                // Adjacent G position (partner to a C variant)
                rastair_types::Strand::OB
            } else if r#ref == "C" {
                // Adjacent C position (partner to a G variant)
                rastair_types::Strand::OT
            } else {
                rastair_types::Strand::Unknown
            }
        } else {
            rastair_types::Strand::Unknown
        };

        let bed = Rastair1BedFormat {
            contig: contig.clone(),
            pos: r.pos() as usize,
            r#ref,
            beta,
            strand,
            unmod: count.unmod,
            r#mod: count.modified,
            no_snp: count.no_snp,
            snp: count.snp,
            coverage: read_depth as usize,
            genotype,
            genotype_likelihood,
            genotype_confidence,
            de_novo: !in_cpg && de_novo,
        };

        if cfg!(debug_assertions)
            && let Some(err) = bed.sanity_check()
        {
            Err(eyre!("invalid bed record"))
                .section(err.header("BED errors"))
                .with_note(|| format!("Position {contig}:{}", r.pos()))
                .with_note(|| format!("CPG={in_cpg}, CPGnovo={de_novo}"))
                .this_is_a_bug()?;
        }

        Ok(Some(bed))
    }
}

impl MethylationEvidenceStrandInfo {
    #[instrument(level = "trace", skip_all)]
    fn from_vcf(r: &HtslibRecord) -> Result<Self> {
        let nums = r
            .info(MethylationEvidenceStrandInfo::ID.as_bytes())
            .integer()
            .wrap_err("Failed to fetch field")?
            .wrap_err("field not set")?;
        ensure!(nums.len() == 4, "field has invalid length");

        #[expect(clippy::get_first, reason = "consistency")]
        Ok(MethylationEvidenceStrandInfo {
            unmod: nums.get(0).copied().wrap_err("missing unmod count")? as u32,
            modified: nums.get(1).copied().wrap_err("missing modified count")? as u32,
            no_snp: nums.get(2).copied().wrap_err("missing no_snp count")? as u32,
            snp: nums.get(3).copied().wrap_err("missing snp count")? as u32,
        })
    }
}
