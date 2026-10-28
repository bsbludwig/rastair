use crate::{
    bed::rastair1::{BedRecordsFilterParams, Rastair1BedFormat, format::GenotypeString},
    call::variant_calling::GenotypeTag,
    metrics::MethylationEvidenceStrandInfo,
    utils::logging::ThisIsABug,
};
use color_eyre::{
    Result, Section as _, SectionExt as _,
    eyre::{Context as _, ContextCompat as _, ensure, eyre},
};
use rust_htslib::bcf::{Record as HtslibRecord, record::Numeric as _};
use seqair_types::SmolStr;
use seqair_types::{Base, SmallVec};
use seqair_types::{Phred, Probability};
use tracing::{instrument, trace};

impl Rastair1BedFormat {
    #[allow(clippy::cast_possible_truncation, reason = "htslib likes i64")]
    #[instrument(level = "trace", skip_all, fields(pos=%r.pos()))]
    pub fn from_vcf(r: &HtslibRecord, params: &BedRecordsFilterParams) -> Result<Option<Self>> {
        let in_cpg = r.info(b"CPG").flag().unwrap_or(false);
        let de_novo = r.info(b"CPGnovo").flag().unwrap_or(false);
        let is_pass = r.has_filter("PASS".as_bytes());

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

        // `call --bed`'s β: the reference CpG's, else the de-novo one's,
        // whichever a read shows first.
        let beta = if let Ok(buffer) = r.format(b"M5mC").float()
            && let Some(betas) = buffer.first()
            && let Some(beta) = betas.iter().find(|beta| !beta.is_missing())
        {
            Some(f64::from(*beta))
        } else {
            None
        };

        let read_depth = if let Some(buffer) =
            r.info(b"DP").integer().wrap_err("Could not fetch read depth from record")?
            && let Some(depth) = buffer.first()
        {
            *depth
        } else {
            0
        };

        // Skip positions without evidence
        if !params.include_empty && read_depth == 0 {
            return Ok(None);
        }

        let count = MethylationEvidenceStrandInfo::from_vcf(r)
            .wrap_err("Failed to read methylation evidence strand info")?;

        let genotype_alleles: SmallVec<_, 2> = if let Ok(gs) = r.genotypes() {
            // Taken as written: a record never mixes '.' with real alleles.
            gs.get(0).iter().copied().collect()
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
        let genotype_confidence = if let Ok(buffer) = r.format(b"GQ").integer()
            && let Some(first) = buffer.first()
            && let Some(&val) = first.first()
            // htslib's missing integer is `i32::MIN`, which is `.` in the VCF.
            && val != i32::MIN
        {
            Some(Phred::from_phred(u8::try_from(val.clamp(0, 255)).unwrap_or(0)))
        } else {
            trace!(?genotype, "No genotype confidence field found");
            None
        };
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
            genotype_likelihood: genotype_confidence,
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
            .info(b"M5mC_Strands")
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
