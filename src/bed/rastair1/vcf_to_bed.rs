use crate::{
    bed::rastair1::{BedRecordsConvertParams, Rastair1BedFormat, format::GenotypeString},
    call::variant_calling::GenotypeTag,
    metrics::MethylationEvidenceStrandInfo,
    utils::logging::ThisIsABug,
};
use color_eyre::{
    Result, Section as _, SectionExt as _,
    eyre::{Context as _, ContextCompat as _, ensure, eyre},
};
use rust_htslib::bcf::Record as HtslibRecord;
use seqair_types::SmolStr;
use seqair_types::{Base, SmallVec};
use seqair_types::{Phred, Probability};
use tracing::{instrument, trace};

impl Rastair1BedFormat {
    #[allow(clippy::cast_possible_truncation, reason = "htslib likes i64")]
    #[instrument(level = "trace", skip_all, fields(pos=%r.pos()))]
    pub fn from_vcf(r: &HtslibRecord, params: &BedRecordsConvertParams) -> Result<Option<Self>> {
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

        let beta = if let Ok(buffer) = r.format(b"M5mC").float()
            && let Some(betas) = buffer.first()
            && let Some(beta) = betas.first()
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
        if !params.filters.include_empty && read_depth == 0 {
            return Ok(None);
        }

        // The first `DPM5mC`/`ADM5mC` value belongs to the same CpG as the first `M5mC`.
        let total = first_format_count(r, b"DPM5mC");
        let r#mod = first_format_count(r, b"ADM5mC");
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
        let genotype_likelihood = if let Ok(buffer) = r.format(b"GL").float()
            && let Some(first) = buffer.first()
            && let Some(val) = first.first()
        {
            Phred::from_phred((*val).clamp(0.0, 255.0) as u8)
        } else {
            trace!(?genotype, "No genotype likelihood field found");
            Phred::from_phred(0_u8)
        };
        let genotype_confidence = if let Ok(buffer) = r.format(b"GC").float()
            && let Some(first) = buffer.first()
            && let Some(val) = first.first()
        {
            Phred::from_phred((*val).clamp(0.0, 255.0) as u8)
        } else {
            trace!(?genotype, "No genotype confidence field found");
            Phred::from_phred(0_u8)
        };

        let beta = if let Some(beta) = beta {
            Some(Probability::new(beta).wrap_err("Beta value out of range").this_is_a_bug()?)
        } else {
            trace!(pos=%contig, pos=r.pos(), ?in_cpg, ?genotype, "why no beta?");
            Some(Probability::ZERO)
        };

        let strand = if in_cpg {
            if r#ref == "C" { seqair_types::Strand::OT } else { seqair_types::Strand::OB }
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
                seqair_types::Strand::OT
            } else if has_alt_g {
                // This position has a variant creating a G
                seqair_types::Strand::OB
            } else if r#ref == "G" {
                // Adjacent G position (partner to a C variant)
                seqair_types::Strand::OB
            } else if r#ref == "C" {
                // Adjacent C position (partner to a G variant)
                seqair_types::Strand::OT
            } else {
                seqair_types::Strand::Unknown
            }
        } else {
            seqair_types::Strand::Unknown
        };

        let bed = Rastair1BedFormat {
            contig: contig.clone(),
            pos: r.pos() as usize,
            r#ref,
            beta,
            strand,
            unmod: total.saturating_sub(r#mod),
            r#mod,
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

/// The first value of an integer FORMAT field, or 0 where the record does not carry one.
fn first_format_count(r: &HtslibRecord, tag: &[u8]) -> u32 {
    r.format(tag)
        .integer()
        .ok()
        .and_then(|buffer| buffer.first().and_then(|values| values.first()).copied())
        .and_then(|n| u32::try_from(n).ok())
        .unwrap_or(0)
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
