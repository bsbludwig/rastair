use crate::{
    call::methylation::{ThresholdParams, filters::add_filters},
    utils::Base::*,
    vcf::{self, Methylated, utils::NoStrandBiasForBaseErrorExt as _},
};
use color_eyre::{
    Result, Section,
    eyre::{Context, eyre},
};
use tracing::{Level, debug, instrument, trace, warn};

#[instrument(
    level="debug",
    skip_all,
    fields(chr = %record.main.chrom, pos = record.main.pos),
    name = "methylation_call"
)]
pub fn call(
    config: &ThresholdParams,
    record: &mut vcf::Record,
    before: Option<&vcf::Record>,
    after: Option<&vcf::Record>,
) -> Result<()> {
    if *record.info.in_cp_g || *record.info.de_novo_cp_g_candidate {
        let res = call_methylation(config, record, before, after)
            .wrap_err("Failed to call CpG methylation")?;
        if let Some(beta) = res.beta()
            && !beta.is_finite()
        {
            warn!(?res, "Methylation calling resulted in non-finite beta value")
        }
        record.samples[0].methylated = res;
        add_filters(config, record).wrap_err("Failed to add filters for CpG methylation")?;
    } else {
        trace!("Not a CpG site, skipping");
        return Ok(());
    };

    Ok(())
}

fn call_methylation(
    config: &ThresholdParams,
    record: &vcf::Record,
    _before: Option<&vcf::Record>,
    _after: Option<&vcf::Record>,
) -> Result<Methylated> {
    let ref_base = record.main.r#ref.clone();
    let sequence_context = &record.info.sequence_context;
    let ref_before = sequence_context.before_1;
    let ref_after = sequence_context.after_1;

    if record.has_alt(C) && ref_after == G {
        // creating new CpG
        if record.main.r#ref == T {
            ref_t_to_c(config, record)
        } else {
            ref_not_t_to_c(config, record)
        }
    } else if record.has_alt(G) && ref_before == C {
        // creating new CpG
        if record.main.r#ref == A {
            ref_a_to_g(config, record)
        } else {
            ref_not_a_to_g(config, record)
        }
    } else if ref_base == C {
        ref_c(config, record)
    } else if ref_base == G {
        ref_g(config, record)
    } else {
        // Getting here should be impossible by construction
        Err(eyre!("Neither C nor G as ref, but also not a SNP")).note("This is a programming error")
    }
}

fn ref_t_to_c(config: &ThresholdParams, record: &vcf::Record) -> Result<Methylated> {
    // T > C case: need to use strand to distinguish mod from unmod
    let c_counts = record.strand_count(C)?;
    let t_counts = record.strand_count(T)?;

    // If there's 2+ reads evidence for T on OB, assume het SNP and adjust beta
    // Note that T is the _ref_ here
    // TODO: some more sophisticated SNP calling here, taking into account baseq, mapq etc
    if t_counts.ob >= config.m_min_denovo_depth {
        // mod (reads showing T) are the ref here
        // divide by 2 assuming diploid genome
        let mod_count = f(t_counts.ot) / 2.;
        let total = f(c_counts.ot) + mod_count;
        if total > 0. {
            Ok(Methylated::DeNovoCpG { beta: mod_count / total })
        } else {
            Ok(Methylated::NoEvidence)
        }
    } else {
        let mod_count = t_counts.ot;
        let total = c_counts.ot + t_counts.ot;
        if total > 0 {
            Ok(Methylated::DeNovoCpG { beta: f(mod_count) / f(total) })
        } else {
            Ok(Methylated::NoEvidence)
        }
    }
}

fn ref_not_t_to_c(_config: &ThresholdParams, record: &vcf::Record) -> Result<Methylated> {
    // Ref is not T: count alt == T and alt == C separately
    let mod_count = record.strand_count(T).or_empty().ot;
    let unmod =
        record.strand_count(C).wrap_err("No evidence for C").note("This is a programming error")?;
    let unmod_count = unmod.ot;

    // Check if there's evidence for T on the OB, which would be very
    // weird, ie a multi-allelic site (X->C _and_ X->T ?!)
    if let Ok(t_counts) = record.strand_count(T)
        && t_counts.ob > 0
    {
        debug!(?t_counts, "Evidence for multi-allelic SNP at het D/C site");
    }

    let total = mod_count + unmod_count;
    if total > 0 {
        Ok(Methylated::DeNovoCpG { beta: f(mod_count) / f(total) })
    } else {
        Ok(Methylated::NoEvidence)
    }
}

fn ref_a_to_g(config: &ThresholdParams, record: &vcf::Record) -> Result<Methylated> {
    // A > G case: similar logic but for OB strand
    let g_counts = record.strand_count(G)?;
    let a_counts = record.strand_count(A)?;

    // If there's 2+ reads evidence for A on OT, assume het SNP and adjust beta
    if a_counts.ot >= config.m_min_denovo_depth {
        // divide by 2 assuming diploid genome
        let mod_count = f(a_counts.ob) / 2.;
        let total = f(g_counts.ob) + mod_count;
        if total > 0. {
            Ok(Methylated::DeNovoCpG { beta: mod_count / total })
        } else {
            Ok(Methylated::NoEvidence)
        }
    } else {
        let mod_count = a_counts.ob;
        let total = g_counts.ob + mod_count;
        if total > 0 {
            Ok(Methylated::DeNovoCpG { beta: f(mod_count) / f(total) })
        } else {
            Ok(Methylated::NoEvidence)
        }
    }
}

fn ref_not_a_to_g(_config: &ThresholdParams, record: &vcf::Record) -> Result<Methylated> {
    // Ref is not A: count alt == A and alt == G separately
    let mod_count = record.strand_count(A).or_empty().ob;
    let unmod =
        record.strand_count(G).wrap_err("No evidence for G").note("This is a programming error")?;
    let unmod_count = unmod.ob;

    if tracing::enabled!(Level::DEBUG)
        && let Ok(a_counts) = record.strand_count(A)
        && a_counts.ot > 0
    {
        debug!(?a_counts, "Evidence for multi-allelic SNP at het H/G site");
    }

    let total = mod_count + unmod_count;
    if total > 0 {
        Ok(Methylated::DeNovoCpG { beta: f(mod_count) / f(total) })
    } else {
        Ok(Methylated::NoEvidence)
    }
}

fn ref_c(_config: &ThresholdParams, record: &vcf::Record) -> Result<Methylated> {
    // Check for non-T alternatives (possible C->N SNP)
    if tracing::enabled!(Level::TRACE)
        && *record.info.de_novo_cp_g_candidate
        && record.has_alts_other_than(T)
    {
        trace!(
            chr = %record.main.chrom,
            pos = record.main.pos,
            "Possible C->N SNP next to a de-novo G"
        );
    }

    if record.has_alt(T) {
        let t_counts = record.strand_count(T)?;
        let c_counts = record.strand_count(C)?;

        let mut mod_count = f(t_counts.ot);
        let unmod_count = f(c_counts.ot);

        if record.samples[0].genotype.heterozygous() {
            // divide by 2 assuming diploid genome
            mod_count /= 2.;
        }

        let total = mod_count + unmod_count;
        if total > 0. {
            Ok(Methylated::OriginalCpG { beta: mod_count / total })
        } else {
            Ok(Methylated::NoEvidence)
        }
    } else {
        Ok(Methylated::NoEvidence)
    }
}

fn ref_g(_config: &ThresholdParams, record: &vcf::Record) -> Result<Methylated> {
    // Check for non-A alternatives (possible G->N SNP)
    if tracing::enabled!(Level::TRACE)
        && *record.info.de_novo_cp_g_candidate
        && record.has_alts_other_than(A)
    {
        trace!(
            chr = %record.main.chrom,
            pos = record.main.pos,
            "Possible G->N SNP next to a de-novo C"
        );
    }

    if record.has_alt(A) {
        let a_counts = record.strand_count(A)?;
        let g_counts = record.strand_count(G)?;

        let mut mod_count = f(a_counts.ob);
        let unmod_count = f(g_counts.ob);

        if record.samples[0].genotype.heterozygous() {
            // divide by 2 assuming diploid genome
            mod_count /= 2.;
        }

        let total = mod_count + unmod_count;
        if total > 0. {
            Ok(Methylated::OriginalCpG { beta: mod_count / total })
        } else {
            Ok(Methylated::NoEvidence)
        }
    } else {
        Ok(Methylated::NoEvidence)
    }
}

fn f(x: impl Into<f64>) -> f64 {
    x.into()
}
