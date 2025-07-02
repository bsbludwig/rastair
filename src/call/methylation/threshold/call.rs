use crate::{
    call::methylation::threshold::{
        ThresholdConfig, filters::add_filters, utils::NoStrandBiasForBaseErrorExt as _,
    },
    utils::Base::*,
    vcf::{self, Methylated},
};
use color_eyre::{Result, Section, eyre::Context};
use smol_str::SmolStr;
use tracing::{Level, debug, instrument, trace, warn};

#[instrument(
    level="debug",
    skip_all,
    fields(chr = %record.main.chrom, pos = record.main.pos),
    name = "methylation_call"
)]
pub fn call(
    config: &ThresholdConfig,
    record: &mut vcf::Record,
    _before: Option<&vcf::Record>,
    _after: Option<&vcf::Record>,
) -> Result<()> {
    if *record.info.in_cp_g || *record.info.de_novo_cp_g_candidate {
        record.samples[0].methylated = call_methylation(record, record.main.r#ref.clone())
            .wrap_err("Failed to call de novo CpG methylation")?;
        add_filters(config, record).wrap_err("Failed to add filters for CpG methylation")?;
    } else {
        trace!("Not a CpG site, skipping");
        return Ok(());
    };

    Ok(())
}

fn call_methylation(record: &vcf::Record, ref_base: SmolStr) -> Result<Methylated> {
    let sequence_context = &record.info.sequence_context;
    let ref_before = sequence_context.before_1;
    let ref_after = sequence_context.after_1;

    // Check if alt contains "C" and ref_after is "G" (creating new CpG)
    if record.has_alt(C) && ref_after == G {
        if record.main.r#ref == T { ref_t_to_c(record) } else { ref_not_t_to_c(record) }
    }
    // Check if alt contains "G" and ref_before is "C" (creating new CpG)
    else if record.has_alt(G) && ref_before == C {
        if record.main.r#ref == A { ref_a_to_g(record) } else { ref_not_a_to_g(record) }
    }
    // Handle C or G positions next to variants
    else if ref_base == C {
        ref_c(record)
    } else if ref_base == G {
        ref_g(record)
    } else {
        warn!(
            chr = %record.main.chrom,
            pos = record.main.pos,
            ref_base = ?ref_base,
            "Neither C nor G as ref, but also not a SNP - this should be impossible"
        );
        Ok(Methylated::NoEvidence)
    }
}

fn ref_t_to_c(record: &vcf::Record) -> Result<Methylated> {
    // T > C case: need to use strand to distinguish mod from unmod
    let c_counts = record.strand_count(C)?;
    let t_counts = record.strand_count(T)?;

    // If there's 2+ reads evidence for T on OB, assume het SNP and adjust beta
    // Note that T is the _ref_ here
    // TODO: some more sophisticated SNP calling here, taking into account baseq, mapq etc
    if t_counts.ob >= 2 {
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
        let total = c_counts.ot + t_counts.ot;
        if total > 0 {
            Ok(Methylated::DeNovoCpG { beta: f(t_counts.ot) / f(total) })
        } else {
            Ok(Methylated::NoEvidence)
        }
    }
}

fn ref_not_t_to_c(record: &vcf::Record) -> Result<Methylated> {
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
        warn!(?t_counts, "Evidence for multi-allelic SNP at het D/C site");
    }

    let total = mod_count + unmod_count;
    Ok(Methylated::DeNovoCpG { beta: f(mod_count) / f(total) })
}

fn ref_a_to_g(record: &vcf::Record) -> Result<Methylated> {
    // A > G case: similar logic but for OB strand
    let g_counts = record.strand_count(G)?;
    let a_counts = record.strand_count(A)?;

    // If there's 2+ reads evidence for A on OT, assume het SNP and adjust beta
    if a_counts.ot >= 2 {
        // divide by 2 assuming diploid genome
        let mod_count = f(a_counts.ob) / 2.;
        let total = f(g_counts.ob) + mod_count;
        if total > 0. {
            Ok(Methylated::DeNovoCpG { beta: mod_count / total })
        } else {
            Ok(Methylated::NoEvidence)
        }
    } else {
        let total = g_counts.ob + a_counts.ob;
        if total > 0 {
            Ok(Methylated::DeNovoCpG { beta: f(a_counts.ob) / f(total) })
        } else {
            Ok(Methylated::NoEvidence)
        }
    }
}

fn ref_not_a_to_g(record: &vcf::Record) -> Result<Methylated> {
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
    Ok(Methylated::DeNovoCpG { beta: f(mod_count) / f(total) })
}

fn ref_c(record: &vcf::Record) -> Result<Methylated> {
    // Check for non-T alternatives (possible C->N SNP)
    if tracing::enabled!(Level::TRACE) && record.has_alts_other_than(T) {
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

        let snp_count = t_counts.ob;
        if snp_count > 1 {
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

fn ref_g(record: &vcf::Record) -> Result<Methylated> {
    // Check for non-A alternatives (possible G->N SNP)
    if tracing::enabled!(Level::TRACE) && record.has_alts_other_than(A) {
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
        let snp_count = f(a_counts.ot);

        if snp_count > 1. {
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
