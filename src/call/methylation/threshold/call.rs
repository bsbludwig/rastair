use crate::{
    call::methylation::threshold::{
        ThresholdConfig, filters::add_filters, utils::NoStrandBiasForBaseErrorExt as _,
    },
    utils::Base::*,
    vcf::{self, Methylated},
};
use color_eyre::{Result, eyre::Context};
use smol_str::SmolStr;
use tracing::{Level, debug, instrument, trace, warn};

#[instrument(level="trace", skip(record, config), fields(
chr = %record.main.chrom,
pos = record.main.pos,
), name = "methylation_call")]
pub fn call(
    config: &ThresholdConfig,
    record: &mut vcf::Record,
    before: Option<&vcf::Record>,
    after: Option<&vcf::Record>,
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
    if record.has_alt(C) && ref_after == Some(G) {
        if record.main.r#ref == T {
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
        } else {
            // Ref is not T: count alt == T and alt == C separately
            let mod_count = record.strand_count(T).or_default().ot;
            let unmod_count = record.strand_count(C).or_default().ot;

            // Check if there's evidence for T on the OB, which would be very
            // weird, ie a multi-allelic site (X->C _and_ X->T ?!)
            if let Ok(t_counts) = record.strand_count(T)
                && t_counts.ob > 0
            {
                warn!(
                    chr = %record.main.chrom,
                    pos = record.main.pos,
                    ?t_counts,
                    "Evidence for multi-allelic SNP at het D/C site"
                );
            }

            if unmod_count == 0 {
                warn!(
                    chr = %record.main.chrom,
                    pos = record.main.pos,
                    "No evidence for C - this should be impossible"
                );
                Ok(Methylated::NoEvidence)
            } else {
                let total = mod_count + unmod_count;
                Ok(Methylated::DeNovoCpG { beta: f(mod_count) / f(total) })
            }
        }
    }
    // Check if alt contains "G" and ref_before is "C" (creating new CpG)
    else if record.has_alt(G) && ref_before == Some(C) {
        if record.main.r#ref == A {
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
        } else {
            // Ref is not A: count alt == A and alt == G separately
            let mod_count = record.strand_count(A).or_default().ob;
            let unmod_count = record.strand_count(G).or_default().ob;

            if tracing::enabled!(Level::DEBUG)
                && let Ok(a_counts) = record.strand_count(A)
                && a_counts.ot > 0
            {
                debug!(
                    chr = %record.main.chrom,
                    pos = record.main.pos,
                    ?a_counts,
                    "Evidence for multi-allelic SNP at het H/G site"
                );
            }

            if unmod_count == 0 {
                trace!(
                    chr = %record.main.chrom,
                    pos = record.main.pos,
                    "No evidence for G"
                );
                Ok(Methylated::NoEvidence)
            } else {
                let total = mod_count + unmod_count;
                Ok(Methylated::DeNovoCpG { beta: f(mod_count) / f(total) })
            }
        }
    }
    // Handle C or G positions next to variants
    else if ref_base == C {
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
    } else if ref_base == G {
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

fn f(x: impl Into<f64>) -> f64 {
    x.into()
}
