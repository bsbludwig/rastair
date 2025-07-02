use crate::{
    call::methylation::threshold::{ThresholdConfig, filters::add_filters},
    utils::Base,
    vcf::{self, Methylated},
};
use color_eyre::{
    Result,
    eyre::{Context, ContextCompat},
};
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

    // Helper function to get strand bias counts for a specific base
    let strand_count = |base: Base| -> Option<&vcf::ByStrand<u32>> {
        record.info.allele_specific_strand_bias.iter().find(|x| x.base == base)
    };

    // Check if alt contains "C" and ref_after is "G" (creating new CpG)
    if record.main.alt.iter().any(|alt| alt == "C") && ref_after == Some(Base::G) {
        if record.main.r#ref == "T" {
            // T > C case: need to use strand to distinguish mod from unmod
            let c_counts = strand_count(Base::C).wrap_err("Missing C counts in strand bias")?;
            let t_counts = strand_count(Base::T).wrap_err("Missing T counts in strand bias")?;

            // If there's 2+ reads evidence for T on OB, assume het SNP and adjust beta
            // Note that T is the _ref_ here
            // TODO: some more sophisticated SNP calling here, taking into account baseq, mapq etc
            if t_counts.ob >= 2 {
                // mod (reads showing T) are the ref here
                // divide by 2 assuming diploid genome
                let mod_count = f64::from(t_counts.ot) / 2.;
                let total = f64::from(c_counts.ot) + mod_count;
                if total > 0. {
                    Ok(Methylated::DeNovoCpG { beta: mod_count / total })
                } else {
                    Ok(Methylated::NoEvidence)
                }
            } else {
                let total = f64::from(c_counts.ot + t_counts.ot);
                if total > 0. {
                    Ok(Methylated::DeNovoCpG { beta: f64::from(t_counts.ot) / total })
                } else {
                    Ok(Methylated::NoEvidence)
                }
            }
        } else {
            // Ref is not T: count alt == T and alt == C separately
            let mod_count = if record.main.alt.iter().any(|alt| alt == "T") {
                strand_count(Base::T).map(|counts| counts.ot).unwrap_or(0)
            } else {
                0
            };

            let unmod_count = strand_count(Base::C).map(|counts| counts.ot).unwrap_or(0);

            // Check if there's evidence for T on the OB, which would be very
            // weird, ie a multi-allelic site (X->C _and_ X->T ?!)
            if let Some(t_counts) = strand_count(Base::T)
                && t_counts.ob > 0
            {
                warn!(
                    chr = %record.main.chrom,
                    pos = record.main.pos,
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
                Ok(Methylated::DeNovoCpG { beta: f64::from(mod_count) / f64::from(total) })
            }
        }
    }
    // Check if alt contains "G" and ref_before is "C" (creating new CpG)
    else if record.main.alt.iter().any(|alt| alt == "G") && ref_before == Some(Base::C) {
        if record.main.r#ref == "A" {
            // A > G case: similar logic but for OB strand
            let g_counts = strand_count(Base::G).wrap_err("Missing G counts in strand bias")?;
            let a_counts = strand_count(Base::A).wrap_err("Missing A counts in strand bias")?;

            // If there's 2+ reads evidence for A on OT, assume het SNP and adjust beta
            if a_counts.ot >= 2 {
                let mod_count = a_counts.ob / 2;
                let total = g_counts.ob + mod_count;
                if total > 0 {
                    Ok(Methylated::DeNovoCpG { beta: f64::from(mod_count) / f64::from(total) })
                } else {
                    Ok(Methylated::NoEvidence)
                }
            } else {
                let total = g_counts.ob + a_counts.ob;
                if total > 0 {
                    Ok(Methylated::DeNovoCpG { beta: f64::from(a_counts.ob) / f64::from(total) })
                } else {
                    Ok(Methylated::NoEvidence)
                }
            }
        } else {
            // Ref is not A: count alt == A and alt == G separately
            let mod_count = if record.main.alt.iter().any(|alt| alt == "A") {
                strand_count(Base::A).map(|counts| counts.ob).unwrap_or(0)
            } else {
                0
            };

            let unmod_count = strand_count(Base::G).map(|counts| counts.ob).unwrap_or(0);

            if tracing::enabled!(Level::DEBUG)
                && let Some(a_counts) = strand_count(Base::A)
                && a_counts.ot > 0
            {
                debug!(
                    chr = %record.main.chrom,
                    pos = record.main.pos,
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
                Ok(Methylated::DeNovoCpG { beta: f64::from(mod_count) / f64::from(total) })
            }
        }
    }
    // Handle C or G positions next to variants
    else if ref_base == "C" {
        // Check for non-T alternatives (possible C->N SNP)
        if tracing::enabled!(Level::TRACE) && record.main.alt.iter().any(|alt| alt != "T") {
            trace!(
                chr = %record.main.chrom,
                pos = record.main.pos,
                "Possible C->N SNP next to a de-novo G"
            );
        }

        if record.main.alt.iter().any(|alt| alt == "T") {
            let t_counts = strand_count(Base::T).wrap_err("Missing T counts in strand bias")?;
            let c_counts = strand_count(Base::C).wrap_err("Missing C counts in strand bias")?;

            let mut mod_count = t_counts.ot;
            let unmod_count = c_counts.ot;
            let snp_count = t_counts.ob;

            if snp_count > 1 {
                mod_count /= 2;
            }

            let total = mod_count + unmod_count;
            if total > 0 {
                Ok(Methylated::OriginalCpG { beta: f64::from(mod_count) / f64::from(total) })
            } else {
                Ok(Methylated::NoEvidence)
            }
        } else {
            Ok(Methylated::NoEvidence)
        }
    } else if ref_base == "G" {
        // Check for non-A alternatives (possible G->N SNP)
        if tracing::enabled!(Level::TRACE) && record.main.alt.iter().any(|alt| alt != "A") {
            trace!(
                chr = %record.main.chrom,
                pos = record.main.pos,
                "Possible G->N SNP next to a de-novo C"
            );
        }

        if record.main.alt.iter().any(|alt| alt == "A") {
            let a_counts = strand_count(Base::A).wrap_err("Missing A counts in strand bias")?;
            let g_counts = strand_count(Base::G).wrap_err("Missing G counts in strand bias")?;

            let mut mod_count = a_counts.ob;
            let unmod_count = g_counts.ob;
            let snp_count = a_counts.ot;

            if snp_count > 1 {
                mod_count /= 2;
            }

            let total = mod_count + unmod_count;
            if total > 0 {
                Ok(Methylated::OriginalCpG { beta: f64::from(mod_count) / f64::from(total) })
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
