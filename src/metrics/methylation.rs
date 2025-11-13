use crate::{
    call::{methylation::ThresholdParams, variant_calling::GenotypeTag},
    metrics::PileupMetrics,
    utils::{Base::*, IntoF64, logging::ThisIsABug},
    vcf::Methylated,
};
use color_eyre::{
    Result,
    eyre::{Context, ContextCompat, eyre},
};
use rastair_types::Probability;
use tracing::{Level, debug, instrument, trace, warn};

#[instrument(
    level="debug",
    skip_all,
    fields(contig = %current.contig(), pos = current.pos()),
    name = "methylation_call"
)]
pub fn call(config: &ThresholdParams, current: &PileupMetrics) -> Result<Option<Methylated>> {
    if *current.pos_metrics.cpg || current.forms_denovo() {
        let res = call_methylation(config, current).wrap_err("Failed to call methylation")?;
        if let Some(beta) = res.beta() {
            if beta.is_finite() {
                Ok(Some(res))
            } else {
                Err(eyre!("Methylation calling resulted in non-finite beta value `{beta}`"))
                    .this_is_a_bug()
            }
        } else {
            warn!(?res, "Methylation calling resulted in no beta value");
            Ok(None)
        }
    } else {
        // Not a CpG site, skipping
        Ok(None)
    }
}

fn call_methylation(config: &ThresholdParams, p: &PileupMetrics) -> Result<Methylated> {
    let ref_base = p.ref_base();
    let sequence_context = &p.pileup.context;
    let ref_before = sequence_context.before_1;
    let ref_after = sequence_context.after_1;

    if ref_base == C {
        ref_c(config, p)
    } else if ref_base == G {
        ref_g(config, p)
    } else if p.alt(C).is_some() && ref_after == G {
        // creating new CpG
        if ref_base == T { ref_t_to_c(config, p) } else { ref_not_t_to_c(config, p) }
    } else if p.alt(G).is_some() && ref_before == C {
        // creating new CpG
        if ref_base == A { ref_a_to_g(config, p) } else { ref_not_a_to_g(config, p) }
    } else {
        // Getting here should be impossible by construction
        warn!(
            ?ref_base,
            ?ref_before,
            ?ref_after,
            cpg=?p.pos_metrics.cpg,
            denovo=?p.pos_metrics.denovo_adj,
            "Position is neither an original CpG nor a de-novo CpG site"
        );
        Ok(Methylated::NoEvidence)
    }
}

fn ref_t_to_c(config: &ThresholdParams, record: &PileupMetrics) -> Result<Methylated> {
    assert_eq!(record.ref_base(), T);
    let t = &record.ref_metrics;
    let c = record.alt(C).wrap_err("Expected alt C at T->C denovo CpG site").this_is_a_bug()?;

    // T > C case: need to use strand to distinguish mod from unmod
    let c_counts = c.strand_count;
    let t_counts = t.strand_count;

    // If there's 2+ reads evidence for T on OB, assume het SNP and adjust beta
    // Note that T is the _ref_ here
    // TODO: some more sophisticated SNP calling here, taking into account baseq, mapq etc
    if t_counts.ob >= config.m_min_denovo_depth {
        // mod (reads showing T) are the ref here
        // divide by 2 assuming diploid genome
        let mod_count = t_counts.ot.f() / 2.;
        let total = c_counts.ot.f() + mod_count;
        if total > 0. {
            Ok(Methylated::DeNovoCpG { beta: Probability::new(mod_count / total).this_is_a_bug()? })
        } else {
            Ok(Methylated::NoEvidence)
        }
    } else {
        let mod_count = t_counts.ot;
        let total = c_counts.ot + t_counts.ot;
        if total > 0 {
            Ok(Methylated::DeNovoCpG {
                beta: Probability::new(mod_count.f() / total.f()).this_is_a_bug()?,
            })
        } else {
            Ok(Methylated::NoEvidence)
        }
    }
}

fn ref_not_t_to_c(_config: &ThresholdParams, record: &PileupMetrics) -> Result<Methylated> {
    assert_ne!(record.ref_base(), T);

    let t = record.alt(T);
    let c = record.alt(C).wrap_err("Expected alt C at non-T->C denovo CpG site").this_is_a_bug()?;

    // Get counts for `alt == T`, default to 0
    let t_counts = t.map(|a| a.strand_count).unwrap_or_default();

    // Ref is not T: count alt == T and alt == C separately
    let mod_count = t_counts.ot;
    let unmod = c.strand_count;
    let unmod_count = unmod.ot;

    // Check if there's evidence for T on the OB, which would be very
    // weird, ie a multi-allelic site (X->C _and_ X->T ?!)
    if t_counts.ob > 0 {
        debug!(?t_counts, "Evidence for multi-allelic SNP at het D/C site");
    }

    let total = mod_count + unmod_count;
    if total > 0 {
        Ok(Methylated::DeNovoCpG {
            beta: Probability::new(mod_count.f() / total.f()).this_is_a_bug()?,
        })
    } else {
        Ok(Methylated::NoEvidence)
    }
}

fn ref_a_to_g(config: &ThresholdParams, record: &PileupMetrics) -> Result<Methylated> {
    assert_eq!(record.ref_base(), A);
    let a = &record.ref_metrics;
    let g = record.alt(G).wrap_err("Expected alt G at A->G denovo CpG site").this_is_a_bug()?;

    // A > G case: similar logic but for OB strand
    let g_counts = g.strand_count;
    let a_counts = a.strand_count;

    // If there's 2+ reads evidence for A on OT, assume het SNP and adjust beta
    if a_counts.ot >= config.m_min_denovo_depth {
        // divide by 2 assuming diploid genome
        let mod_count = a_counts.ob.f() / 2.;
        let total = g_counts.ob.f() + mod_count;
        if total > 0. {
            Ok(Methylated::DeNovoCpG { beta: Probability::new(mod_count / total).this_is_a_bug()? })
        } else {
            Ok(Methylated::NoEvidence)
        }
    } else {
        let mod_count = a_counts.ob;
        let total = g_counts.ob + mod_count;
        if total > 0 {
            Ok(Methylated::DeNovoCpG {
                beta: Probability::new(mod_count.f() / total.f()).this_is_a_bug()?,
            })
        } else {
            Ok(Methylated::NoEvidence)
        }
    }
}

fn ref_not_a_to_g(_config: &ThresholdParams, record: &PileupMetrics) -> Result<Methylated> {
    assert_ne!(record.ref_base(), A);

    let a = record.alt(A);
    let g = record.alt(G).wrap_err("Expected alt G at non-A->G denovo CpG site").this_is_a_bug()?;

    let a_counts = a.map(|a| a.strand_count).unwrap_or_default();

    // Ref is not A: count alt == A and alt == G separately
    let mod_count = a_counts.ob;
    let unmod = g.strand_count;
    let unmod_count = unmod.ob;

    if a_counts.ot > 0 {
        debug!(?a_counts, "Evidence for multi-allelic SNP at het H/G site");
    }

    let total = mod_count + unmod_count;
    if total > 0 {
        Ok(Methylated::DeNovoCpG {
            beta: Probability::new(mod_count.f() / total.f()).this_is_a_bug()?,
        })
    } else {
        Ok(Methylated::NoEvidence)
    }
}

fn ref_c(_config: &ThresholdParams, record: &PileupMetrics) -> Result<Methylated> {
    assert_eq!(record.ref_base(), C);
    let c = &record.ref_metrics;

    // Check for non-T alternatives (possible C->N SNP)
    if tracing::enabled!(Level::TRACE)
        && *record.pos_metrics.denovo_adj
        && record.alts().iter().any(|b| *b != T)
    {
        trace!(
            pos = %record.pos(),
            "Possible C->N SNP next to a de-novo G"
        );
    }

    if let Some(t) = record.alt(T) {
        let t_counts = t.strand_count;
        let c_counts = c.strand_count;

        let mut mod_count = t_counts.ot.f();
        let unmod_count = c_counts.ot.f();

        if let Some(gt) = record.pos_metrics.genotype
            && gt.genotype == GenotypeTag::CT
        {
            // divide by 2 assuming diploid genome
            mod_count /= 2.;
        }

        let total = mod_count + unmod_count;
        if total > 0. {
            Ok(Methylated::OriginalCpG {
                beta: Probability::new(mod_count / total).this_is_a_bug()?,
            })
        } else {
            Ok(Methylated::NoEvidence)
        }
    } else {
        Ok(Methylated::NoEvidence)
    }
}

fn ref_g(_config: &ThresholdParams, record: &PileupMetrics) -> Result<Methylated> {
    assert_eq!(record.ref_base(), G);
    let g = &record.ref_metrics;

    // Check for non-A alternatives (possible G->N SNP)
    if tracing::enabled!(Level::TRACE)
        && *record.pos_metrics.denovo_adj
        && record.alts().iter().any(|b| *b != A)
    {
        trace!(
            pos = %record.pos(),
            "Possible G->N SNP next to a de-novo C"
        );
    }

    if let Some(a) = record.alt(A) {
        let a_counts = a.strand_count;
        let g_counts = g.strand_count;

        let mut mod_count = a_counts.ob.f();
        let unmod_count = g_counts.ob.f();

        if let Some(gt) = record.pos_metrics.genotype
            && gt.genotype == GenotypeTag::CT
        {
            // divide by 2 assuming diploid genome
            mod_count /= 2.;
        }

        let total = mod_count + unmod_count;
        if total > 0. {
            Ok(Methylated::OriginalCpG {
                beta: Probability::new(mod_count / total).this_is_a_bug()?,
            })
        } else {
            Ok(Methylated::NoEvidence)
        }
    } else {
        Ok(Methylated::NoEvidence)
    }
}

// TODO: Rewrite methylation tests
// #[cfg(test)]
// mod tests {
//     use color_eyre::eyre::ContextCompat;

//     use super::*;
//     use crate::call::{test_helpers::variant_pileup, variant_calling::VariantCallingParams};

//     #[test]
//     fn test_beta_value_c_to_t() -> Result<()> {
//         let known_c_to_t_pos = variant_pileup("bacteriophage_lambda_CpG", 47482)?
//             .variant_metrics(&VariantCallingParams::default())?;
//         assert_eq!("C", known_c_to_t_pos.main.r#ref);
//         let known_g_pos = variant_pileup("bacteriophage_lambda_CpG", 47483)?
//             .variant_metrics(&VariantCallingParams::default())?;
//         assert_eq!("G", known_g_pos.main.r#ref);

//         let methylation = call_methylation(
//             &ThresholdParams::default(),
//             &known_c_to_t_pos,
//             None,
//             Some(&known_g_pos),
//         )?;

//         let mod_count = 12. / 2.;
//         let unmod_count = 1.;
//         let expected_beta = mod_count / (mod_count + unmod_count);
//         let actual_beta = methylation.beta().wrap_err("No beta value")?;

//         assert_eq!(expected_beta, actual_beta);

//         Ok(())
//     }

//     #[test]
//     fn test_beta_value_all_mod() -> Result<()> {
//         let known_c_to_t_pos = variant_pileup("bacteriophage_lambda_CpG", 42236)?
//             .variant_metrics(&VariantCallingParams::default())?;
//         assert_eq!("C", known_c_to_t_pos.main.r#ref);
//         let known_g_pos = variant_pileup("bacteriophage_lambda_CpG", 42237)?
//             .variant_metrics(&VariantCallingParams::default())?;
//         assert_eq!("G", known_g_pos.main.r#ref);

//         let methylation = call_methylation(&ThresholdParams::default(), &known_c_to_t_pos)?;

//         let mod_count = 7.;
//         let unmod_count = 0.;
//         let expected_beta = mod_count / (mod_count + unmod_count);
//         let actual_beta = methylation.beta().wrap_err("No beta value")?;

//         assert_eq!(expected_beta, actual_beta);

//         Ok(())
//     }

//     #[test]
//     fn test_beta_value_all_mod_a() -> Result<()> {
//         let known_c_pos = variant_pileup("bacteriophage_lambda_CpG", 14987)?
//             .variant_metrics(&VariantCallingParams::default())?;
//         assert_eq!("C", known_c_pos.main.r#ref);
//         let known_g_to_a_pos = variant_pileup("bacteriophage_lambda_CpG", 14988)?
//             .variant_metrics(&VariantCallingParams::default())?;
//         assert_eq!("G", known_g_to_a_pos.main.r#ref);

//         let methylation = call_methylation(
//             &ThresholdParams::default(),
//             &known_g_to_a_pos,
//             Some(&known_c_pos),
//             None,
//         )?;

//         let mod_count = 3.;
//         let unmod_count = 0.;
//         let expected_beta = mod_count / (mod_count + unmod_count);
//         let actual_beta = methylation.beta().wrap_err("No beta value")?;

//         assert_eq!(expected_beta, actual_beta);

//         Ok(())
//     }
// }
