use crate::{
    utils::Base,
    vcf::{self, Methylated},
};
use color_eyre::{
    Result, Section,
    eyre::{Context, ContextCompat},
};
use rastair2_vcf::standard_fields::PASS;
use smol_str::SmolStr;
use tracing::{instrument, trace, warn};

#[derive(Debug, Clone, clap::Args)]
pub struct ThresholdConfig {
    /// The minimum variant allele frequency
    #[clap(long, default_value_t = 0.)]
    pub vaf_min: f64,

    /// The minimum number of reads to call a variant
    #[clap(long, default_value_t = 3)]
    pub reads_min: usize,
}

#[instrument(level="trace", skip(record, config), fields(
    chr = %record.fixed_fields.chrom,
    pos = record.fixed_fields.pos,
), name = "methylation_call")]
pub fn call(
    config: &ThresholdConfig,
    record: &mut vcf::Record,
    before: Option<&vcf::Record>,
    after: Option<&vcf::Record>,
) -> Result<()> {
    if record.info.in_cp_g.0 || record.info.de_novo_cp_g_candidate.0 {
        call_cpg(config, record, before, after).wrap_err("Failed to call CpG methylation")?;
    } else {
        trace!("Not a CpG site, skipping");
        return Ok(());
    };

    Ok(())
}

fn add_filters(
    record: &mut vcf::Record,
    config: &ThresholdConfig,
    _ref_base: Base,
    _alt_base: Base,
) -> Result<()> {
    if *record.info.read_depth < config.reads_min {
        record.filters.add(vcf::lowDP);
    }

    Ok(())
}

#[derive(Debug)]
enum MethylationEvent {
    /// `CpG` methylation event found, give beta value
    ///
    /// `alt_count/(alt_count+ref_count)` for OT (in case of ref `C`) or OB (in case of ref `G`)
    Found {
        beta: f64,
        others: u32,
    },
    NotACpG,
    NoEvidenceByAlt,
}

pub fn call_cpg(
    config: &ThresholdConfig,
    record: &mut vcf::Record,
    _before: Option<&vcf::Record>,
    _after: Option<&vcf::Record>,
) -> Result<()> {
    let ref_base = record.fixed_fields.r#ref.clone();

    // todo: add filters
    let beta =
        call_methylation(record, ref_base).wrap_err("Failed to call de novo CpG methylation")?;

    if let Some(beta_value) = beta {
        trace!(beta = beta_value, "De novo CpG methylation event found");
        record.samples[0].methylated = Methylated(Some(beta_value));

        // TODO: add filters along the way
        // if record.filters.is_empty() {
        //     record.filters.add(PASS);
        // }
    } else {
        trace!("No de novo CpG methylation detected");
        record.samples[0].methylated = Methylated(Some(0.0));
    }

    Ok(())
}

fn call_methylation(record: &vcf::Record, ref_base: SmolStr) -> Result<Option<f64>> {
    let sequence_context = &record.info.sequence_context;
    let ref_before = sequence_context.before_1;
    let ref_after = sequence_context.after_1;

    // Helper function to get strand bias counts for a specific base
    let strand_count = |base: Base| -> Option<&vcf::ByStrand<u32>> {
        record.info.allele_specific_strand_bias.iter().find(|x| x.base == base)
    };

    // Check if alt contains "C" and ref_after is "G" (creating new CpG)
    if record.fixed_fields.alt.iter().any(|alt| alt == "C") && ref_after == Some(Base::G) {
        if record.fixed_fields.r#ref == "T" {
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
                if total > 0. { Ok(Some(mod_count / total)) } else { Ok(Some(0.0)) }
            } else {
                let total = f64::from(c_counts.ot + t_counts.ot);
                if total > 0. { Ok(Some(f64::from(t_counts.ot) / total)) } else { Ok(Some(0.0)) }
            }
        } else {
            // Ref is not T: count alt == T and alt == C separately
            let mod_count = if record.fixed_fields.alt.iter().any(|alt| alt == "T") {
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
                    chr = %record.fixed_fields.chrom,
                    pos = record.fixed_fields.pos,
                    "Evidence for multi-allelic SNP at het D/C site"
                );
            }

            if unmod_count == 0 {
                warn!(
                    chr = %record.fixed_fields.chrom,
                    pos = record.fixed_fields.pos,
                    "No evidence for C - this should be impossible"
                );
                Ok(None)
            } else {
                let total = mod_count + unmod_count;
                Ok(Some(f64::from(mod_count) / f64::from(total)))
            }
        }
    }
    // Check if alt contains "G" and ref_before is "C" (creating new CpG)
    else if record.fixed_fields.alt.iter().any(|alt| alt == "G") && ref_before == Some(Base::C) {
        if record.fixed_fields.r#ref == "A" {
            // A > G case: similar logic but for OB strand
            let g_counts = strand_count(Base::G).wrap_err("Missing G counts in strand bias")?;
            let a_counts = strand_count(Base::A).wrap_err("Missing A counts in strand bias")?;

            // If there's 2+ reads evidence for A on OT, assume het SNP and adjust beta
            if a_counts.ot >= 2 {
                let mod_count = a_counts.ob / 2;
                let total = g_counts.ob + mod_count;
                if total > 0 {
                    Ok(Some(f64::from(mod_count) / f64::from(total)))
                } else {
                    Ok(Some(0.0))
                }
            } else {
                let total = g_counts.ob + a_counts.ob;
                if total > 0 {
                    Ok(Some(f64::from(a_counts.ob) / f64::from(total)))
                } else {
                    Ok(Some(0.0))
                }
            }
        } else {
            // Ref is not A: count alt == A and alt == G separately
            let mod_count = if record.fixed_fields.alt.iter().any(|alt| alt == "A") {
                strand_count(Base::A).map(|counts| counts.ob).unwrap_or(0)
            } else {
                0
            };

            let unmod_count = strand_count(Base::G).map(|counts| counts.ob).unwrap_or(0);

            if let Some(a_counts) = strand_count(Base::A)
                && a_counts.ot > 0
            {
                warn!(
                    chr = %record.fixed_fields.chrom,
                    pos = record.fixed_fields.pos,
                    "Evidence for multi-allelic SNP at het H/G site"
                );
            }

            if unmod_count == 0 {
                warn!(
                    chr = %record.fixed_fields.chrom,
                    pos = record.fixed_fields.pos,
                    "No evidence for G"
                );
                Ok(None)
            } else {
                let total = mod_count + unmod_count;
                Ok(Some(f64::from(mod_count) / f64::from(total)))
            }
        }
    }
    // Handle C or G positions next to variants
    else if ref_base == "C" {
        // Check for non-T alternatives (possible C->N SNP)
        if record.fixed_fields.alt.iter().any(|alt| alt != "T") {
            warn!(
                chr = %record.fixed_fields.chrom,
                pos = record.fixed_fields.pos,
                "Possible C->N SNP next to a de-novo G"
            );
        }

        if record.fixed_fields.alt.iter().any(|alt| alt == "T") {
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
                Ok(Some(f64::from(mod_count) / f64::from(total)))
            } else {
                Ok(Some(0.0))
            }
        } else {
            Ok(Some(0.0))
        }
    } else if ref_base == "G" {
        // Check for non-A alternatives (possible G->N SNP)
        if record.fixed_fields.alt.iter().any(|alt| alt != "A") {
            warn!(
                chr = %record.fixed_fields.chrom,
                pos = record.fixed_fields.pos,
                "Possible G->N SNP next to a de-novo C"
            );
        }

        if record.fixed_fields.alt.iter().any(|alt| alt == "A") {
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
                Ok(Some(f64::from(mod_count) / f64::from(total)))
            } else {
                Ok(Some(0.0))
            }
        } else {
            Ok(Some(0.0))
        }
    } else {
        warn!(
            chr = %record.fixed_fields.chrom,
            pos = record.fixed_fields.pos,
            ref_base = ?ref_base,
            "Neither C nor G as ref, but also not a SNP - this should be impossible"
        );
        Ok(None)
    }
}

// todo: test this!
// - create builder for records, maybe not here
// - test: only consider CpG sites
// - test: only consider T base in alts
// - test: check VAF threshold
// - test: check read depth threshold
// - test: check A and G bases vs C and T bases
#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::{test_helpers::variant_pileup, variant_calling::VariantCallingParams};
    use color_eyre::Result;
    use insta::assert_debug_snapshot;

    #[test]
    fn cpg_c_methylated() -> Result<()> {
        let pileup = variant_pileup("chr19", 6105084)?;
        let config = ThresholdConfig { vaf_min: 0.1, reads_min: 5 };
        let mut metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        call(&config, &mut metrics, None, None)?;

        assert_debug_snapshot!((
            pileup.chrom(),
            pileup.pos,
            pileup.reference_base,
            &metrics.info.in_cp_g,
            &metrics.info.allele_specific_strand_bias,
            &metrics.samples[0].methylated,
        ), @r#"
        (
            "chr19",
            6105084,
            C,
            InCpG(
                true,
            ),
            AlleleSpecificStrandBias(
                [
                    ByStrand {
                        base: C,
                        ot: 4,
                        ob: 19,
                    },
                    ByStrand {
                        base: T,
                        ot: 14,
                        ob: 0,
                    },
                ],
            ),
            Methylated(
                Some(
                    0.7777777777777778,
                ),
            ),
        )
        "#);
        Ok(())
    }

    #[test]
    fn cpg_g_methylated() -> Result<()> {
        let pileup = variant_pileup("chr19", 6105085)?;
        let config = ThresholdConfig { vaf_min: 0.1, reads_min: 5 };
        let mut metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        call(&config, &mut metrics, None, None)?;

        assert_debug_snapshot!((
            pileup.chrom(),
            pileup.pos,
            pileup.reference_base,
            &metrics.info.in_cp_g,
            &metrics.info.allele_specific_strand_bias,
            &metrics.samples[0].methylated,
        ), @r#"
        (
            "chr19",
            6105085,
            G,
            InCpG(
                true,
            ),
            AlleleSpecificStrandBias(
                [
                    ByStrand {
                        base: G,
                        ot: 18,
                        ob: 4,
                    },
                    ByStrand {
                        base: A,
                        ot: 0,
                        ob: 15,
                    },
                ],
            ),
            Methylated(
                Some(
                    0.7894736842105263,
                ),
            ),
        )
        "#);
        Ok(())
    }

    #[test]
    fn c_but_not_cpg() -> Result<()> {
        let pileup = variant_pileup("chr19", 6105197)?;
        let config = ThresholdConfig { vaf_min: 0.1, reads_min: 5 };
        let mut metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        call(&config, &mut metrics, None, None)?;

        assert_debug_snapshot!((
            pileup.chrom(),
            pileup.pos,
            pileup.reference_base,
            &metrics.info.in_cp_g,
            &metrics.info.allele_specific_strand_bias,
            &metrics.samples[0].methylated,
        ), @r#"
        (
            "chr19",
            6105197,
            C,
            InCpG(
                false,
            ),
            AlleleSpecificStrandBias(
                [
                    ByStrand {
                        base: C,
                        ot: 20,
                        ob: 13,
                    },
                    ByStrand {
                        base: T,
                        ot: 1,
                        ob: 0,
                    },
                ],
            ),
            Methylated(
                None,
            ),
        )
        "#);
        Ok(())
    }

    #[test]
    fn random_other_variant() -> Result<()> {
        let pileup = variant_pileup("chr19", 6105114)?;
        let config = ThresholdConfig { vaf_min: 0.1, reads_min: 5 };
        let mut metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        call(&config, &mut metrics, None, None)?;

        assert_debug_snapshot!((
            pileup.chrom(),
            pileup.pos,
            pileup.reference_base,
            &metrics.info.in_cp_g,
            &metrics.info.allele_specific_strand_bias,
            &metrics.samples[0].methylated,
        ), @r#"
        (
            "chr19",
            6105114,
            A,
            InCpG(
                false,
            ),
            AlleleSpecificStrandBias(
                [
                    ByStrand {
                        base: A,
                        ot: 20,
                        ob: 18,
                    },
                    ByStrand {
                        base: G,
                        ot: 0,
                        ob: 1,
                    },
                ],
            ),
            Methylated(
                None,
            ),
        )
        "#);
        Ok(())
    }

    #[test]
    fn methylatable_position_not_methylated() -> Result<()> {
        let pileup = variant_pileup("chr19", 6115809)?;
        let config = ThresholdConfig { vaf_min: 0.1, reads_min: 5 };
        let mut metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        call(&config, &mut metrics, None, None)?;

        assert_debug_snapshot!((
            pileup.chrom(),
            pileup.pos,
            pileup.reference_base,
            &metrics.info.in_cp_g,
            &metrics.info.allele_specific_strand_bias,
            &metrics.samples[0].methylated,
        ), @r#"
        (
            "chr19",
            6115809,
            C,
            InCpG(
                true,
            ),
            AlleleSpecificStrandBias(
                [
                    ByStrand {
                        base: C,
                        ot: 5,
                        ob: 4,
                    },
                    ByStrand {
                        base: A,
                        ot: 1,
                        ob: 0,
                    },
                ],
            ),
            Methylated(
                Some(
                    0.0,
                ),
            ),
        )
        "#);

        Ok(())
    }
}
