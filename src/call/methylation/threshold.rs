use crate::{
    utils::Base,
    vcf::{self, Methylated},
};
use color_eyre::{Result, Section, eyre::ContextCompat};
use smol_str::SmolStr;
use tracing::{instrument, trace};

#[derive(Debug, Clone, clap::Args)]
pub struct ThresholdConfig {
    /// The minimum VAF to call a variant
    #[clap(long, default_value_t = 0.)]
    pub vaf_min: f64,

    /// The minimum number of reads to call a variant
    #[clap(long, default_value_t = 0)]
    pub reads_min: usize,
}

#[instrument(level="trace", skip(record, config), fields(
    chr = %record.fixed_fields.chrom,
    pos = record.fixed_fields.pos,
), name = "methylation_call")]
pub fn call(mut record: vcf::Record, config: &ThresholdConfig) -> Result<vcf::Record> {
    match call_methylation(&record, config)? {
        MethylationEvent::NotACpG => {
            // trace!("Not a CpG site, skipping");
        }
        MethylationEvent::Found(beta) => {
            trace!(beta, "Methylation event found");
            record.samples[0].methylated = Methylated(Some(beta));
            if beta > 1.0 {
                // fixme: this disables it right now
                record.fixed_fields.alt = smallvec::smallvec![".".into()];
            }
        }
        MethylationEvent::NotFound { failed_at } => {
            trace!(%failed_at, "Not methylated");
            record.samples[0].methylated = Methylated(Some(0.));
        }
    }
    Ok(record)
}

fn call_methylation(record: &vcf::Record, config: &ThresholdConfig) -> Result<MethylationEvent> {
    if record.fixed_fields.r#ref == "C" {
        call_position(record, config, Base::C, Base::T)
    } else if record.fixed_fields.r#ref == "G" {
        call_position(record, config, Base::G, Base::A)
    } else {
        Ok(MethylationEvent::NotACpG)
    }
}

#[derive(Debug)]
enum MethylationEvent {
    NotACpG,
    NotFound {
        failed_at: SmolStr,
    },
    /// `CpG` methylation event found, give beta value
    ///
    /// `alt_count/(alt_count+ref_count)` for OT (in case of ref `C`) or OB (in case of ref `G`)
    Found(f64),
}

impl MethylationEvent {
    /// Returns true if the event is a methylation event
    pub fn no(t: impl Into<SmolStr>) -> MethylationEvent {
        MethylationEvent::NotFound { failed_at: t.into() }
    }
}

fn call_position(
    record: &vcf::Record,
    config: &ThresholdConfig,
    ref_base: Base,
    alt_base: Base,
) -> Result<MethylationEvent> {
    if !record.info.in_cp_g.0 {
        return Ok(MethylationEvent::no("Not a CpG site"));
    }

    if !record.fixed_fields.alt.iter().any(|alt| alt == alt_base.as_str()) {
        return Ok(MethylationEvent::no("No T base in alts"));
    }

    if *record.info.read_depth < config.reads_min {
        return Ok(MethylationEvent::no("Not enough reads"));
    }

    // Check if the VAF is above the minimum threshold
    let t_alt_idx = record
        .fixed_fields
        .alt
        .iter()
        .position(|b| b == alt_base.as_str())
        .wrap_err("Alt base should be present in alts after previous checks")
        .note("This is a program error")?;
    if *record
        .info
        .allel_frequency
        .get(t_alt_idx)
        .wrap_err("Failed to get alt base in VAF")
        .note("This is a program error")?
        < config.vaf_min
    {
        return Ok(MethylationEvent::no("VAF below minimum threshold"));
    }

    let beta = {
        let refs = record
            .info
            .allele_specific_strand_bias
            .iter()
            .find(|x| x.base == ref_base)
            .wrap_err("allele specific strand bias should have ref allele")
            .note("This is a program error")?;

        let alts = record
            .info
            .allele_specific_strand_bias
            .iter()
            .find(|x| x.base == alt_base)
            .wrap_err("Failed to get alt base in allele specific strand bias")
            .note("This is a program error")?;

        let (refs, alts) = if ref_base == Base::C {
            // ref C: alt is T and we need to look at OT only
            (refs.ot, alts.ot)
        } else {
            // ref G: alt is A and we need to look at OB only
            (refs.ob, alts.ob)
        };

        f64::from(alts) / f64::from(alts + refs)
    };
    Ok(MethylationEvent::Found(beta))
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
        let metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        let config = ThresholdConfig { vaf_min: 0.1, reads_min: 5 };

        assert_debug_snapshot!((
            pileup.chrom(),
            pileup.pos,
            pileup.reference_base,
            &metrics.info.allele_specific_strand_bias,
            call_methylation(&metrics, &config),
        ), @r#"
        (
            "chr19",
            6105084,
            C,
            AlleleSpecificStrandBias(
                [
                    StrandCounts {
                        base: C,
                        ot: 4,
                        ob: 19,
                    },
                    StrandCounts {
                        base: T,
                        ot: 14,
                        ob: 0,
                    },
                ],
            ),
            Ok(
                Found(
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
        let metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        let config = ThresholdConfig { vaf_min: 0.1, reads_min: 5 };

        assert_debug_snapshot!((
            pileup.chrom(),
            pileup.pos,
            pileup.reference_base,
            &metrics.info.allele_specific_strand_bias,
            call_methylation(&metrics, &config),
        ), @r#"
        (
            "chr19",
            6105085,
            G,
            AlleleSpecificStrandBias(
                [
                    StrandCounts {
                        base: G,
                        ot: 18,
                        ob: 4,
                    },
                    StrandCounts {
                        base: A,
                        ot: 0,
                        ob: 15,
                    },
                ],
            ),
            Ok(
                Found(
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
        let metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        let config = ThresholdConfig { vaf_min: 0.1, reads_min: 5 };

        assert_debug_snapshot!((
            pileup.chrom(),
            pileup.pos,
            pileup.reference_base,
            &metrics.info.allele_specific_strand_bias,
            call_methylation(&metrics, &config),
        ), @r#"
        (
            "chr19",
            6105197,
            C,
            AlleleSpecificStrandBias(
                [
                    StrandCounts {
                        base: C,
                        ot: 19,
                        ob: 13,
                    },
                    StrandCounts {
                        base: T,
                        ot: 1,
                        ob: 0,
                    },
                ],
            ),
            Ok(
                NotFound {
                    failed_at: "Not a CpG site",
                },
            ),
        )
        "#);
        Ok(())
    }

    #[test]
    fn random_other_variant() -> Result<()> {
        let pileup = variant_pileup("chr19", 6105114)?;
        let metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        let config = ThresholdConfig { vaf_min: 0.1, reads_min: 5 };

        assert_debug_snapshot!((
            pileup.chrom(),
            pileup.pos,
            pileup.reference_base,
            &metrics.info.allele_specific_strand_bias,
            call_methylation(&metrics, &config),
        ), @r#"
        (
            "chr19",
            6105114,
            A,
            AlleleSpecificStrandBias(
                [
                    StrandCounts {
                        base: A,
                        ot: 20,
                        ob: 18,
                    },
                    StrandCounts {
                        base: G,
                        ot: 0,
                        ob: 1,
                    },
                ],
            ),
            Ok(
                NotACpG,
            ),
        )
        "#);
        Ok(())
    }

    #[test]
    fn methylatable_position_not_methylated() -> Result<()> {
        let pileup = variant_pileup("chr19", 6115809)?;
        let metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        let config = ThresholdConfig { vaf_min: 0.1, reads_min: 5 };

        let after = call(metrics, &config)?;

        assert_debug_snapshot!((
            pileup.chrom(),
            pileup.pos,
            pileup.reference_base,
            &after.info.allele_specific_strand_bias,
            &after.samples[0].methylated,
        ), @r#"
        (
            "chr19",
            6115809,
            C,
            AlleleSpecificStrandBias(
                [
                    StrandCounts {
                        base: C,
                        ot: 5,
                        ob: 4,
                    },
                    StrandCounts {
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
