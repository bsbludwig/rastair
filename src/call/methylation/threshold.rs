use crate::{
    utils::Base,
    vcf::{self, Methylated},
};
use color_eyre::{
    Result, Section,
    eyre::{Context, ContextCompat},
};
use rastair2_vcf::standard_fields::PASS;
use tracing::{instrument, trace};

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
pub fn call(mut record: vcf::Record, config: &ThresholdConfig) -> Result<vcf::Record> {
    if !record.info.in_cp_g.0 {
        trace!("Not a CpG site, skipping");
        return Ok(record);
    }

    let (ref_, alt) = if record.fixed_fields.r#ref == "C" {
        (Base::C, Base::T)
    } else if record.fixed_fields.r#ref == "G" {
        (Base::G, Base::A)
    } else {
        // trace!("Not a CpG site, skipping");
        return Ok(record);
    };

    match call_position(&record, ref_, alt).wrap_err("Failed to call position")? {
        MethylationEvent::Found { beta, others } => {
            trace!(beta, "Methylation event found");
            record.samples[0].methylated = Methylated(Some(beta));
            if record.fixed_fields.alt == [alt.as_str()] && others == 0 {
                record.fixed_fields.alt = smallvec::smallvec!["<*>".into()];
            }
            add_filters(&mut record, config, ref_, alt).wrap_err("Failed to add filters")?;
            if record.filters.is_empty() {
                record.filters.add(PASS);
            }
        }
        other => {
            trace!(?other, "Not methylated");
            record.samples[0].methylated = Methylated(Some(0.));
        }
    }

    Ok(record)
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

fn call_position(record: &vcf::Record, ref_base: Base, alt_base: Base) -> Result<MethylationEvent> {
    if !record.info.in_cp_g.0 {
        return Ok(MethylationEvent::NotACpG);
    }

    if !record.fixed_fields.alt.iter().any(|alt| alt == alt_base.as_str()) {
        return Ok(MethylationEvent::NoEvidenceByAlt);
    }

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

    let beta = {
        let (refs, alts) = if ref_base == Base::C {
            // ref C: alt is T and we need to look at OT only
            (refs.ot, alts.ot)
        } else {
            // ref G: alt is A and we need to look at OB only
            (refs.ob, alts.ob)
        };
        f64::from(alts) / f64::from(alts + refs)
    };

    let others = {
        let (_other_refs, other_alts) =
            if ref_base == Base::C { (refs.ob, alts.ob) } else { (refs.ot, alts.ot) };
        other_alts
    };

    Ok(MethylationEvent::Found { beta, others })
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
        let metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        let metrics = call(metrics, &config)?;

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
        let metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        let metrics = call(metrics, &config)?;

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
        let metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        let metrics = call(metrics, &config)?;

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
                    StrandCounts {
                        base: C,
                        ot: 20,
                        ob: 13,
                    },
                    StrandCounts {
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
        let metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        let metrics = call(metrics, &config)?;

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
        let metrics = pileup.variant_metrics(&VariantCallingParams::default())?;
        let metrics = call(metrics, &config)?;

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
