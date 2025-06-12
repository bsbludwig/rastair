use crate::{
    utils::Base,
    vcf::{self, Methylated, Record},
};
use color_eyre::{Result, Section, eyre::ContextCompat};
use smallvec::SmallVec;
use smol_str::SmolStr;
use tracing::{instrument, trace};

#[derive(Debug, Clone, clap::Args)]
pub struct ThresholdConfig {
    /// The minimum VAF to call a variant
    #[clap(long, default_value_t = 0.1)]
    pub vaf_min: f64,

    /// The minimum number of reads to call a variant
    #[clap(long, default_value_t = 5)]
    pub reads_min: usize,
}

#[instrument(level="trace", skip(record, config), fields(
    chr = %record.fixed_fields.chrom,
    pos = record.fixed_fields.pos,
))]
pub fn call(mut record: vcf::Record, config: &ThresholdConfig) -> Result<vcf::Record> {
    match call_cpg(&record, config)? {
        MethylationEvent::Found(beta) => {
            record.samples[0].methylated = Methylated(beta);
            trace!(beta, "CpG methylation event found");
            return update_record(record);
        }
        MethylationEvent::NotFound { failed_at } => {
            trace!(failed = failed_at, "Not methylated");
        }
    }
    match call_gpc(&record, config)? {
        MethylationEvent::Found(beta) => {
            record.samples[0].methylated = Methylated(beta);
            trace!(beta, "CpG methylation event found");
            return update_record(record);
        }
        MethylationEvent::NotFound { failed_at } => {
            trace!(failed = failed_at, "Not methylated");
        }
    }

    // It's not a methylation event
    Ok(record)
}

#[derive(Debug)]
enum MethylationEvent {
    NotFound {
        failed_at: &'static str,
    },
    /// `CpG` methylation event found, give beta value
    ///
    /// `alt_count/(alt_count+ref_count)` for OT (in case of ref `C`) or OB (in case of ref `G`)
    Found(f64),
}

impl MethylationEvent {
    /// Returns true if the event is a methylation event
    pub fn no(failed_at: &'static str) -> MethylationEvent {
        MethylationEvent::NotFound { failed_at }
    }
}

struct Alleles(SmallVec<SmolStr, 4>);

impl Alleles {
    fn pos(&self, base: &str) -> Result<usize> {
        self.0
            .iter()
            .position(|x| x == base)
            .wrap_err_with(|| format!("Base {} not found in alleles", base))
    }
}

impl Record {
    fn allels(&self) -> Alleles {
        Alleles(
            std::iter::once(self.fixed_fields.r#ref.clone())
                .chain(self.fixed_fields.alt.iter().cloned())
                .collect::<SmallVec<SmolStr, 4>>(),
        )
    }
}

fn call_cpg(record: &vcf::Record, config: &ThresholdConfig) -> Result<MethylationEvent> {
    let context = &record.info.sequence_context;
    let is_cpg =
        &record.fixed_fields.r#ref == "C" && context.after_1.filter(|x| *x == Base::G).is_some();
    if !is_cpg {
        return Ok(MethylationEvent::no("Not a CpG site"));
    }

    if is_cpg && !record.fixed_fields.alt.iter().any(|alt| alt == "T") {
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
        .position(|b| b == "T")
        .wrap_err("T base should be present in alts after previous checks")
        .note("This is a program error")?;
    if *record
        .info
        .allel_frequency
        .get(t_alt_idx)
        .wrap_err("Failed to get T base in VAF")
        .note("This is a program error")?
        < config.vaf_min
    {
        return Ok(MethylationEvent::no("VAF below minimum threshold"));
    }

    // Check if more A and G bases than C and T
    let allels = record.allels();
    let read_depth_c =
        allels.pos("C").map(|idx| record.info.read_depth_per_allel[idx]).unwrap_or(0);
    let read_depth_t =
        allels.pos("T").map(|idx| record.info.read_depth_per_allel[idx]).unwrap_or(0);
    let read_depth_a =
        allels.pos("A").map(|idx| record.info.read_depth_per_allel[idx]).unwrap_or(0);
    let read_depth_g =
        allels.pos("G").map(|idx| record.info.read_depth_per_allel[idx]).unwrap_or(0);
    if read_depth_a + read_depth_g >= read_depth_c + read_depth_t {
        return Ok(MethylationEvent::no("More A and G bases than C and T, not likely methylation"));
    }

    // we're looking at ref C, so alt is T and we need to look at OT only
    let beta = {
        let ref_count = record
            .info
            .allele_specific_strand_bias
            .iter()
            .find(|x| x.base.as_str() == record.fixed_fields.r#ref) // the first allele is the reference
            .wrap_err("allele specific strand bias should have ref allele")
            .note("This is a program error")?
            .ot;

        let alt_count = record
            .info
            .allele_specific_strand_bias
            .iter()
            .find(|x| x.base.as_str() == "T")
            .wrap_err("Failed to get T base in allele specific strand bias")
            .note("This is a program error")?
            .ot;
        f64::from(alt_count) / f64::from(alt_count + ref_count)
    };
    Ok(MethylationEvent::Found(beta))
}

fn call_gpc(record: &vcf::Record, config: &ThresholdConfig) -> Result<MethylationEvent> {
    let context = &record.info.sequence_context;
    let is_reverse_cpg =
        &record.fixed_fields.r#ref == "G" && context.before_1.filter(|x| *x == Base::C).is_some();
    if !is_reverse_cpg {
        return Ok(MethylationEvent::no("Not a CpG site (reverse strand)"));
    }

    if !record.fixed_fields.alt.iter().any(|alt| alt == "A") {
        return Ok(MethylationEvent::no("No A base in alts"));
    }

    if *record.info.read_depth < config.reads_min {
        return Ok(MethylationEvent::no("Not enough reads"));
    }

    // Check if the VAF is above the minimum threshold
    let t_alt_idx = record
        .fixed_fields
        .alt
        .iter()
        .position(|b| b == "A")
        .wrap_err("T base should be present in alts after previous checks")
        .note("This is a program error")?;
    if *record
        .info
        .allel_frequency
        .get(t_alt_idx)
        .wrap_err("Failed to get T base in VAF")
        .note("This is a program error")?
        < config.vaf_min
    {
        return Ok(MethylationEvent::no("VAF below minimum threshold"));
    }

    // Check if more A and G bases than C and T
    let allels = record.allels();
    let read_depth_c =
        allels.pos("C").map(|idx| record.info.read_depth_per_allel[idx]).unwrap_or(0);
    let read_depth_t =
        allels.pos("T").map(|idx| record.info.read_depth_per_allel[idx]).unwrap_or(0);
    let read_depth_a =
        allels.pos("A").map(|idx| record.info.read_depth_per_allel[idx]).unwrap_or(0);
    let read_depth_g =
        allels.pos("G").map(|idx| record.info.read_depth_per_allel[idx]).unwrap_or(0);
    if read_depth_c + read_depth_t >= read_depth_a + read_depth_g {
        return Ok(MethylationEvent::no(
            "More C and T bases than A and G, not likely methylation on opposite strand",
        ));
    }

    // we're looking at ref G, so alt is A and we need to look at OB only
    let beta = {
        let ref_count = record
            .info
            .allele_specific_strand_bias
            .iter()
            .find(|x| x.base.as_str() == record.fixed_fields.r#ref)
            .wrap_err("allele specific strand bias should have ref allele")
            .note("This is a program error")?
            .ob;

        let alt_count = record
            .info
            .allele_specific_strand_bias
            .iter()
            .find(|x| x.base.as_str() == "A")
            .wrap_err("Failed to get A base in allele specific strand bias")
            .note("This is a program error")?
            .ob;
        trace!(?record
            .info
            .allele_specific_strand_bias, ?alt_count, ?ref_count, "beta");
        f64::from(alt_count) / f64::from(alt_count + ref_count)
    };
    Ok(MethylationEvent::Found(beta))
}

/// Update the record to reflect a methylation event.
fn update_record(mut record: vcf::Record) -> Result<vcf::Record> {
    // Set alts to "missing"
    record.fixed_fields.alt = smallvec::smallvec![".".into()];

    Ok(record)
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
    use crate::call::test_helpers::variant_pileup;
    use color_eyre::Result;
    use insta::assert_debug_snapshot;

    #[test]
    fn test_cpg_position() -> Result<()> {
        let pileup = variant_pileup("chr19", 6105084)?;
        let metrics = pileup.variant_metrics()?;
        let config = ThresholdConfig { vaf_min: 0.1, reads_min: 5 };

        assert_debug_snapshot!((
            pileup.chrom(),
            pileup.pos,
            pileup.reference_base,
            &pileup.bases,
            &metrics.info.allele_specific_strand_bias,
            call_cpg(&metrics, &config),
        ), @r#"
        (
            "chr19",
            6105084,
            C,
            [
                C OB Q36 MQ60,
                C OB Q36 MQ60,
                T OT Q36 MQ60,
                C OT Q32 MQ60,
                C OB Q36 MQ60,
                C OT Q36 MQ60,
                T OT Q36 MQ60,
                C OB Q36 MQ60,
                T OT Q36 MQ60,
                C OB Q36 MQ60,
                C OT Q36 MQ60,
                T OT Q36 MQ60,
                C OB Q36 MQ60,
                T OT Q36 MQ60,
                T OT Q36 MQ60,
                T OT Q36 MQ60,
                C OB Q36 MQ60,
                C OB Q36 MQ60,
                T OT Q36 MQ60,
                C OB Q36 MQ60,
                T OT Q36 MQ60,
                C OB Q36 MQ60,
                T OT Q36 MQ60,
                C OB Q36 MQ60,
                C OB Q36 MQ60,
                C OB Q36 MQ60,
                T OT Q36 MQ60,
                C OB Q36 MQ60,
                T OT Q36 MQ60,
                C OB Q36 MQ60,
                C OB Q36 MQ60,
                C OB Q36 MQ60,
                T OT Q36 MQ60,
                C OB Q36 MQ60,
                C OB Q36 MQ60,
                T OT Q36 MQ60,
                C OT Q36 MQ60,
            ],
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
}
