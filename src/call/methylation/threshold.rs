use crate::call::vcf::{self, Record};
use color_eyre::{
    Result, Section,
    eyre::{ContextCompat, bail},
};
use smallvec::SmallVec;
use smol_str::SmolStr;

#[derive(Debug, Clone, clap::Args)]
pub struct ThresholdConfig {
    /// The minimum VAF to call a variant
    #[clap(long, default_value_t = 0.1)]
    pub vaf_min: f64,

    /// The maximum binomial p-value to call a variant
    #[clap(long, default_value_t = 0.08)]
    pub binomial_max: f64,

    /// The minimum number of reads to call a variant
    #[clap(long, default_value_t = 5)]
    pub reads_min: usize,
}

pub fn call(record: vcf::Record, config: &ThresholdConfig) -> Result<vcf::Record> {
    let is_cpg =
        record.fixed_fields.r#ref == "C" && record.base_after()?.filter(|x| x == "G").is_some();
    if !is_cpg {
        // Not a CpG site, cannot be a methylation event
        return Ok(record);
    }

    let could_be_methylation_event = record.fixed_fields.alt.iter().any(|alt| alt == "T");
    if !could_be_methylation_event {
        // No T base found in alts, cannot be a methylation event
        return Ok(record);
    }

    if *record.info.read_depth < config.reads_min {
        // Not enough evidence for methylation
        return Ok(record);
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
        // VAF is below the minimum threshold
        return Ok(record);
    }

    // Check if more A and G bases than C and T
    let allels = std::iter::once(record.fixed_fields.r#ref.clone())
        .chain(record.fixed_fields.alt.iter().cloned())
        .collect::<SmallVec<SmolStr, 4>>();
    let read_depth_c = record.info.read_depth_per_allel
        [allels.iter().position(|x| x == "C").expect("C should be present in allels")];
    let read_depth_t = record.info.read_depth_per_allel
        [allels.iter().position(|x| x == "T").expect("C should be present in allels")];
    let read_depth_a = record.info.read_depth_per_allel
        [allels.iter().position(|x| x == "A").expect("C should be present in allels")];
    let read_depth_g = record.info.read_depth_per_allel
        [allels.iter().position(|x| x == "G").expect("C should be present in allels")];
    if read_depth_a + read_depth_g >= read_depth_c + read_depth_t {
        // More A and G bases than C and T, not likely methylation
        return Ok(record);
    }

    Ok(record)
}

impl Record {
    fn base_after(&self) -> Result<Option<SmolStr>> {
        let me = &self.fixed_fields.r#ref.as_bytes()[0];
        let context = &self.info.sequence_context[0];
        match context.as_bytes() {
            [_p2, _p1, mid, _n1, _n2] if mid == me => {
                Ok(Some(context.get(2..3).wrap_err("index exists")?.into()))
            }
            [mid, _n1, _n2] if mid == me => {
                Ok(Some(context.get(1..2).wrap_err("index exists")?.into()))
            }
            [_p2, _p1, mid] if mid == me => {
                // we are at the end of the sequence, no base after
                Ok(None)
            }
            _ => bail!(
                "Sequence context with unexpected length {}: {:?}",
                self.info.sequence_context.0.len(),
                self.info.sequence_context
            ),
        }
    }
}

// todo: test this!
// - create builder for records, maybe not here
// - test: only consider CpG sites
//   - test: Record::base_after
// - test: only consider T base in alts
// - test: check VAF threshold
// - test: check read depth threshold
// - test: check A and G bases vs C and T bases
