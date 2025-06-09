use crate::call::vcf::{self, Methylated};
use color_eyre::{Result, Section, eyre::ContextCompat};
use smallvec::SmallVec;
use smol_str::SmolStr;
use tracing::warn;

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
    let base_after =
        base_after(record.fixed_fields.r#ref.as_bytes()[0], &record.info.sequence_context[0]);
    let is_cpg = &record.fixed_fields.r#ref == "C" && base_after.filter(|x| x == "G").is_some();
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
    let read_depth_c = allels
        .iter()
        .position(|x| x == "C")
        .map(|idx| record.info.read_depth_per_allel[idx])
        .unwrap_or(0);
    let read_depth_t = allels
        .iter()
        .position(|x| x == "T")
        .map(|idx| record.info.read_depth_per_allel[idx])
        .unwrap_or(0);
    let read_depth_a = allels
        .iter()
        .position(|x| x == "A")
        .map(|idx| record.info.read_depth_per_allel[idx])
        .unwrap_or(0);
    let read_depth_g = allels
        .iter()
        .position(|x| x == "G")
        .map(|idx| record.info.read_depth_per_allel[idx])
        .unwrap_or(0);
    if read_depth_a + read_depth_g >= read_depth_c + read_depth_t {
        // More A and G bases than C and T, not likely methylation
        return Ok(record);
    }

    // It's a methylation event!
    update_record(record)
}

/// Update the record to reflect a methylation event.
fn update_record(mut record: vcf::Record) -> Result<vcf::Record> {
    // Set the methylation event flag
    record.samples[0].methylated = Methylated(1.);

    // Set alts to "missing"
    record.fixed_fields.alt = smallvec::smallvec![".".into()];

    Ok(record)
}

fn base_after(me: u8, context: &str) -> Option<SmolStr> {
    fn s(b: &u8) -> Option<SmolStr> {
        Some(SmolStr::new_inline(std::str::from_utf8(&[*b]).expect("constructed from `&str`")))
    }

    match context.as_bytes() {
        // usual case, 5-base context
        [_p2, _p1, mid, n1, _n2] if *mid == me => s(n1),
        // start of the sequence, no bases before
        [mid, n1, _n2] if *mid == me => s(n1),
        // start of the sequence, only one base before
        [_p1, mid, n1, _n2] if *mid == me => s(n1),
        // almost at the end of the sequence, only one base after
        [_p2, _p1, mid, n1] if *mid == me => s(n1),
        [_p2, _p1, mid] if *mid == me => {
            // we are at the end of the sequence, no base after
            None
        }
        _ => {
            warn!(len = context.len(), ?context, "Sequence context with unexpected length");
            None
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
