use crate::{
    call::methylation::threshold::ThresholdConfig,
    vcf::{self},
};
use color_eyre::Result;
use rastair2_vcf::standard_fields::PASS;
use tracing::trace;

pub fn add_filters(config: &ThresholdConfig, record: &mut vcf::Record) -> Result<()> {
    vcf::lowDP::apply(config, record);
    vcf::m_vaf::apply(config, record);

    let did_we_implement_all_filters = false; // just to keep track during development
    if did_we_implement_all_filters && record.filters.is_empty() {
        record.filters.add(PASS);
    }

    Ok(())
}

/// Check if this filter applies to the record.
trait CheckFilter: rastair2_vcf::VcfFilter {
    /// Check if the filter condition is met for the given record.
    fn check(_config: &ThresholdConfig, _record: &mut vcf::Record) -> bool;

    /// Apply the filter to the record if the condition is met.
    fn apply(config: &ThresholdConfig, record: &mut vcf::Record) {
        if Self::check(config, record) {
            record.filters.add(Self::default());
        }
    }
}

impl CheckFilter for vcf::lowDP {
    fn check(config: &ThresholdConfig, record: &mut vcf::Record) -> bool {
        *record.info.read_depth < config.m_min_depth
    }
}

impl CheckFilter for vcf::m_vaf {
    fn check(config: &ThresholdConfig, record: &mut vcf::Record) -> bool {
        let Some(alt) = record.m_base() else {
            return false;
        };
        let Some(alt_index) = record.fixed_fields.alt.iter().position(|a| a == alt.as_str()) else {
            return false;
        };
        let Some(vaf) = record.info.allele_frequency.get(alt_index) else {
            trace!("Alt allele {alt} not found in VAF info -- possibly a bug");
            return false;
        };

        *vaf < config.m_vaf_min
    }
}
