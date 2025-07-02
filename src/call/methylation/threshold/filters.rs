use crate::{call::methylation::threshold::ThresholdConfig, vcf};
use color_eyre::Result;
use rastair2_vcf::standard_fields::PASS;

pub fn add_filters(config: &ThresholdConfig, record: &mut vcf::Record) -> Result<()> {
    vcf::lowDP::apply(config, record);

    Ok(())
}

trait CheckFilter: rastair2_vcf::VcfFilter + Default {
    fn check(_config: &ThresholdConfig, _record: &mut vcf::Record) -> bool {
        false
    }

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
