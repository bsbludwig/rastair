use crate::{
    call::methylation::threshold::ThresholdParams,
    vcf::{self},
};
use color_eyre::Result;
use tracing::trace;

/// Add filters
pub fn add_filters(config: &ThresholdParams, record: &mut vcf::Record) -> Result<()> {
    if !(*record.info.in_cp_g || *record.info.de_novo_cp_g_candidate) {
        // Not a CpG site, skipping filters
        return Ok(());
    }

    vcf::lowDp::apply(config, record);
    vcf::m_vaf::apply(config, record);
    vcf::m_bq_ratio::apply(config, record);
    vcf::m_pos::apply(config, record);
    vcf::m_highDp::apply(config, record);

    // If no filters were added, we're gonna call it
    if record.filters.is_empty() {
        record.filters.add(rastair2_vcf::standard_fields::PASS);
    }

    Ok(())
}

/// Check if this filter applies to the record.
trait CheckFilter: rastair2_vcf::VcfFilter {
    /// Check if the filter condition is met for the given record.
    fn check(_config: &ThresholdParams, _record: &mut vcf::Record) -> bool;

    /// Apply the filter to the record if the condition is met.
    fn apply(config: &ThresholdParams, record: &mut vcf::Record) {
        if Self::check(config, record) {
            record.filters.add(Self::default());
        }
    }
}

impl CheckFilter for vcf::lowDp {
    fn check(config: &ThresholdParams, record: &mut vcf::Record) -> bool {
        *record.info.read_depth < config.m_min_depth
    }
}

impl CheckFilter for vcf::m_vaf {
    fn check(config: &ThresholdParams, record: &mut vcf::Record) -> bool {
        let Some(alt) = record.m_base() else {
            return false;
        };
        let Some(alt_index) = record.main.alt.iter().position(|a| a == alt.as_str()) else {
            return false;
        };
        let Some(vaf) = record.info.allele_frequency.get(alt_index) else {
            trace!("Alt allele {alt} not found in VAF info -- possibly a bug");
            return false;
        };

        *vaf < config.m_vaf_min
    }
}

impl CheckFilter for vcf::m_bq_ratio {
    fn check(config: &ThresholdParams, record: &mut vcf::Record) -> bool {
        let Some(alt) = record.m_base() else {
            return false;
        };
        let Some(alt_index) = record.main.alt.iter().position(|a| a == alt.as_str()) else {
            return false;
        };

        // Get allele read depths (ref is at index 0, alts follow)
        let Some(ad_ref) = record.info.allele_read_depth.first() else {
            trace!("Missing ref allele depth");
            return false;
        };
        let Some(ad_alt) = record.info.allele_read_depth.get(alt_index + 1) else {
            trace!("Missing alt allele depth for index {}", alt_index);
            return false;
        };

        let Some(bq_ref) = record.info.allele_base_quality.first() else {
            trace!("Missing ref allele base quality");
            return false;
        };
        let Some(bq_alt) = record.info.allele_base_quality.get(alt_index + 1) else {
            trace!("Missing alt allele base quality for index {}", alt_index);
            return false;
        };

        let quality_ratio = ((*ad_alt as f64) * bq_alt + 1.0) / ((*ad_ref as f64) * bq_ref + 1.0);

        quality_ratio < config.m_bq_ratio_min
    }
}

impl CheckFilter for vcf::m_pos {
    fn check(config: &ThresholdParams, record: &mut vcf::Record) -> bool {
        let Some(alt) = record.m_base() else {
            return false;
        };
        let Some(alt_index) = record.main.alt.iter().position(|a| a == alt.as_str()) else {
            return false;
        };
        let Some(pos_in_read) = record.info.position_in_read.get(alt_index + 1) else {
            trace!("Missing position in read for alt allele at index {}", alt_index);
            return false;
        };

        // Check if position is outside the acceptable range
        *pos_in_read < config.m_read_position_min || *pos_in_read > config.m_read_position_max
    }
}

impl CheckFilter for vcf::m_highDp {
    fn check(config: &ThresholdParams, record: &mut vcf::Record) -> bool {
        *record.info.read_depth > config.m_max_coverage
    }
}
