use crate::utils::ByStrand;
use color_eyre::{Result, eyre::Context as _};
use rastair_vcf::{HeaderField, InfoField, InfoFieldNumber, VcfField};
use rust_htslib::bcf::Record;
use smallvec::SmallVec;
use std::ops::Deref;

/// Allele-specific RMS base quality by strand
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct StrandSpecificBaseQuality(pub SmallVec<ByStrand<f64>, 4>);

mod as_ss_bq {
    use super::*;
    use crate::utils::ByStrand;
    use rastair_vcf::StrandSpecificInfoField;

    impl Deref for StrandSpecificBaseQuality {
        type Target = SmallVec<ByStrand<f64>, 4>;

        fn deref(&self) -> &Self::Target {
            &self.0
        }
    }

    impl VcfField for StrandSpecificBaseQuality {
        const ID: &'static cstr8::CStr8 = cstr8::cstr8!("AS_SS_BQ");
    }

    impl HeaderField for StrandSpecificBaseQuality {
        const DESCRIPTION: &'static str = "Strand-specific RMS of base quality per allele";
    }

    impl InfoField for StrandSpecificBaseQuality {
        /// One tuple of two integers for each allele
        const NUMBER: InfoFieldNumber = InfoFieldNumber::OnePerAltAndRef;
        type Type = f32;

        fn write_header(header: &mut rust_htslib::bcf::Header) -> Result<()> {
            <Self as StrandSpecificInfoField>::write_header(header)
        }

        #[allow(clippy::cast_possible_truncation)]
        fn write(&self, record: &mut Record) -> Result<()> {
            {
                // Write OT field
                let tag = Self::ID_OT;
                let counts: SmallVec<f32, 8> = self.0.iter().map(|c| c.ot as f32).collect();
                record.push_info_float(tag, &counts).wrap_err("Failed to set AS_SS_BQ_OT field")?;
            }

            {
                // Write OB field
                let tag = Self::ID_OB;
                let counts: SmallVec<f32, 8> = self.0.iter().map(|c| c.ob as f32).collect();
                record.push_info_float(tag, &counts).wrap_err("Failed to set AS_SS_BQ_OB field")?;
            }
            Ok(())
        }

        fn description() -> Vec<rastair_vcf::reflect::Info> {
            Self::descriptions()
        }
    }

    impl StrandSpecificInfoField for StrandSpecificBaseQuality {
        const ID_OT: &'static cstr8::CStr8 = cstr8::cstr8!("AS_SS_BQ_OT");
        const DESCRIPTION_OT: &'static str =
            "Strand-specific RMS of base quality per allele on the original top strand";
        const ID_OB: &'static cstr8::CStr8 = cstr8::cstr8!("AS_SS_BQ_OB");
        const DESCRIPTION_OB: &'static str =
            "Strand-specific RMS of base quality per allele on the original bottom strand";
    }
}

/// Allele-specific RMS mapping quality by strand
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct StrandSpecificMappingQuality(pub SmallVec<ByStrand<f64>, 4>);

mod as_ss_mq {
    use super::*;

    impl Deref for StrandSpecificMappingQuality {
        type Target = SmallVec<ByStrand<f64>, 4>;

        fn deref(&self) -> &Self::Target {
            &self.0
        }
    }

    impl VcfField for StrandSpecificMappingQuality {
        const ID: &'static cstr8::CStr8 = cstr8::cstr8!("AS_SS_MQ");
    }

    impl HeaderField for StrandSpecificMappingQuality {
        const DESCRIPTION: &'static str = "Strand-specific RMS of mapping quality per allele (tuples of [reads_ot, reads_ob] for each allele)";
    }

    impl InfoField for StrandSpecificMappingQuality {
        /// One tuple of two integers for each allele
        const NUMBER: InfoFieldNumber = InfoFieldNumber::Dot;
        type Type = f32;

        #[allow(clippy::cast_possible_truncation)]
        fn write(&self, record: &mut Record) -> Result<()> {
            let counts: SmallVec<f32, 8> =
                self.0.iter().flat_map(|c| [c.ot as f32, c.ob as f32]).collect();
            record.push_info_float(Self::ID, &counts).wrap_err("Failed to set field")
        }
    }
}
