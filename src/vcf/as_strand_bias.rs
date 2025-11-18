use crate::utils::ByStrand;
use color_eyre::eyre::{Context as _, Result};
use rastair_types::SmallVec;
use rastair_vcf::{HeaderField, InfoField, InfoFieldNumber, StrandSpecificInfoField, VcfField};
use rust_htslib::bcf::Record;
use std::ops::Deref;

/// Allele-specific strand bias information for a variant
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct AlleleSpecificStrandBias(pub SmallVec<ByStrand<u32>, 4>);

impl Deref for AlleleSpecificStrandBias {
    type Target = SmallVec<ByStrand<u32>, 4>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl VcfField for AlleleSpecificStrandBias {
    const ID: &'static cstr8::CStr8 = cstr8::cstr8!("AS_SB");
}

impl HeaderField for AlleleSpecificStrandBias {
    const DESCRIPTION: &'static str = "strand bias per allele";
}

impl StrandSpecificInfoField for AlleleSpecificStrandBias {
    const ID_OT: &'static cstr8::CStr8 = cstr8::cstr8!("AS_SB_OT");
    const DESCRIPTION_OT: &'static str = "OT counts per allele";

    const ID_OB: &'static cstr8::CStr8 = cstr8::cstr8!("AS_SB_OB");
    const DESCRIPTION_OB: &'static str = "OB counts per allele";
}

/// This writes two fields to the VCF INFO section, one for OB strand, one for OT
impl InfoField for AlleleSpecificStrandBias {
    const NUMBER: InfoFieldNumber = InfoFieldNumber::OnePerAltAndRef;
    type Type = u32;

    fn write_header(header: &mut rust_htslib::bcf::Header) -> Result<()> {
        <Self as StrandSpecificInfoField>::write_header(header)
    }

    fn write(&self, record: &mut Record) -> Result<()> {
        {
            // Write OT field
            let tag = Self::ID_OT;
            let counts: SmallVec<i32, 8> = self
                .0
                .iter()
                .map(|c| c.ot)
                .map(i32::try_from)
                .collect::<Result<_, _>>()
                .wrap_err("strand counts should fit in i32")?;
            record.push_info_integer(tag, &counts).wrap_err("Failed to set AS_SB_OT field")?;
        }

        {
            // Write OB field
            let tag = Self::ID_OB;
            let counts: SmallVec<i32, 8> = self
                .0
                .iter()
                .map(|c| c.ob)
                .map(i32::try_from)
                .collect::<Result<_, _>>()
                .wrap_err("strand counts should fit in i32")?;
            record.push_info_integer(tag, &counts).wrap_err("Failed to set AS_SB_OB field")?;
        }

        Ok(())
    }

    fn description() -> Vec<rastair_vcf::reflect::Info> {
        Self::descriptions()
    }
}
