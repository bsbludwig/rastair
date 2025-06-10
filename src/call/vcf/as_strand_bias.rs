use color_eyre::eyre::{Context as _, Result};
use rastair2_vcf::{HeaderField, InfoField, InfoFieldNumber, VcfField};
use rust_htslib::bcf::Record;
use smallvec::SmallVec;

/// Allele-specific strand bias information for a variant
pub struct AlleleSpecificStrandBias(pub SmallVec<StrandCounts, 4>);

/// Counts of reads supporting an allele on the original top and bottom strands
pub struct StrandCounts {
    /// Number of reads from original top
    pub ot: u32,
    /// Number of reads from original bottom
    pub ob: u32,
}

impl VcfField for AlleleSpecificStrandBias {
    const ID: &'static str = "AS_SB";
}

impl HeaderField for AlleleSpecificStrandBias {
    const DESCRIPTION: &'static str =
        "Strand bias per allele: tuples of [reads_ot, reads_ob] for each allele";
}

impl InfoField for AlleleSpecificStrandBias {
    /// One tuple of two integers for each allele
    const NUMBER: InfoFieldNumber = InfoFieldNumber::Dot;
    type Type = u32;

    fn write(&self, record: &mut Record) -> Result<()> {
        let tag = Self::ID;
        record.clear_info_integer(tag.as_bytes()).wrap_err("Failed to clear field")?;
        let counts: SmallVec<i32, 8> = self
            .0
            .iter()
            .flat_map(|c| [c.ot, c.ob])
            .map(|count| i32::try_from(count).wrap_err("strand counts should fit in i32"))
            .collect::<Result<_>>()?;

        record.push_info_integer(tag.as_bytes(), &counts).wrap_err("Failed to set field")
    }
}
