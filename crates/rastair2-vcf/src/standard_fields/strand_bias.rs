use crate::{
    FormatField, FormatFieldNumber, HeaderField, InfoField, InfoFieldNumber, InfoFieldValue as _,
    VcfField,
};
use color_eyre::eyre::Context as _;

/// Strand bias information for a variant
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct StrandBias {
    /// Number of reads supporting the reference allele on the forward strand
    pub reads_ref_fwd: usize,
    /// Number of reads supporting the reference allele on the reverse strand
    pub reads_ref_rev: usize,
    /// Number of reads supporting the alternative allele on the forward strand
    pub reads_alt_fwd: usize,
    /// Number of reads supporting the alternative allele on the reverse strand
    pub reads_alt_rev: usize,
}

impl VcfField for StrandBias {
    const ID: &'static str = "SB";
}

impl HeaderField for StrandBias {
    const DESCRIPTION: &'static str =
        "Strand bias: counts of [reads_ref_fwd, reads_ref_rev, reads_alt_fwd, reads_alt_rev]";
}

impl InfoField for StrandBias {
    const NUMBER: InfoFieldNumber = InfoFieldNumber::Num(4);
    type Type = usize;

    fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::eyre::Result<()> {
        let tag = Self::ID;
        record.clear_info_integer(tag.as_bytes()).wrap_err_with(|| {
            format!("Failed to clear info field {tag} ({})", Self::Type::TYPE_NAME)
        })?;
        record
            .push_info_integer(
                tag.as_bytes(),
                &[
                    self.reads_ref_fwd as i32,
                    self.reads_ref_rev as i32,
                    self.reads_alt_fwd as i32,
                    self.reads_alt_rev as i32,
                ],
            )
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::Type::TYPE_NAME))
    }
}

impl FormatField for StrandBias {
    const NUMBER: FormatFieldNumber = FormatFieldNumber::Num(4);
    type Type = usize;

    fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::eyre::Result<()> {
        let tag = Self::ID;
        record
            .push_format_integer(
                tag.as_bytes(),
                &[
                    self.reads_ref_fwd as i32,
                    self.reads_ref_rev as i32,
                    self.reads_alt_fwd as i32,
                    self.reads_alt_rev as i32,
                ],
            )
            .wrap_err_with(|| format!("Failed to set info field {tag} ({})", Self::Type::TYPE_NAME))
    }
}
