use crate::utils::cli;
use better_default::Default;
use clap_num::maybe_hex;
use rust_htslib::bam::Record;

use flags::*;

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct ReadFlags {
    /// Include reads that match all of these bit-flags
    #[arg(
        short = 'f', long,
        value_parser=maybe_hex::<u16>,
        default_value_t = IS_PAIRED | IS_PROPERLY_PAIRED
    )]
    #[default(IS_PAIRED | IS_PROPERLY_PAIRED)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub include_flags: u16,

    /// Exclude reads that match any of these bit-flags
    #[arg(
        short = 'F', long,
        value_parser=maybe_hex::<u16>,
        default_value_t = IS_FAILED | IS_NOT_PRIMARY | IS_UNMAPPED | MATE_IS_UNMAPPED | IS_DUPLICATE | IS_SUPPLEMENTAL
    )]
    #[default(IS_FAILED | IS_NOT_PRIMARY | IS_UNMAPPED | MATE_IS_UNMAPPED | IS_DUPLICATE | IS_SUPPLEMENTAL)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub exclude_flags: u16,
}

impl ReadFlags {
    /// Check if the read matches the include flags and does not match the exclude flags
    pub fn filter(&self, record: &Record) -> bool {
        let flags = record.flags();
        let include = flags & self.include_flags == self.include_flags;
        let exclude = flags & self.exclude_flags != 0;
        include && !exclude
    }
}

#[allow(unused)]
mod flags {
    pub const IS_PAIRED: u16 = 0x1;
    pub const IS_PROPERLY_PAIRED: u16 = 0x2;
    pub const IS_UNMAPPED: u16 = 0x4;
    pub const MATE_IS_UNMAPPED: u16 = 0x8;
    pub const IS_REVERSE_STRAND: u16 = 0x10;
    pub const MATE_IS_REVERSE_STRAND: u16 = 0x20;
    pub const IS_FIRST_IN_PAIR: u16 = 0x40;
    pub const IS_SECOND_IN_PAIR: u16 = 0x80;
    pub const IS_NOT_PRIMARY: u16 = 0x100;
    pub const IS_FAILED: u16 = 0x200;
    pub const IS_DUPLICATE: u16 = 0x400;
    pub const IS_SUPPLEMENTAL: u16 = 0x800;
}
