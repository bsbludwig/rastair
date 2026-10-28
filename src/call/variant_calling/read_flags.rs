use crate::utils::cli;
use better_default::Default;
use clap_num::maybe_hex;

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

/// Whether the library's reads come in pairs (`--unpaired` says they do not).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Pairing {
    #[default]
    Paired,
    Unpaired,
}

impl Pairing {
    #[must_use]
    pub const fn from_unpaired_flag(unpaired: bool) -> Self {
        if unpaired { Self::Unpaired } else { Self::Paired }
    }
}

impl ReadFlags {
    /// Check if raw flags match the include flags and do not match the exclude flags.
    /// Unpaired, pairing-related flags are masked out before comparison.
    pub fn filter_flags(&self, flags: u16, pairing: Pairing) -> bool {
        let (flags, include_flags, exclude_flags) = match pairing {
            Pairing::Unpaired => (
                flags & !UNPAIRED_IGNORED_FLAGS,
                self.include_flags & !UNPAIRED_IGNORED_FLAGS,
                self.exclude_flags & !UNPAIRED_IGNORED_FLAGS,
            ),
            Pairing::Paired => (flags, self.include_flags, self.exclude_flags),
        };
        let include = flags & include_flags == include_flags;
        let exclude = flags & exclude_flags != 0;
        include && !exclude
    }
}

#[allow(unused, reason = "list all for documentation purposes")]
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

    pub const UNPAIRED_IGNORED_FLAGS: u16 = IS_PAIRED
        | IS_PROPERLY_PAIRED
        | MATE_IS_REVERSE_STRAND
        | IS_FIRST_IN_PAIR
        | IS_SECOND_IN_PAIR;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_filter_rejects_single_end_reads() {
        let read_flags = ReadFlags::default();
        assert!(!read_flags.filter_flags(0x0, Pairing::Paired)); // single-end forward
        assert!(!read_flags.filter_flags(0x10, Pairing::Paired)); // single-end reverse
    }

    #[test]
    fn unpaired_mode_accepts_single_end_reads() {
        let read_flags = ReadFlags::default();
        assert!(read_flags.filter_flags(0x0, Pairing::Unpaired)); // single-end forward
        assert!(read_flags.filter_flags(0x10, Pairing::Unpaired)); // single-end reverse
    }

    #[test]
    fn default_filter_keeps_paired_requirements_for_paired_reads() {
        let read_flags = ReadFlags::default();
        assert!(!read_flags.filter_flags(IS_PAIRED, Pairing::Paired)); // missing proper pair
        assert!(read_flags.filter_flags(IS_PAIRED | IS_PROPERLY_PAIRED, Pairing::Paired));
    }

    #[test]
    fn unpaired_mode_ignores_pair_related_flags() {
        let read_flags = ReadFlags::default();

        // In unpaired mode, required 0x1/0x2 are ignored
        assert!(read_flags.filter_flags(IS_PAIRED, Pairing::Unpaired));

        // 0x20 and 0x40 are ignored as well
        assert!(read_flags.filter_flags(
            IS_PAIRED | IS_PROPERLY_PAIRED | MATE_IS_REVERSE_STRAND | IS_FIRST_IN_PAIR,
            Pairing::Unpaired
        ));
    }

    #[test]
    fn exclude_flags_still_apply_to_single_end_reads() {
        let read_flags = ReadFlags::default();
        assert!(!read_flags.filter_flags(IS_UNMAPPED, Pairing::Unpaired));
    }
}
