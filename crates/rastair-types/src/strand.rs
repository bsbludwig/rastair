use rust_htslib::bam::Record;
use std::fmt;

/// Original top or bottom strand of a read
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[must_use]
pub enum Strand {
    /// Original top
    OT,
    /// Original bottom
    OB,
    /// Unknown
    Unknown,
}

impl Strand {
    /// Return `Some(self)` if strand is known, otherwise `None`
    pub fn ok(self) -> Option<Self> {
        match self {
            Strand::OT | Strand::OB => Some(self),
            Strand::Unknown => None,
        }
    }

    /// As symbol `+` (OT) or `-` (OB), or `.` (unknown)
    pub fn as_symbol(&self) -> &'static str {
        match self {
            Strand::OT => "+",
            Strand::OB => "-",
            Strand::Unknown => ".",
        }
    }
}

impl AsRef<str> for Strand {
    fn as_ref(&self) -> &str {
        match self {
            Strand::OT => "OT",
            Strand::OB => "OB",
            Strand::Unknown => "NA",
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_ref())
    }
}

/// Extension trait to get strand information from a BCF record
pub trait StrandFromRecord {
    /// Create a `Strand` from a BCF record
    fn strand(&self) -> Strand;

    /// Create a `Strand` from a record, optionally allowing single-end mode.
    fn strand_with_single_strand_mode(&self, single_strand_mode: bool) -> Strand;
}

/// Get strand from `htslib` record
///
/// # Flags used
///
/// |   Flag | Decimal | Meaning                              | `htslib` method        |
/// | -----: | ------: | ------------------------------------ | ---------------------- |
/// | `0x10` |      16 | Read is mapped to the reverse strand | `is_reverse`           |
/// | `0x20` |      32 | Mate is mapped to the reverse strand | `is_mate_reverse`      |
/// | `0x40` |      64 | Read is first in pair                | `is_first_in_template` |
/// | `0x80` |     128 | Read is second in pair               | `is_last_in_template`  |
impl StrandFromRecord for Record {
    /// Get strand from record
    #[allow(clippy::collapsible_else_if)] // clearer
    fn strand(&self) -> Strand {
        if self.is_first_in_template() {
            if self.is_reverse() {
                Strand::OB // Original bottom
            } else {
                Strand::OT // Original top
            }
        } else if self.is_last_in_template() {
            if self.is_mate_reverse() {
                Strand::OB // Original bottom
            } else {
                Strand::OT // Original top
            }
        } else {
            Strand::Unknown // Not a paired read or no flags set
        }
    }

    fn strand_with_single_strand_mode(&self, single_strand_mode: bool) -> Strand {
        if single_strand_mode {
            if self.is_reverse() { Strand::OB } else { Strand::OT }
        } else {
            StrandFromRecord::strand(self)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_various_records() {
        let mut record = Record::default();
        record.set_flags(0x40 | 0x10); // First in pair, reverse strand
        assert_eq!(StrandFromRecord::strand(&record), Strand::OB);

        record.set_flags(0x80 | 0x20); // Second in pair, reverse strand
        assert_eq!(StrandFromRecord::strand(&record), Strand::OB);

        record.set_flags(0x40 | 0x20); // First in pair, mate reverse strand
        assert_eq!(StrandFromRecord::strand(&record), Strand::OT);

        record.set_flags(0x80 | 0x10); // Second in pair, mate reverse strand
        assert_eq!(StrandFromRecord::strand(&record), Strand::OT);

        record.set_flags(0x40); // First in pair, forward strand
        assert_eq!(StrandFromRecord::strand(&record), Strand::OT);

        record.set_flags(0x80); // Second in pair, forward strand
        assert_eq!(StrandFromRecord::strand(&record), Strand::OT);

        record.set_flags(0x00); // No flags set
        assert_eq!(StrandFromRecord::strand(&record), Strand::Unknown);

        record.set_flags(0x10); // No pairing flags
        assert_eq!(StrandFromRecord::strand(&record), Strand::Unknown);

        record.set_flags(0x01); // Paired but no first/second information
        assert_eq!(StrandFromRecord::strand(&record), Strand::Unknown);
    }

    #[test]
    fn test_single_strand_mode() {
        let mut record = Record::default();

        record.set_flags(0x00); // Single-end, forward
        assert_eq!(StrandFromRecord::strand_with_single_strand_mode(&record, true), Strand::OT);

        record.set_flags(0x10); // Single-end, reverse
        assert_eq!(StrandFromRecord::strand_with_single_strand_mode(&record, true), Strand::OB);

        record.set_flags(0x40 | 0x10); // Pair flags ignored in single-strand mode
        assert_eq!(StrandFromRecord::strand_with_single_strand_mode(&record, true), Strand::OB);

        record.set_flags(0x80 | 0x20); // Pair/mate flags ignored, only 0x10 matters
        assert_eq!(StrandFromRecord::strand_with_single_strand_mode(&record, true), Strand::OT);

        record.set_flags(0x00); // Without single-strand mode, still unknown
        assert_eq!(
            StrandFromRecord::strand_with_single_strand_mode(&record, false),
            Strand::Unknown
        );
    }

    #[test]
    fn test_single_strand_mode_ignores_paired_flags() {
        let mut record = Record::default();
        record.set_flags(0x01 | 0x02 | 0x20 | 0x40);
        assert_eq!(StrandFromRecord::strand_with_single_strand_mode(&record, true), Strand::OT);

        record.set_flags(0x01 | 0x02 | 0x20 | 0x40 | 0x10);
        assert_eq!(StrandFromRecord::strand_with_single_strand_mode(&record, true), Strand::OB);
    }
}
