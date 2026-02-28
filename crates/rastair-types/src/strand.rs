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

/// Determine strand from raw BAM flags.
///
/// # Flags used
///
/// |   Flag | Decimal | Meaning                              |
/// | -----: | ------: | ------------------------------------ |
/// | `0x10` |      16 | Read is mapped to the reverse strand |
/// | `0x20` |      32 | Mate is mapped to the reverse strand |
/// | `0x40` |      64 | Read is first in pair                |
/// | `0x80` |     128 | Read is second in pair               |
#[allow(clippy::collapsible_else_if)] // clearer
pub fn strand_from_flags(flags: u16) -> Strand {
    if flags & 0x1 == 0 {
        // Unpaired read: strand determined solely by alignment direction
        if flags & 0x10 != 0 { Strand::OB } else { Strand::OT }
    } else if flags & 0x40 != 0 {
        // First in pair
        if flags & 0x10 != 0 { Strand::OB } else { Strand::OT }
    } else if flags & 0x80 != 0 {
        // Last in pair
        if flags & 0x20 != 0 { Strand::OB } else { Strand::OT }
    } else {
        Strand::Unknown
    }
}

/// Extension trait to get strand information from a BCF record
pub trait StrandFromRecord {
    /// Create a `Strand` from a BCF record
    fn strand(&self) -> Strand;
}

impl StrandFromRecord for Record {
    fn strand(&self) -> Strand {
        strand_from_flags(self.flags())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_various_records() {
        let mut record = Record::default();
        record.set_flags(0x1 | 0x2 | 0x40 | 0x10); // First in pair, reverse strand
        assert_eq!(StrandFromRecord::strand(&record), Strand::OB);

        record.set_flags(0x1 | 0x2 | 0x80 | 0x20); // Second in pair, reverse strand
        assert_eq!(StrandFromRecord::strand(&record), Strand::OB);

        record.set_flags(0x1 | 0x2 | 0x40 | 0x20); // First in pair, mate reverse strand
        assert_eq!(StrandFromRecord::strand(&record), Strand::OT);

        record.set_flags(0x1 | 0x2 | 0x80 | 0x10); // Second in pair, mate reverse strand
        assert_eq!(StrandFromRecord::strand(&record), Strand::OT);

        record.set_flags(0x1 | 0x2 | 0x40); // First in pair, forward strand
        assert_eq!(StrandFromRecord::strand(&record), Strand::OT);

        record.set_flags(0x1 | 0x2 | 0x80); // Second in pair, forward strand
        assert_eq!(StrandFromRecord::strand(&record), Strand::OT);

        record.set_flags(0x00); // No flags set, ie top strand
        assert_eq!(StrandFromRecord::strand(&record), Strand::OT);

        record.set_flags(0x10); // No pairing flags, but read reverse strand -> OB
        assert_eq!(StrandFromRecord::strand(&record), Strand::OB);

        record.set_flags(0x01); // Paired but no first/second information
        assert_eq!(StrandFromRecord::strand(&record), Strand::Unknown);
    }

    #[test]
    fn test_unpaired_mode() {
        let mut record = Record::default();

        record.set_flags(0x00); // Single-end, forward
        assert_eq!(StrandFromRecord::strand(&record), Strand::OT);

        record.set_flags(0x10); // Single-end, reverse
        assert_eq!(StrandFromRecord::strand(&record), Strand::OB);

        record.set_flags(0x40 | 0x10); // Pair flags ignored in unpaired mode
        assert_eq!(StrandFromRecord::strand(&record), Strand::OB);

        record.set_flags(0x80 | 0x20); // Pair/mate flags ignored, only 0x10 matters
        assert_eq!(StrandFromRecord::strand(&record), Strand::OT);
    }

    #[test]
    fn test_unpaired_mode_ignores_paired_flags() {
        let mut record = Record::default();
        record.set_flags(0x01 | 0x02 | 0x20 | 0x40);
        assert_eq!(StrandFromRecord::strand(&record), Strand::OT);

        record.set_flags(0x01 | 0x02 | 0x20 | 0x40 | 0x10);
        assert_eq!(StrandFromRecord::strand(&record), Strand::OB);
    }
}
