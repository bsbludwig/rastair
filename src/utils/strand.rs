use std::fmt;

use color_eyre::eyre::{Result, bail};
use rust_htslib::bam::Record;

/// Original top or bottom strand of a read
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    /// Original top
    OT,
    /// Original bottom
    OB,
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::OT => write!(f, "OT"),
            Strand::OB => write!(f, "OB"),
        }
    }
}

pub trait StrandFromRecord {
    /// Create a `Strand` from a BCF record
    fn strand(&self) -> Result<Strand>;
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
    fn strand(&self) -> Result<Strand> {
        if self.is_first_in_template() {
            if self.is_reverse() {
                Ok(Strand::OB) // Original bottom
            } else {
                Ok(Strand::OT) // Original top
            }
        } else if self.is_last_in_template() {
            if self.is_mate_reverse() {
                Ok(Strand::OB) // Original bottom
            } else {
                Ok(Strand::OT) // Original top
            }
        } else {
            bail!(
                "Record {} is not first or last in template, cannot determine strand from flags {}",
                String::from_utf8_lossy(self.qname()),
                self.flags()
            );
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
        assert_eq!(StrandFromRecord::strand(&record).unwrap(), Strand::OB);

        record.set_flags(0x80 | 0x20); // Second in pair, reverse strand
        assert_eq!(StrandFromRecord::strand(&record).unwrap(), Strand::OB);

        record.set_flags(0x40 | 0x20); // First in pair, mate reverse strand
        assert_eq!(StrandFromRecord::strand(&record).unwrap(), Strand::OT);

        record.set_flags(0x80 | 0x10); // Second in pair, mate reverse strand
        assert_eq!(StrandFromRecord::strand(&record).unwrap(), Strand::OT);

        record.set_flags(0x40); // First in pair, forward strand
        assert_eq!(StrandFromRecord::strand(&record).unwrap(), Strand::OT);

        record.set_flags(0x80); // Second in pair, forward strand
        assert_eq!(StrandFromRecord::strand(&record).unwrap(), Strand::OT);

        record.set_flags(0x00); // No flags set
        assert!(StrandFromRecord::strand(&record).is_err(), "Should fail with no flags set");
    }
}
