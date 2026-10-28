//! What rastair makes of a MAPQ of 255, "mapping quality not available" in
//! the SAM spec. Code that needs a read's MAPQ as a number, or filters on
//! it, goes through [`MappingQualityExt`], so a 255 is never taken as a
//! near-certain mapping.

use seqair_types::MappingQuality;

/// What an aligner such as bwa writes for a confidently unique read.
pub const MAPQ_UNIQUE: u8 = 60;

pub trait MappingQualityExt {
    /// The MAPQ, an unavailable one read as a unique mapping. MAPQ 0 would be
    /// the other reading, but an aligner that writes 255 throughout would then
    /// leave rastair no read it trusts.
    fn or_unique(self) -> u8;

    /// Whether the read passes a minimum-MAPQ filter. An unavailable MAPQ
    /// passes whatever the minimum, as with `samtools view -q`, and a minimum
    /// of 255 is above every available one.
    fn at_least(self, min: MappingQuality) -> bool;
}

impl MappingQualityExt for MappingQuality {
    fn or_unique(self) -> u8 {
        self.get().unwrap_or(MAPQ_UNIQUE)
    }

    fn at_least(self, min: MappingQuality) -> bool {
        self.get().is_none_or(|mapq| mapq >= min.as_byte())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_unavailable_mapq_is_a_unique_mapping_that_passes_any_filter() {
        assert_eq!(MappingQuality::UNAVAILABLE.or_unique(), MAPQ_UNIQUE);
        assert!(MappingQuality::UNAVAILABLE.at_least(MappingQuality::UNAVAILABLE));
        assert!(!MappingQuality::from_byte(254).at_least(MappingQuality::UNAVAILABLE));
        assert_eq!(MappingQuality::from_byte(254).or_unique(), 254);
        assert_eq!(MappingQuality::from_byte(0).or_unique(), 0);
        assert!(MappingQuality::from_byte(20).at_least(MappingQuality::from_byte(20)));
        assert!(!MappingQuality::from_byte(19).at_least(MappingQuality::from_byte(20)));
    }
}
