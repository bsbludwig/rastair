use std::num::{NonZeroU32, ParseIntError};

use smol_str::SmolStr;

/// A struct representing a genomic region string.
///
/// # Examples
///
/// ```rust
/// use std::str::FromStr;
///
/// let entire_chr17 = RegionString::from_str("chr17").unwrap();
/// let chr17_from_100 = RegionString::from_str("chr17:100").unwrap();
/// let chr17_from_100_to_200 = RegionString::from_str("chr17:100-200").unwrap();
/// ```
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct RegionString {
    /// The chromosome name, includes the "chr" prefix.
    // todo: is there always a chr prefix? should we filter it out?
    pub chromosome: SmolStr,
    /// The start position of the region, inclusive.
    pub start: Option<NonZeroU32>,
    /// The end position of the region, inclusive.
    pub end: Option<NonZeroU32>,
}

impl std::str::FromStr for RegionString {
    type Err = RegionStringError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        if s.is_empty() {
            return Err(RegionStringError::EmptyInput);
        }
        if !s.is_ascii() {
            return Err(RegionStringError::InvalidAscii);
        }

        let mut parts = s.split(':');
        let chromosome = parts
            .next()
            .filter(|c| !c.trim().is_empty())
            .ok_or(RegionStringError::InvalidChromosome)?
            .into();

        let Some(range) = parts.next().filter(|r| !r.trim().is_empty()) else {
            return Ok(Self { chromosome, start: None, end: None });
        };

        let mut range_parts = range.split('-');
        let start = range_parts
            .next()
            .expect("split always returns given string")
            .parse()
            .map_err(|e| RegionStringError::InvalidStartPosition(e))?;

        let Some(end) = range_parts.next() else {
            return Ok(Self { chromosome, start: Some(start), end: None });
        };
        let end: NonZeroU32 = end.parse().map_err(|e| RegionStringError::InvalidEndPosition(e))?;
        if start.get() > end.get() {
            return Err(RegionStringError::StartGreaterThanEnd);
        }
        Ok(Self { chromosome, start: Some(start), end: Some(end) })
    }
}

#[derive(Debug, Clone, thiserror::Error)]
#[non_exhaustive]
pub enum RegionStringError {
    #[error("Empty region string")]
    EmptyInput,
    #[error("Invalid ASCII string")]
    InvalidAscii,
    #[error("Invalid chromosome name")]
    InvalidChromosome,
    #[error("Invalid start position")]
    InvalidStartPosition(ParseIntError),
    #[error("Invalid end position")]
    InvalidEndPosition(ParseIntError),
    #[error("Start position cannot be greater than end position")]
    StartGreaterThanEnd,
}

impl std::fmt::Display for RegionString {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.chromosome)?;
        if let Some(start) = self.start {
            write!(f, ":{}", start)?;
        }
        if let Some(end) = self.end {
            write!(f, "-{}", end)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn test_valid_region_strings() {
        // full chromosome
        let region = RegionString::from_str("chr1").unwrap();
        assert_eq!(region.chromosome, "chr1");
        assert_eq!(region.start, None);
        assert_eq!(region.end, None);

        // chromosome with start position
        let region = RegionString::from_str("chr2:100").unwrap();
        assert_eq!(region.chromosome, "chr2");
        assert_eq!(region.start, Some(NonZeroU32::new(100).unwrap()));
        assert_eq!(region.end, None);

        // chromosome with start and end positions
        let region = RegionString::from_str("chr3:100-200").unwrap();
        assert_eq!(region.chromosome, "chr3");
        assert_eq!(region.start, Some(NonZeroU32::new(100).unwrap()));
        assert_eq!(region.end, Some(NonZeroU32::new(200).unwrap()));
    }

    #[test]
    fn test_error_cases() {
        // empty input
        let err = RegionString::from_str("").unwrap_err();
        assert!(matches!(err, RegionStringError::EmptyInput));

        // invalid chromosome
        let err = RegionString::from_str(":100").unwrap_err();
        assert!(matches!(err, RegionStringError::InvalidChromosome));

        // invalid start position
        let err = RegionString::from_str("chr1:invalid").unwrap_err();
        assert!(matches!(err, RegionStringError::InvalidStartPosition(_)));

        // invalid end position
        let err = RegionString::from_str("chr1:100-invalid").unwrap_err();
        assert!(matches!(err, RegionStringError::InvalidEndPosition(_)));

        // start greater than end
        let err = RegionString::from_str("chr1:200-100").unwrap_err();
        assert!(matches!(err, RegionStringError::StartGreaterThanEnd));
    }

    #[test]
    fn test_edge_cases() {
        // with whitespace
        let region = RegionString::from_str(" chr4:150-250 ").unwrap();
        assert_eq!(region.chromosome, "chr4");
        assert_eq!(region.start, Some(NonZeroU32::new(150).unwrap()));
        assert_eq!(region.end, Some(NonZeroU32::new(250).unwrap()));

        // only whitespace in chromosome part
        let err = RegionString::from_str("  :100").unwrap_err();
        assert!(matches!(err, RegionStringError::InvalidChromosome));

        // empty range part
        let region = RegionString::from_str("chr1:  ").unwrap();
        assert_eq!(region.chromosome, "chr1");
        assert_eq!(region.start, None);
        assert_eq!(region.end, None);

        // invalid characters in start
        let err = RegionString::from_str("chr1:xxx").unwrap_err();
        assert!(matches!(err, RegionStringError::InvalidStartPosition(_)));

        // invalid characters in end
        let err = RegionString::from_str("chr1:100-xxx").unwrap_err();
        assert!(matches!(err, RegionStringError::InvalidEndPosition(_)));

        // non-ascii characters
        let err = RegionString::from_str("chrü1:100-200").unwrap_err();
        assert!(matches!(err, RegionStringError::InvalidAscii));
    }

    #[test]
    fn test_display() {
        // Test full chromosome
        let region = RegionString { chromosome: "chr1".into(), start: None, end: None };
        assert_eq!(region.to_string(), "chr1");

        // Test chromosome with start position
        let region = RegionString {
            chromosome: "chr2".into(),
            start: Some(NonZeroU32::new(100).unwrap()),
            end: None,
        };
        assert_eq!(region.to_string(), "chr2:100");

        // Test chromosome with start and end positions
        let region = RegionString {
            chromosome: "chr3".into(),
            start: Some(NonZeroU32::new(100).unwrap()),
            end: Some(NonZeroU32::new(200).unwrap()),
        };
        assert_eq!(region.to_string(), "chr3:100-200");
    }

    proptest::proptest! {
        #[test]
        fn proptest_roundtrip_region_string(
            // Generate chromosome names with "chr" prefix and some alphanumeric characters
            chrom in "chr[0-9A-Za-z]{1,10}",
            // Ensure start is between 1 and 1,000,000
            start in 1u32..1_000_000u32,
            // Ensure end is >= start and within reasonable bounds
            end_offset in 0u32..1_000_000u32
        ) {
            let end = start + end_offset;

            // Test chromosome only
            let region_str = chrom.clone();
            let parsed = RegionString::from_str(&region_str)?;
            assert_eq!(parsed.chromosome, chrom);
            assert_eq!(parsed.start, None);
            assert_eq!(parsed.end, None);
            assert_eq!(parsed.to_string(), region_str);

            // Test chromosome with start
            let region_str = format!("{}:{}", chrom, start);
            let parsed = RegionString::from_str(&region_str)?;
            assert_eq!(parsed.chromosome, chrom);
            assert_eq!(parsed.start, NonZeroU32::new(start));
            assert_eq!(parsed.end, None);
            assert_eq!(parsed.to_string(), region_str);

            // Test chromosome with start and end
            let region_str = format!("{}:{}-{}", chrom, start, end);
            let parsed = RegionString::from_str(&region_str)?;
            assert_eq!(parsed.chromosome, chrom);
            assert_eq!(parsed.start, NonZeroU32::new(start));
            assert_eq!(parsed.end, NonZeroU32::new(end));
            assert_eq!(parsed.to_string(), region_str);
        }

        #[test]
        fn proptest_roundtrip_random_string(
            // Generate random strings with alphanumeric characters
            random_str in r"\PC{0,100}"
        ) {
            let Ok(parsed) = RegionString::from_str(&random_str) else {
                // We're just checking that there is no panic, but errors are fine!
                return Ok(());
            };
            // Check that the parsed string is the same as the original
            assert_eq!(parsed.to_string(), random_str);
        }

    }
}
