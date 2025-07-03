//! Mask specific parts of reads at the edges (mbias)
//!
//! There is a bias towards low-quality methylation calls caused by evidence
//! that comes from the fringes of the read. You can mask out specific areas
//! from reads to reduce this.

use crate::{call::variants::SeenBase, utils::Strand};
use smallvec::SmallVec;
use std::{num::ParseIntError, str::FromStr};

#[derive(Debug, Clone, Default, clap::Args)]
pub struct ReadMaskParams {
    /// For OT reads, exclude `[r1_start, r1_end, r2_start, r2_end]` bases from counting.
    ///
    /// The coordinates are relative to the read, so start is the distance from
    /// the 5' of the read, the end is the distance to the 3', irrespective of
    /// which way around the read aligns to the reference.
    ///
    /// Also note that the distance is relative to read length, not alignment
    /// length, so soft-clipped bases count, too!
    #[arg(long = "nOT", default_value = "0,0,0,0")]
    n_ot: ReadMaskSetting,

    /// For OB reads, exclude `[r1_start, r1_end, r2_start, r2_end]` bases from counting.
    ///
    /// The coordinates are relative to the read, so start is the distance from
    /// the 5' of the read, the end is the distance to the 3', irrespective of
    /// which way around the read aligns to the reference.
    ///
    /// Also note that the distance is relative to read length, not alignment
    /// length, so soft-clipped bases count, too!
    #[arg(long = "nOB", default_value = "0,0,0,0")]
    n_ob: ReadMaskSetting,
}

impl ReadMaskParams {
    pub fn filter(&self, read: &SeenBase) -> bool {
        let len = read.position.read_length;
        let pos = read.position.pos;

        match (read.strand, read.reverse) {
            (Strand::OT, true) => {
                let mask = self.n_ot.r2;

                let too_small = len < mask.from_start + mask.from_end + 1;

                // flipped end/start mask, the read is mapped in reverse
                let masked_start = pos < mask.from_start;
                let masked_end = pos > len - mask.from_start - 1;

                !too_small && !masked_start && !masked_end
            }
            (Strand::OT, false) => {
                let mask = self.n_ot.r1;

                let too_small = len < mask.from_start + mask.from_end + 1;

                // normal start/end mask, the read is mapped in forward
                let masked_start = pos < mask.from_start;
                let masked_end = pos > len - mask.from_end - 1;

                !too_small && !masked_start && !masked_end
            }
            (Strand::OB, true) => {
                let mask = self.n_ob.r1;

                let too_small = len < mask.from_start + mask.from_end + 1;

                // I'm flipping the start/end here, because the R1 of the OB is reversed but
                // samtools reports it in ref direction, so if I want to remove 5 bases from the start
                // of the read, that's actually the "end" in the coordinate system that htslib provides
                let masked_start = pos < mask.from_end;
                let masked_end = pos > len - mask.from_start - 1;

                !too_small && !masked_start && !masked_end
            }
            (Strand::OB, false) => {
                let mask = self.n_ob.r2;

                let too_small = len < mask.from_start + mask.from_end + 1;

                // normal start/end mask, the read is mapped in forward
                let masked_start = pos < mask.from_start;
                let masked_end = pos > len - mask.from_end - 1;

                !too_small && !masked_start && !masked_end
            }
        }
    }
}

/// Represent a read softmask, to exclude certain portions of the read
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ReadMask {
    /// Number of bases to exclude from the start of the read
    pub from_start: u32,
    /// Number of bases to exclude from the end of the read
    pub from_end: u32,
}

/// Mask settings for both read orientations
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ReadMaskSetting {
    /// Mask for reads in forward orientation
    pub r1: ReadMask,
    /// Mask for reads in reverse orientation
    pub r2: ReadMask,
}

#[derive(Debug, PartialEq, Eq, thiserror::Error)]
pub enum ParseMaskError {
    #[error("Invalid read mask format, expected four comma-separated integers")]
    InvalidFormat,
    #[error("Failed to parse read mask integer")]
    ParseIntError(#[source] ParseIntError),
}

impl FromStr for ReadMaskSetting {
    type Err = ParseMaskError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let values = s
            .trim()
            .splitn(5, ',') // if there are 5 values, it's an error anyway
            .map(|i| i.parse::<u32>().map_err(ParseMaskError::ParseIntError))
            .collect::<Result<SmallVec<u32, 5>, _>>()?;

        match values[..] {
            [r1_left, r1_right, r2_left, r2_right] => Ok(ReadMaskSetting {
                r1: ReadMask { from_start: r1_left, from_end: r1_right },
                r2: ReadMask { from_start: r2_left, from_end: r2_right },
            }),
            _ => Err(ParseMaskError::InvalidFormat),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use proptest::prelude::*;

    proptest! {
        #[test]
        fn test_parse_read_mask_setting_valid(
            r1_left in 0u32..1000,
            r1_right in 0u32..1000,
            r2_left in 0u32..1000,
            r2_right in 0u32..1000
        ) {
            let input = format!("{r1_left},{r1_right},{r2_left},{r2_right}");
            let result = ReadMaskSetting::from_str(&input).unwrap();

            assert_eq!(result.r1.from_start, r1_left);
            assert_eq!(result.r1.from_end, r1_right);
            assert_eq!(result.r2.from_start, r2_left);
            assert_eq!(result.r2.from_end, r2_right);
        }

        #[test]
        fn test_parse_read_mask_setting_invalid_count(
            values in prop::collection::vec(0u32..1000, 0..10)
        ) {
            prop_assume!(values.len() != 4);

            let input = values.iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(",");

            let result = ReadMaskSetting::from_str(&input);
            assert!(result.is_err());
        }

        #[test]
        fn test_parse_read_mask_setting_invalid_integer(
            invalid_str in "[a-zA-Z!@#$%^&*()]+",
            valid_nums in prop::collection::vec(0u32..1000, 0..3)
        ) {
            let mut parts = valid_nums.iter().map(|v| v.to_string()).collect::<Vec<_>>();
            parts.push(invalid_str);

            // Pad with more valid numbers if needed to get exactly 4 elements
            while parts.len() < 4 {
                parts.push("0".to_string());
            }
            parts.truncate(4);

            let input = parts.join(",");
            let result = ReadMaskSetting::from_str(&input);
            assert!(result.is_err());
            assert!(matches!(result.unwrap_err(), ParseMaskError::ParseIntError(_)));
        }
    }

    #[test]
    fn test_parse_read_mask_setting_basic_cases() {
        // Test default case
        let result = ReadMaskSetting::from_str("0,0,0,0").unwrap();
        assert_eq!(result, ReadMaskSetting::default());

        // Test non-zero values
        let result = ReadMaskSetting::from_str("1,2,3,4").unwrap();
        assert_eq!(result.r1.from_start, 1);
        assert_eq!(result.r1.from_end, 2);
        assert_eq!(result.r2.from_start, 3);
        assert_eq!(result.r2.from_end, 4);

        // Test error cases
        assert!(ReadMaskSetting::from_str("1,2,3").is_err());
        assert!(ReadMaskSetting::from_str("1,2,3,4,5").is_err());
        assert!(ReadMaskSetting::from_str("1,2,x,4").is_err());
        assert!(ReadMaskSetting::from_str("").is_err());
    }
}
