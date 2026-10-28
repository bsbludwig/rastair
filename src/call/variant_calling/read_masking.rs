//! Mask specific parts of reads at the edges (mbias)
//!
//! There is a bias towards low-quality methylation calls caused by evidence
//! that comes from the fringes of the read. You can mask out specific areas
//! from reads to reduce this.

#[cfg(test)]
use crate::call::tests::pileup::SimpleRead;
use crate::utils::{Strand, cli};
use seqair_types::{BamFlags, SmallVec};
use std::{num::ParseIntError, ops::Range, str::FromStr};

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
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
    #[arg(help_heading = cli::sections::FILTER)]
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
    #[arg(help_heading = cli::sections::FILTER)]
    n_ob: ReadMaskSetting,
}

/// Which way a read aligns to the reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ReadDirection {
    #[default]
    Forward,
    Reverse,
}

impl ReadDirection {
    #[must_use]
    pub const fn of(flags: BamFlags) -> Self {
        if flags.is_reverse() { Self::Reverse } else { Self::Forward }
    }
}

impl ReadMaskParams {
    pub fn new(n_ot: ReadMaskSetting, n_ob: ReadMaskSetting) -> Self {
        Self { n_ot, n_ob }
    }

    #[cfg(test)]
    pub fn filter(&self, read: &SimpleRead) -> bool {
        self.filter_fields(
            read.strand,
            read.direction,
            read.position.pos,
            read.position.read_length,
        )
    }

    pub fn filter_fields(
        &self,
        strand: Strand,
        direction: ReadDirection,
        pos: u32,
        len: u32,
    ) -> bool {
        self.keep_range(strand, direction, len).is_some_and(|keep| keep.contains(&pos))
    }

    /// The `[start, end)` window of read positions this read contributes, or
    /// `None` when the mask covers the whole read (including reads too short
    /// for it, which contribute nothing anywhere).
    ///
    /// Computed once per read at fetch time so the per-column check is a range
    /// comparison instead of four branches and two subtractions per base.
    pub fn keep_range(
        &self,
        strand: Strand,
        direction: ReadDirection,
        len: u32,
    ) -> Option<Range<u32>> {
        // The mask is expressed in read 5'->3' coordinates while `pos` counts
        // in reference direction, so a reverse-aligned mate has its two ends
        // swapped. R1/R2 selection differs between OT and OB for the same
        // reason.
        let mask = match (strand, direction) {
            (Strand::OT, ReadDirection::Reverse) => self.n_ot.r2,
            (Strand::OT, ReadDirection::Forward) => self.n_ot.r1,
            (Strand::OB, ReadDirection::Reverse) => self.n_ob.r1,
            (Strand::OB, ReadDirection::Forward) => self.n_ob.r2,
        };

        let (from_start, from_end) = match direction {
            ReadDirection::Reverse => (mask.from_end, mask.from_start),
            ReadDirection::Forward => (mask.from_start, mask.from_end),
        };

        let keep = len.checked_sub(from_start)?.checked_sub(from_end)?;
        if keep == 0 {
            return None;
        }
        Some(from_start..from_start.saturating_add(keep))
    }
}

/// Represent a read softmask, to exclude certain portions of the read
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ReadMask {
    /// Number of bases to exclude from the start of the read
    pub from_start: u32,
    /// Number of bases to exclude from the end of the read
    pub from_end: u32,
}

/// Mask settings for both mates of a pair
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ReadMaskSetting {
    /// Mask for read 1: the mate aligned forward on OT, reverse on OB
    pub r1: ReadMask,
    /// Mask for read 2: the mate aligned reverse on OT, forward on OB
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
    use crate::utils::default;

    use super::*;

    use hegel::{TestCase, generators as gs};

    /// The masking rule written out position by position, in i64 so the
    /// subtractions cannot wrap: an independent oracle for `keep_range`.
    fn masked_oracle(
        p: &ReadMaskParams,
        strand: Strand,
        direction: ReadDirection,
        pos: u32,
        len: u32,
    ) -> bool {
        let reverse = direction == ReadDirection::Reverse;
        let mask = match (strand, reverse) {
            (Strand::OT, true) => p.n_ot.r2,
            (Strand::OT, false) => p.n_ot.r1,
            (Strand::OB, true) => p.n_ob.r1,
            (Strand::OB, false) => p.n_ob.r2,
        };
        // Reverse-aligned mates have the 5'/3' ends swapped in reference coords.
        let (start, end) = if reverse {
            (i64::from(mask.from_end), i64::from(mask.from_start))
        } else {
            (i64::from(mask.from_start), i64::from(mask.from_end))
        };
        let (pos, len) = (i64::from(pos), i64::from(len));
        let too_small = len < start + end + 1;
        !too_small && pos >= start && pos < len - end
    }

    fn params(ot: [u32; 4], ob: [u32; 4]) -> ReadMaskParams {
        let setting = |v: [u32; 4]| ReadMaskSetting {
            r1: ReadMask { from_start: v[0], from_end: v[1] },
            r2: ReadMask { from_start: v[2], from_end: v[3] },
        };
        ReadMaskParams::new(setting(ot), setting(ob))
    }

    #[test]
    fn keep_range_covers_the_whole_read_when_unmasked() {
        let p = params([0, 0, 0, 0], [0, 0, 0, 0]);
        assert_eq!(p.keep_range(Strand::OT, ReadDirection::Forward, 100), Some(0..100));
        assert_eq!(p.keep_range(Strand::OB, ReadDirection::Reverse, 100), Some(0..100));
    }

    #[test]
    fn keep_range_trims_both_ends_and_flips_for_reverse_mates() {
        // OT R1 masks 5 from the 5' end and 3 from the 3' end; a reverse-aligned
        // OT mate uses R2 with the two ends swapped.
        let p = params([5, 3, 7, 2], [0, 0, 0, 0]);
        assert_eq!(p.keep_range(Strand::OT, ReadDirection::Forward, 100), Some(5..97));
        assert_eq!(p.keep_range(Strand::OT, ReadDirection::Reverse, 100), Some(2..93));
    }

    #[test]
    fn a_read_shorter_than_its_mask_keeps_nothing() {
        let p = params([50, 50, 50, 50], [50, 50, 50, 50]);
        assert_eq!(p.keep_range(Strand::OT, ReadDirection::Forward, 100), None);
        assert_eq!(p.keep_range(Strand::OT, ReadDirection::Forward, 99), None);
        assert_eq!(p.keep_range(Strand::OT, ReadDirection::Forward, 101), Some(50..51));
    }

    fn mask_lengths() -> impl gs::PrintableGenerator<[u32; 4]> {
        gs::arrays(gs::integers::<u32>().max_value(39))
    }

    /// `keep_range` must decide exactly what the original per-position
    /// masking rule decided, for every position of every read.
    #[hegel::test]
    #[hegel::explicit_test_case(
        ot = [0u32; 4],
        ob = [0u32; 4],
        len = 0u32,
        direction = ReadDirection::Forward,
        strand = Strand::OB
    )]
    fn keep_range_agrees_with_the_original_rule(tc: TestCase) {
        let ot = tc.draw(mask_lengths());
        let ob = tc.draw(mask_lengths());
        let len = tc.draw(gs::integers::<u32>().max_value(79));
        let direction =
            tc.draw(gs::sampled_from([ReadDirection::Forward, ReadDirection::Reverse].as_slice()));
        let strand = tc.draw(gs::sampled_from([Strand::OT, Strand::OB].as_slice()));
        let p = params(ot, ob);
        for pos in 0..len.max(1) {
            assert_eq!(
                p.filter_fields(strand, direction, pos, len),
                masked_oracle(&p, strand, direction, pos, len),
                "pos {pos} of {len} ({strand:?}, {direction:?})"
            );
        }
    }

    #[hegel::test]
    fn test_parse_read_mask_setting_valid(tc: TestCase) {
        let [r1_left, r1_right, r2_left, r2_right] =
            tc.draw(gs::arrays(gs::integers::<u32>().max_value(999)));
        let input = format!("{r1_left},{r1_right},{r2_left},{r2_right}");
        let result = ReadMaskSetting::from_str(&input).unwrap();

        assert_eq!(result.r1.from_start, r1_left);
        assert_eq!(result.r1.from_end, r1_right);
        assert_eq!(result.r2.from_start, r2_left);
        assert_eq!(result.r2.from_end, r2_right);
    }

    #[hegel::test]
    fn test_parse_read_mask_setting_invalid_count(tc: TestCase) {
        // Any count up to nine but four: skip over it instead of rejecting it.
        let count = tc.draw(gs::integers::<usize>().max_value(8));
        let count = if count >= 4 { count + 1 } else { count };
        let values: Vec<u32> =
            tc.draw(gs::vecs(gs::integers::<u32>().max_value(999)).min_size(count).max_size(count));

        let input = values.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",");

        let result = ReadMaskSetting::from_str(&input);
        assert!(result.is_err());
    }

    #[hegel::test]
    fn test_parse_read_mask_setting_invalid_integer(tc: TestCase) {
        let invalid_str = tc.draw(gs::from_regex("[a-zA-Z!@#$%^&*()]+"));
        let valid = tc.draw(gs::integers::<usize>().max_value(2));
        let valid_nums: Vec<u32> =
            tc.draw(gs::vecs(gs::integers::<u32>().max_value(999)).min_size(valid).max_size(valid));
        let mut parts = valid_nums.iter().map(|v| v.to_string()).collect::<Vec<_>>();
        parts.push(invalid_str);

        // Pad with more valid numbers if needed to get exactly 4 elements
        while parts.len() < 4 {
            parts.push("0".to_string());
        }
        parts.truncate(4);

        let input = parts.join(",");
        let result = ReadMaskSetting::from_str(&input);
        assert!(matches!(result, Err(ParseMaskError::ParseIntError(_))));
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

    fn read(pos: u32, len: u32, strand: Strand, direction: ReadDirection) -> SimpleRead {
        SimpleRead {
            strand,
            direction,
            position: crate::call::tests::pileup::PositionInRead { pos, read_length: len },
            ..default()
        }
    }

    #[test]
    fn filter_strand_orientation_logic() {
        let mask = ReadMaskParams {
            n_ot: ReadMaskSetting::from_str("5,3,2,4").unwrap(),
            n_ob: ReadMaskSetting::from_str("7,2,3,6").unwrap(),
        };

        // OT forward: normal masking (start=5, end=3)
        assert!(mask.filter(&read(10, 100, Strand::OT, ReadDirection::Forward)));
        assert!(!mask.filter(&read(4, 100, Strand::OT, ReadDirection::Forward))); // < start
        assert!(!mask.filter(&read(97, 100, Strand::OT, ReadDirection::Forward))); // > len-end-1

        // OT reverse: flipped masking (start=from_end=4, end=from_start=2)
        assert!(mask.filter(&read(10, 100, Strand::OT, ReadDirection::Reverse)));
        assert!(!mask.filter(&read(3, 100, Strand::OT, ReadDirection::Reverse))); // < from_end
        assert!(!mask.filter(&read(98, 100, Strand::OT, ReadDirection::Reverse))); // > len-from_start-1

        // OB reverse: flipped masking (start=from_end=2, end=from_start=7)
        assert!(mask.filter(&read(10, 100, Strand::OB, ReadDirection::Reverse)));
        assert!(!mask.filter(&read(1, 100, Strand::OB, ReadDirection::Reverse))); // < from_end
        assert!(!mask.filter(&read(93, 100, Strand::OB, ReadDirection::Reverse))); // > len-from_start-1

        // OB forward: normal masking (start=3, end=6)
        assert!(mask.filter(&read(10, 100, Strand::OB, ReadDirection::Forward)));
        assert!(!mask.filter(&read(2, 100, Strand::OB, ReadDirection::Forward))); // < start
        assert!(!mask.filter(&read(94, 100, Strand::OB, ReadDirection::Forward))); // > len-end-1
    }

    #[test]
    fn filter_size_and_boundary_conditions() {
        let mask =
            ReadMaskParams { n_ot: ReadMaskSetting::from_str("1,1,1,1").unwrap(), ..default() };

        // Too small reads
        assert!(!mask.filter(&read(1, 2, Strand::OT, ReadDirection::Forward))); // 2 < 1+1+1

        // Boundary conditions
        assert!(!mask.filter(&read(0, 10, Strand::OT, ReadDirection::Forward))); // pos < start
        assert!(mask.filter(&read(1, 10, Strand::OT, ReadDirection::Forward))); // pos >= start
        assert!(mask.filter(&read(8, 10, Strand::OT, ReadDirection::Forward))); // pos <= len-end-1
        assert!(!mask.filter(&read(9, 10, Strand::OT, ReadDirection::Forward))); // pos > len-end-1
    }

    #[test]
    fn filter_special_cases() {
        let mask = ReadMaskParams {
            n_ot: ReadMaskSetting::from_str("40,40,40,40").unwrap(),
            n_ob: ReadMaskSetting::from_str("40,40,40,40").unwrap(),
        };

        // No masking passes everything
        let no_mask = ReadMaskParams::default();
        assert!(no_mask.filter(&read(0, 1, Strand::OT, ReadDirection::Forward)));
        assert!(no_mask.filter(&read(99, 100, Strand::OB, ReadDirection::Reverse)));

        // Aggressive masking filters most
        assert!(mask.filter(&read(50, 100, Strand::OT, ReadDirection::Forward))); // 40 <= 50 <= 59
        assert!(!mask.filter(&read(39, 100, Strand::OT, ReadDirection::Forward))); // < 40
        assert!(!mask.filter(&read(60, 100, Strand::OT, ReadDirection::Forward))); // > 59
        assert!(!mask.filter(&read(50, 80, Strand::OT, ReadDirection::Forward))); // 80 < 81 (too small)
    }
}
