use super::Base;
use color_eyre::eyre::{Context, Result};
use rust_htslib::bam::{Record, record::Aux};
use smallvec::SmallVec;
use std::fmt::Write;

/// A struct representing a list of methylated positions in a sequence
pub struct MethylatedPositions {
    /// Unmodified "fundamental" base (on top strand)
    pub base: Base,
    /// List of positions (0-based) of the base in the sequence
    pub positions: SmallVec<u32, 10>,
}

impl MethylatedPositions {
    pub fn apply_to_record(&self, record: &mut Record) -> Result<()> {
        let strand = record.strand();
        record
            .push_aux(b"Mm", Aux::String(&self.to_mod_string(strand.strand_symbol())))
            .wrap_err("could not apply modification to record")
    }

    /// Write the modification string for this base and positions
    ///
    /// See [SAM tags], ch. 1.7 for details
    ///
    /// [SAM tags]: https://samtools.github.io/hts-specs/SAMtags.pdf
    pub fn to_mod_string(&self, strand: &str) -> String {
        let mut mod_string = String::new();
        // unmodified "fundamental" base on top strand
        mod_string.push(self.base.as_char());
        // strand: + or -
        mod_string.push_str(strand);
        // 5-Methylcytosine
        mod_string.push_str("m,");

        let mut prev_pos = None;
        for pos in &self.positions {
            let steps_between_prev_and_this = match prev_pos {
                Some(prev) => pos.saturating_sub(prev).saturating_sub(1),
                None => *pos,
            };
            let _ = prev_pos.insert(*pos);
            write!(&mut mod_string, "{steps_between_prev_and_this},")
                .expect("Write to String failed");
        }

        // replace last comma
        mod_string.pop();
        mod_string.push(';');

        mod_string
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn test_empty_positions() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::new() };
        let expected = "C+m;";
        assert_eq!(expected, input.to_mod_string("+"));
    }

    #[test]
    fn test_single_position() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::from([4]) };
        let expected = "C+m,4;";
        assert_eq!(expected, input.to_mod_string("+"));
    }

    #[test]
    fn test_multiple_positions() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::from([5, 18, 19]) };
        let expected = "C+m,5,12,0;";
        assert_eq!(expected, input.to_mod_string("+"));
    }

    #[test]
    fn test_consecutive_positions() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::from([0, 1, 2]) };
        let expected = "C+m,0,0,0;";
        assert_eq!(expected, input.to_mod_string("+"));
    }

    #[test]
    fn test_apply_to_record() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::from([0, 1, 2]) };

        let mut record = Record::new();
        input.apply_to_record(&mut record).unwrap();
        let output = record.aux(b"Mm").unwrap();
        let Aux::String(output_str) = output else { panic!("Expected string aux type") };

        let expected = "C+m,0,0,0;";
        assert_eq!(expected, output_str);
    }

    proptest! {
        #[test]
        fn test_mod_string_roundtrip(
            base in prop_oneof![Just(Base::A), Just(Base::C), Just(Base::G), Just(Base::T)],
            mut positions in prop::collection::vec(0..100u32, 0..10),
            strand in prop_oneof![Just("+"), Just("-")]
        ) {
            positions.sort();
            let input = MethylatedPositions { base, positions: SmallVec::from(positions) };
            let mod_string = input.to_mod_string(strand);
            assert!(mod_string.starts_with(&format!("{}{}", base.as_char(), strand)));
            assert!(mod_string.ends_with(";"));
            // assert right number of positions
            let positions_in_str = mod_string.trim_end_matches(';')
                .split(',')
                .skip(1)
                .count();
            assert_eq!(input.positions.len(), positions_in_str);
        }
    }
}
