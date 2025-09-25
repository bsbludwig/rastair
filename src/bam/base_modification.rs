use color_eyre::eyre::{Context, Result};
use rastair_types::{Base, Strand};
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
    /// Create methylation positions for CpG sites from a read record
    ///
    /// This determines the correct base and strand for CpG methylation based on
    /// the actual sequence content and read orientation
    // FIMXE: This is wrong
    // - reverse should have seq reversed?
    // - positions are:  comma separated list of how many seq bases of the stated base type to skip, stored as a delta to the last and starting with 0 as the first (or next) base, starting from the uncomplemented 5’ end of the SEQ field.
    // work with modkit code -- MmTagInfo, DeltaListConverter, test_get_base_mod_probs
    pub fn for_cpg_methylation(record: &Record, positions: SmallVec<u32, 10>) -> Self {
        let is_reverse = record.is_reverse();

        // For CpG methylation, we need to determine if we're looking at C or G
        // On forward strand: C gets methylated (5mC)
        // On reverse strand: G represents the C on the opposite strand that got methylated
        let base = if is_reverse { Base::G } else { Base::C };

        let seq = record.seq();

        // Validate that positions actually contain the expected bases
        let valid_positions: SmallVec<u32, 10> = positions
            .into_iter()
            .filter(|&pos| {
                if pos as usize >= seq.len() {
                    return false;
                }
                let observed_base = seq[pos as usize];
                match base {
                    Base::C => observed_base == b'C',
                    Base::G => observed_base == b'G',
                    _ => false,
                }
            })
            .collect();

        if valid_positions.is_empty() {
            Self { base, positions: SmallVec::<u32, 10>::new() }
        } else {
            Self { base, positions: valid_positions }
        }
    }

    /// Apply the modification information to a BAM record
    ///
    /// There are two steps to this:
    ///
    /// 1. Rewrite the sequence to un-modify the methylated bases (T back to C, A back to G)
    /// 2. Add MM and ML tags to the record
    pub fn apply_to_record(&self, record: &mut Record) -> Result<()> {
        // WIP - not yet rewriting the sequence
        let is_reverse = record.is_reverse();
        let strand_symbol = if is_reverse { "-" } else { "+" };

        record
            .push_aux(b"MM", Aux::String(&self.to_mod_string(strand_symbol)))
            .wrap_err("could not apply modification to record")?;

        // also apply ML dummy tag with all 255 (unknown probability)
        record
            .push_aux(b"ML", Aux::ArrayU8(vec![255_u8; self.positions.len()].as_slice().into()))
            .wrap_err("could not apply ML tag to record")?;

        Ok(())
    }

    /// Write the modification string for this base and positions
    ///
    /// See [SAM tags], ch. 1.7 for details
    /// Format: `BASE STRAND MOD_CODE , DELTA , DELTA , ... ;`
    /// Where deltas are distances between consecutive positions of the base type
    ///
    /// [SAM tags]: https://samtools.github.io/hts-specs/SAMtags.pdf
    pub fn to_mod_string(&self, strand: &str) -> String {
        let mut mod_string = String::new();

        // Unmodified "fundamental" base on sequenced strand
        mod_string.push(self.base.as_char());

        // Strand: + or -
        mod_string.push_str(strand);

        // Modification code: 'm' for 5-Methylcytosine
        mod_string.push('m');

        if self.positions.is_empty() {
            // Empty coordinate list - indicates this modification type is not present
            mod_string.push(';');
            return mod_string;
        }

        mod_string.push('.');
        mod_string.push(',');

        // Encode positions as deltas between consecutive occurrences of this base type
        let mut prev_pos = None;
        for (i, &pos) in self.positions.iter().enumerate() {
            if i > 0 {
                mod_string.push(',');
            }

            let delta = match prev_pos {
                Some(prev) => pos.saturating_sub(prev).saturating_sub(1),
                None => pos,
            };

            write!(&mut mod_string, "{delta}").expect("Write to String failed");

            prev_pos = Some(pos);
        }

        mod_string.push(';');
        mod_string
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use insta::assert_snapshot;
    use proptest::prelude::*;

    #[test]
    fn test_empty_positions() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::new() };
        assert_snapshot!(input.to_mod_string("-"), @"C-m;");
    }

    #[test]
    fn test_single_position() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::from([4]) };
        assert_snapshot!(input.to_mod_string("+"), @"C+m.,4;");
    }

    #[test]
    fn test_multiple_positions() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::from([5, 18, 19]) };
        assert_snapshot!(input.to_mod_string("+"), @"C+m.,5,12,0;");
    }

    #[test]
    fn test_consecutive_positions() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::from([0, 1, 2]) };
        assert_snapshot!(input.to_mod_string("+"), @"C+m.,0,0,0;");
    }

    #[test]
    fn test_reverse_strand() {
        let input = MethylatedPositions { base: Base::G, positions: SmallVec::from([10, 15]) };
        assert_snapshot!(input.to_mod_string("-"), @"G-m.,10,4;");
    }

    #[test]
    fn test_apply_to_record() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::from([0, 1, 2]) };

        let mut record = Record::new();
        input.apply_to_record(&mut record).unwrap();
        let output = record.aux(b"MM").unwrap();
        let Aux::String(output_str) = output else { panic!("Expected string aux type") };

        assert_snapshot!(output_str, @"C+m.,0,0,0;");
    }

    proptest! {
        #[test]
        fn test_mod_string_roundtrip(
            base in prop_oneof![Just(Base::C), Just(Base::G)],
            mut positions in prop::collection::vec(0..100u32, 0..10),
            strand in prop_oneof![Just("+"), Just("-")]
        ) {
            positions.sort();
            positions.dedup();
            let input = MethylatedPositions { base, positions: SmallVec::from(positions) };
            let mod_string = input.to_mod_string(strand);

            // Basic format validation
            assert!(mod_string.starts_with(&format!("{}{}", base.as_char(), strand)));
            assert!(mod_string.contains("m"));
            assert!(mod_string.ends_with(";"));

            // Count positions in string
            if input.positions.is_empty() {
                assert_eq!("C+m;" .len(), mod_string.len().min(4)); // Allow for G-m; etc
            } else {
                let positions_in_str = mod_string
                    .trim_end_matches(';')
                    .split(',')
                    .skip(1) // Skip the "Xm" part
                    .count();
                assert_eq!(input.positions.len(), positions_in_str);
            }
        }
    }
}
