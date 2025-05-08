use super::Base;
use color_eyre::eyre::{Context, Result};
use rust_htslib::bam::{Record, record::Aux};
use smallvec::SmallVec;
use std::fmt::Write;

pub struct MethylatedPositions {
    pub base: Base,
    pub positions: SmallVec<u32, 10>,
}

impl MethylatedPositions {
    pub fn apply_to_record(&self, record: &mut Record) -> Result<()> {
        let strand = record.strand();
        record
            .push_aux(b"Mm", Aux::String(&self.to_mod_string(strand.strand_symbol())))
            .wrap_err("could not apply modification to record")
    }

    pub fn to_mod_string(&self, strand: &str) -> String {
        let mut mod_string = String::new();
        // unmodified "fundamental" base on top strand
        mod_string.push(self.base.as_char());
        // strand: + or -
        mod_string.push_str(strand);
        // 5-Methylcytosine
        mod_string.push_str("m,");

        // list modified bases as skip list
        let mut last_pos = 0;
        for pos in &self.positions {
            let index = pos - last_pos - 1;
            last_pos = *pos;
            write!(&mut mod_string, "{index},").unwrap();
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

    #[test]
    fn test_positive_strand_multiple_positions() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::from([6, 19, 20]) };
        let expected = "C+m,5,12,0;";
        assert_eq!(expected, input.to_mod_string("+"));
    }

    #[test]
    fn test_negative_strand_multiple_positions() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::from([6, 19, 20]) };
        let expected = "C-m,5,12,0;";
        assert_eq!(expected, input.to_mod_string("-"));
    }

    #[test]
    fn test_single_position() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::from([5]) };
        let expected = "C+m,4;";
        assert_eq!(expected, input.to_mod_string("+"));
    }

    #[test]
    fn test_consecutive_positions() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::from([1, 2, 3]) };
        let expected = "C+m,0,0,0;";
        assert_eq!(expected, input.to_mod_string("+"));
    }

    #[test]
    fn test_empty_positions() {
        let input = MethylatedPositions { base: Base::C, positions: SmallVec::new() };
        let expected = "C+m;";
        assert_eq!(expected, input.to_mod_string("+"));
    }
}
