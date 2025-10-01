use color_eyre::eyre::{Context, Result};
use rastair_types::{Base, Strand, StrandFromRecord};
use rust_htslib::bam::{Record, record::Aux};
use smallvec::SmallVec;
use std::fmt::Write;
use tracing::debug;

/// A struct representing a list of methylated positions in a sequence
pub struct MethylatedPositions {
    /// Unmodified "fundamental" base (on top strand)
    pub base: Base,
    pub strand: Strand,
    /// List of positions (0-based) of the base in the sequence
    ///
    /// Comma separated list of how many seq bases of the stated base type to
    /// skip, stored as a delta to the last and starting with 0 as the first (or
    /// next) base, starting from the uncomplemented 5’ end of the SEQ field.
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
    pub fn new(strand: Strand, seq: &[u8], methylated_positions: &[u32]) -> Self {
        let (base, strand) = match strand {
            Strand::OT => (Base::C, Strand::OT),
            Strand::OB => (Base::G, Strand::OB),
            Strand::Unknown => {
                debug!("Unknown strand, cannot create MethylatedPositions");
                return Self { base: Base::C, strand: Strand::Unknown, positions: SmallVec::new() };
            }
        };

        // Skip list of positions of the base in the sequence
        //
        // e.g., for ACGTACGTACGT there is a C at positions 1, 5, 9. If the
        // second and third are methylated, we write it as 1,0.
        let positions = {
            let mut base_count = 0;
            let mut skip_list = SmallVec::new();
            for (i, &b) in seq.iter().enumerate() {
                if base == b {
                    if methylated_positions.contains(&(i as u32)) {
                        skip_list.push(base_count);
                        base_count = 0;
                    } else {
                        base_count += 1;
                    }
                }
            }
            skip_list
        };

        Self { base, strand, positions }
    }

    /// Apply the modification information to a BAM record
    ///
    /// There are two steps to this:
    ///
    /// 1. Rewrite the sequence to un-modify the methylated bases (T back to C, A back to G)
    /// 2. Add MM and ML tags to the record
    pub fn apply_to_record(&self, record: &mut Record) -> Result<()> {
        let strand_symbol = match StrandFromRecord::strand(record) {
            Strand::OT => "+",
            Strand::OB => "-",
            Strand::Unknown => {
                debug!(flags = record.flags(), "Unknown strand, not modifying record");
                return Ok(());
            }
        };

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

        // mod_string.push('.');
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
