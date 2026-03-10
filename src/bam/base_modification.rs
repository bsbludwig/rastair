use color_eyre::eyre::{Context, Result};
use rastair_types::SmallVec;
use rastair_types::{Base, Strand};
use rust_htslib::bam::{Record, record::Aux};
use rustc_hash::FxHashMap;
use std::fmt::Write;
use tracing::debug;

/// Methylation context for XM tag annotation.
///
/// Determines which letter pair is used in the Bismark-style XM string:
/// - CpG → `z`/`Z`
/// - CHG → `x`/`X`
/// - CHH → `h`/`H`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MethylationContext {
    CpG,
    CHG,
    CHH,
}

impl MethylationContext {
    fn xm_char(self, methylated: bool) -> char {
        match (self, methylated) {
            (Self::CpG, true) => 'Z',
            (Self::CpG, false) => 'z',
            (Self::CHG, true) => 'X',
            (Self::CHG, false) => 'x',
            (Self::CHH, true) => 'H',
            (Self::CHH, false) => 'h',
        }
    }
}

/// Determine the methylation context for a cytosine at a given reference position.
///
/// For OT strand (C on top strand at `pos`):
/// - ref[pos+1] == G → CpG
/// - ref[pos+1] != G but ref[pos+2] == G → CHG
/// - otherwise → CHH
///
/// For OB strand (G on top strand at `pos`, i.e. C on bottom strand):
/// - ref[pos-1] == C → CpG
/// - ref[pos-1] != C but ref[pos-2] == C → CHG
/// - otherwise → CHH
pub fn determine_context(
    pos: u32,
    strand: Strand,
    ref_base: impl Fn(u32) -> Option<Base>,
) -> MethylationContext {
    match strand {
        Strand::OT => {
            if ref_base(pos + 1) == Some(Base::G) {
                MethylationContext::CpG
            } else if ref_base(pos + 2) == Some(Base::G) {
                MethylationContext::CHG
            } else {
                MethylationContext::CHH
            }
        }
        Strand::OB => {
            if pos > 0 && ref_base(pos - 1) == Some(Base::C) {
                MethylationContext::CpG
            } else if pos > 1 && ref_base(pos - 2) == Some(Base::C) {
                MethylationContext::CHG
            } else {
                MethylationContext::CHH
            }
        }
        Strand::Unknown => MethylationContext::CHH,
    }
}

/// Conversion type for XR/XG tags (CT or GA)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConversionType {
    CT,
    GA,
}

impl ConversionType {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::CT => "CT",
            Self::GA => "GA",
        }
    }
}

/// A struct representing a list of methylated positions in a sequence
#[derive(Debug, Clone)]
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
    /// Positions must be indices into `seq` (in stored/SEQ orientation).
    pub fn new(strand: Strand, seq: &[u8], methylated_positions: &[u32]) -> Self {
        let base = match strand {
            Strand::OT => Base::C,
            Strand::OB => Base::G,
            Strand::Unknown => {
                debug!("Unknown strand, cannot create MethylatedPositions");
                return Self { base: Base::C, strand: Strand::Unknown, positions: SmallVec::new() };
            }
        };

        // Skip list of positions of the base in the sequence
        //
        // e.g., for ACGTACGTACGT there is a C at positions 1, 5, 9. If the
        // second and third are methylated, we write it as 1,0.
        let positions = calculate_mm_skips(seq, base, methylated_positions);

        Self { base, strand, positions }
    }

    /// Apply the modification information to a BAM record
    ///
    /// There are two steps to this:
    ///
    /// 1. Rewrite the sequence to un-modify the methylated bases (T back to C, A back to G)
    /// 2. Add MM and ML tags to the record
    ///
    /// Apply the MM/ML modification tags to a BAM record.
    ///
    /// When there are no methylated positions, neither MM nor ML is written.
    /// An absent MM means "no modification data" per the SAM spec, which is
    /// correct for reads where no CpG positions had methylation evidence.
    ///
    /// Writing `C+m;` (empty position list) or an empty `ML:B:C` array causes
    /// modbedtools to segfault, so we omit the tags entirely in that case.
    /// The downside is that modkit cannot count unmethylated reads — those reads
    /// simply contribute no data to modkit's per-position counts.
    pub fn apply_to_record(&self, record: &mut Record) -> Result<()> {
        if self.positions.is_empty() {
            return Ok(());
        }

        record
            .push_aux(b"MM", Aux::String(&self.to_mod_string()))
            .wrap_err("could not apply modification to record")?;

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
    pub fn to_mod_string(&self) -> String {
        let mut mod_string = String::new();
        let strand = match self.strand {
            Strand::OT => "+",
            Strand::OB => "-",
            Strand::Unknown => {
                debug!("Unknown strand, cannot create mod string");
                return mod_string;
            }
        };

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

        for (i, &pos) in self.positions.iter().enumerate() {
            if i > 0 {
                mod_string.push(',');
            }

            write!(&mut mod_string, "{pos}").expect("Write to String failed");
        }

        mod_string.push(';');
        mod_string
    }
}

/// Calculate the skip list for MM:Z tag given a sequence and methylated positions
///
/// `methylated_indices` contains the 0-based indices of methylated bases in the
/// sequence.
fn calculate_mm_skips(seq: &[u8], base: Base, methylated_indices: &[u32]) -> SmallVec<u32, 10> {
    // Step 1: Find all positions of the target base in the sequence
    let base_positions: SmallVec<usize, 50> =
        seq.iter().enumerate().filter(|(_, b)| base == **b).map(|(idx, _)| idx).collect();

    // Step 2: Map methylated indices to their positions in the base-only list
    let mut methylated_positions: SmallVec<u32, 50> = methylated_indices
        .iter()
        .filter_map(|&meth_idx| {
            // Find where this index appears in the base_positions list
            base_positions.iter().position(|&pos| pos == meth_idx as usize).map(|base_pos_idx| {
                u32::try_from(base_pos_idx).expect("base position index fits in u32")
            })
        })
        .collect();

    // Sort to ensure proper ordering
    methylated_positions.sort_unstable();

    // Step 3: Calculate skip deltas
    let mut skips = SmallVec::new();
    let mut last_position = 0;

    for &pos in &methylated_positions {
        // The skip is the difference from the last position
        // For the first methylation, it's just the position itself
        skips.push(pos - last_position);
        // Update last_position to be one past the current position
        last_position = pos + 1;
    }

    skips
}

/// Generate XM tag string for a rewritten sequence (standard mode)
#[cfg(test)]
fn generate_xm_string(
    seq: &[u8],
    strand: Strand,
    methylated_positions: &[u32],
    methylated_set: &mut rustc_hash::FxHashSet<usize>,
) -> String {
    methylated_set.clear();
    methylated_set.extend(methylated_positions.iter().map(|&pos| pos as usize));

    let target_base = match strand {
        Strand::OT => Base::C,
        Strand::OB => Base::G,
        Strand::Unknown => return ".".repeat(seq.len()),
    };

    seq.iter()
        .enumerate()
        .map(|(i, &b)| {
            let base = Base::from(b);
            if base == target_base {
                if methylated_set.contains(&i) { 'Z' } else { 'z' }
            } else {
                '.'
            }
        })
        .collect()
}

/// Per-position annotation for a CpG site in the XM string.
#[derive(Debug, Clone, Copy)]
pub struct XmAnnotation {
    pub methylated: bool,
    pub context: MethylationContext,
}

/// Generate XM tag string for legacy mode (original, unrewritten sequence).
///
/// In legacy mode the SEQ is not rewritten, so methylated positions still show
/// T (OT) or A (OB). The XM tag uses context-dependent letters:
/// - CpG: `z`/`Z` (unmethylated/methylated)
/// - CHG: `x`/`X`
/// - CHH: `h`/`H`
/// - `.` for non-CpG positions
pub fn generate_xm_string_legacy(
    seq_len: usize,
    annotations: &FxHashMap<usize, XmAnnotation>,
) -> String {
    (0..seq_len)
        .map(|i| match annotations.get(&i) {
            Some(ann) => ann.context.xm_char(ann.methylated),
            None => '.',
        })
        .collect()
}

/// XR/XG/XM tags for legacy methylation format
#[derive(Debug, Clone)]
pub struct XrTags {
    /// XR:Z tag - read conversion (CT or GA)
    pub xr: ConversionType,
    /// XG:Z tag - reference conversion (CT or GA)
    pub xg: ConversionType,
    /// XM:Z tag - methylation call string
    pub xm: String,
}

impl XrTags {
    /// Create XR/XG/XM tags for legacy mode (original, unrewritten seq)
    pub fn new_legacy(
        seq_len: usize,
        strand: Strand,
        is_first_in_pair: bool,
        annotations: &FxHashMap<usize, XmAnnotation>,
    ) -> Self {
        let (xr, xg) = xr_xg(strand, is_first_in_pair);
        let xm = generate_xm_string_legacy(seq_len, annotations);
        Self { xr, xg, xm }
    }

    /// Apply XR/XG/XM tags to a BAM record
    pub fn apply_to_record(&self, record: &mut Record) -> Result<()> {
        record.push_aux(b"XR", Aux::String(self.xr.as_str())).wrap_err("could not add XR tag")?;

        record.push_aux(b"XG", Aux::String(self.xg.as_str())).wrap_err("could not add XG tag")?;

        record.push_aux(b"XM", Aux::String(&self.xm)).wrap_err("could not add XM tag")?;

        Ok(())
    }
}

fn xr_xg(strand: Strand, is_first_in_pair: bool) -> (ConversionType, ConversionType) {
    let xr = if is_first_in_pair { ConversionType::CT } else { ConversionType::GA };
    let xg = match strand {
        Strand::OT => ConversionType::CT,
        Strand::OB => ConversionType::GA,
        Strand::Unknown => {
            debug!("Unknown strand detected, defaulting XG to CT");
            ConversionType::CT
        }
    };
    (xr, xg)
}

#[cfg(test)]
mod tests {
    use super::*;
    use insta::{assert_compact_debug_snapshot, assert_snapshot};
    use rustc_hash::FxHashSet;
    use std::iter::repeat_n;

    #[test]
    fn skip_list() {
        // boring seq
        let skip_list = calculate_mm_skips(b"CCCCC", Base::C, &[1, 2]);
        assert_compact_debug_snapshot!(skip_list, @"[1, 0]");

        // long seq
        let skip_list = calculate_mm_skips(
            &repeat_n(b'A', 10).chain(repeat_n(b'C', 10)).collect::<Vec<u8>>(),
            Base::C,
            &[10, 15],
        );
        assert_compact_debug_snapshot!(skip_list, @"[0, 4]");
    }

    #[test]
    fn skip_list_edge_cases() {
        // not in seq
        let skip_list = calculate_mm_skips(b"AAAAAA", Base::C, &[1, 2]);
        assert_compact_debug_snapshot!(skip_list, @"[]");

        // empty seq
        let skip_list = calculate_mm_skips(b"", Base::C, &[1, 2]);
        assert_compact_debug_snapshot!(skip_list, @"[]");

        // empty list
        let skip_list = calculate_mm_skips(b"CCCCC", Base::C, &[]);
        assert_compact_debug_snapshot!(skip_list, @"[]");
    }

    #[test]
    fn mod_string() {
        let meth_pos = MethylatedPositions::new(Strand::OT, b"CCCCCCC", &[3, 4]);
        let mod_str = meth_pos.to_mod_string();
        assert_snapshot!(mod_str, @"C+m,3,0;");
    }

    #[test]
    fn test_xm_string_generation() {
        // Test XM tag generation for OT strand
        // Sequence: ACGTCGATCG with Cs at positions 1, 4, 8
        // Methylated positions: 1, 8 (C at index 1 and 8)
        let seq = b"ACGTCGATCG";
        let methylated_positions = vec![1, 8];
        let mut methylated_set = FxHashSet::default();
        let xm = generate_xm_string(seq, Strand::OT, &methylated_positions, &mut methylated_set);
        // Expected: .Z..z...Z. (Z = methylated CpG C, z = unmethylated CpG C, . = non-C)
        assert_snapshot!(xm, @".Z..z...Z.");
    }

    #[test]
    fn test_xm_string_ob_strand() {
        // Test XM tag generation for OB strand
        // Sequence: GCGATATGCG with Gs at positions 0, 2, 7, 9
        // Methylated positions: 0, 7 (G at index 0 and 7)
        let seq = b"GCGATATGCG";
        let methylated_positions = vec![0, 7];
        let mut methylated_set = FxHashSet::default();
        let xm = generate_xm_string(seq, Strand::OB, &methylated_positions, &mut methylated_set);
        // Expected: Z.z...Z.z. (Z = methylated CpG G, z = unmethylated CpG G, . = non-G)
        assert_snapshot!(xm, @"Z.z....Z.z");
    }

    #[test]
    fn test_xm_string_no_methylation() {
        // Test XM tag generation with no methylation
        let seq = b"ACGTCGATCG";
        let methylated_positions = vec![];
        let mut methylated_set = FxHashSet::default();
        let xm = generate_xm_string(seq, Strand::OT, &methylated_positions, &mut methylated_set);
        // All Cs should be lowercase z (unmethylated)
        assert_snapshot!(xm, @".z..z...z.");
    }

    #[test]
    fn test_xm_string_all_methylated() {
        // Test XM tag generation with all positions methylated
        let seq = b"CGCGCG";
        let methylated_positions = vec![0, 2, 4];
        let mut methylated_set = FxHashSet::default();
        let xm = generate_xm_string(seq, Strand::OT, &methylated_positions, &mut methylated_set);
        // All Cs should be uppercase Z (methylated)
        assert_snapshot!(xm, @"Z.Z.Z.");
    }

    #[test]
    fn context_ot_cpg() {
        // Reference: C at pos 10, G at pos 11 → CpG
        use Base::*;
        let ref_bases = FxHashMap::from_iter([(10, C), (11, G), (12, A)]);
        let lookup = |p: u32| ref_bases.get(&p).copied();
        assert_eq!(determine_context(10, Strand::OT, lookup), MethylationContext::CpG);
    }

    #[test]
    fn context_ot_chg() {
        // Reference: C at pos 10, A at pos 11, G at pos 12 → CHG
        use Base::*;
        let ref_bases = FxHashMap::from_iter([(10, C), (11, A), (12, G)]);
        let lookup = |p: u32| ref_bases.get(&p).copied();
        assert_eq!(determine_context(10, Strand::OT, lookup), MethylationContext::CHG);
    }

    #[test]
    fn context_ot_chh() {
        // Reference: C at pos 10, A at pos 11, T at pos 12 → CHH
        use Base::*;
        let ref_bases = FxHashMap::from_iter([(10, C), (11, A), (12, T)]);
        let lookup = |p: u32| ref_bases.get(&p).copied();
        assert_eq!(determine_context(10, Strand::OT, lookup), MethylationContext::CHH);
    }

    #[test]
    fn context_ob_cpg() {
        // Reference: G at pos 10, C at pos 9 → CpG on bottom strand
        use Base::*;
        let ref_bases = FxHashMap::from_iter([(9, C), (10, G)]);
        let lookup = |p: u32| ref_bases.get(&p).copied();
        assert_eq!(determine_context(10, Strand::OB, lookup), MethylationContext::CpG);
    }

    #[test]
    fn context_ob_chg() {
        // Reference: G at pos 10, T at pos 9, C at pos 8 → CHG on bottom strand
        use Base::*;
        let ref_bases = FxHashMap::from_iter([(8, C), (9, T), (10, G)]);
        let lookup = |p: u32| ref_bases.get(&p).copied();
        assert_eq!(determine_context(10, Strand::OB, lookup), MethylationContext::CHG);
    }

    #[test]
    fn context_ob_chh() {
        // Reference: G at pos 10, T at pos 9, A at pos 8 → CHH on bottom strand
        use Base::*;
        let ref_bases = FxHashMap::from_iter([(8, A), (9, T), (10, G)]);
        let lookup = |p: u32| ref_bases.get(&p).copied();
        assert_eq!(determine_context(10, Strand::OB, lookup), MethylationContext::CHH);
    }

    #[test]
    fn xm_legacy_context_letters() {
        let mut annotations = FxHashMap::default();
        // pos 0: methylated CpG → Z
        annotations.insert(0, XmAnnotation { methylated: true, context: MethylationContext::CpG });
        // pos 2: unmethylated CpG → z
        annotations.insert(2, XmAnnotation { methylated: false, context: MethylationContext::CpG });
        // pos 4: methylated CHG → X
        annotations.insert(4, XmAnnotation { methylated: true, context: MethylationContext::CHG });
        // pos 6: unmethylated CHH → h
        annotations.insert(6, XmAnnotation { methylated: false, context: MethylationContext::CHH });
        // pos 8: methylated CHH → H
        annotations.insert(8, XmAnnotation { methylated: true, context: MethylationContext::CHH });

        let xm = generate_xm_string_legacy(10, &annotations);
        assert_snapshot!(xm, @"Z.z.X.h.H.");
    }
}
