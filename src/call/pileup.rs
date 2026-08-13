use crate::{sequence::ChunkRegion, utils::Base, vcf::SequenceContext};
use seqair_types::{SmallVec, SmolStr};

pub mod indels;
mod read;
pub use read::*;
pub(crate) mod from_hts;
pub(crate) mod overlapping_reads;

/// Reference bases kept upstream / downstream of the anchor for indel slippage
/// detection. Downstream must span the indel plus a few repeat units; a little
/// upstream covers reads that are not left-aligned.
pub(crate) const INDEL_REF_WINDOW_UP: usize = 8;
pub(crate) const INDEL_REF_WINDOW_DOWN: usize = 24;
/// Maximum window length (upstream + anchor + downstream). The inline capacity
/// of [`Pileup::indel_ref_window`] is sized to this so a populated window never
/// spills to the heap.
pub(crate) const INDEL_REF_WINDOW_LEN: usize = INDEL_REF_WINDOW_UP + 1 + INDEL_REF_WINDOW_DOWN;

/// Rastair's representation of a pileup at a specific position in the genome
///
/// Deliberately carries no `#[serde(default)]`: the `.mpk` payload encodes structs
/// as positional arrays, so a field added anywhere but the end shifts every field
/// after it and a default would decode neighbouring data instead of failing.
/// [`crate::io::mpk::format::MPK_FORMAT_VERSION`] is the compatibility mechanism —
/// bump it when this struct changes.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Pileup {
    /// Region of the chunk this pileup belongs to
    pub region: ChunkRegion,
    /// Sequence context around the position in the reference
    pub context: SequenceContext,
    /// Position in the sequence, 0-based
    pub pos: u32,
    // if this becomes slow, consider boxing or Rc-ing this
    pub reads: SimpleReads,
    /// Reference base at this position
    pub reference_base: Base,
    /// Indel observations collected from reads at this position.
    /// Empty at most positions — `SmallVec<_, 0>` avoids heap allocation when empty.
    pub indel_observations: SmallVec<indels::IndelObservation, 0>,
    /// Non-indel-supporting reads with a terminal tandem repeat, subtracted from
    /// depth by the ML indel pathway.
    pub depth_offset: u32,
    pub homopolymer_run: u8,
    pub dinucleotide_run: u8,
    pub soft_clip_count: u32,
    /// Non-indel-supporting reads the hard-filter indel pathway excludes as noise
    /// (terminal tandem repeat or soft-clip). Paired with the per-allele
    /// `IndelObservation::noisy` so the exclusion applies to both sides of the
    /// alt/depth ratio; see [`indels::IndelCounts::noisy_ref_count`].
    pub noisy_ref_count: u32,
    /// Reference bases around the anchor (`indel_ref_anchor` is the anchor's
    /// index), used for tandem-repeat / slippage detection of indel alleles.
    /// Only populated when indel observations are present.
    pub indel_ref_window: SmallVec<Base, INDEL_REF_WINDOW_LEN>,
    /// Index of the anchor (pileup `pos`) base within [`Pileup::indel_ref_window`].
    pub indel_ref_anchor: u8,
}

impl Pileup {
    /// Chromosome name of the segment
    pub fn contig(&self) -> SmolStr {
        self.region.contig.clone()
    }

    /// Position in the segment sequence, 0-based
    pub fn idx(&self) -> usize {
        self.region.pos_to_idx(self.pos).expect("valid position")
    }

    /// Reference base right before the variant position
    pub fn ref_before(&self) -> Option<Base> {
        self.context.before_1
    }

    /// Reference base right after the variant position
    pub fn ref_after(&self) -> Option<Base> {
        self.context.after_1
    }

    pub fn alleles(&self) -> SmallVec<Base, 4> {
        let mut res = SmallVec::new();
        res.push(self.reference_base);
        self.reads.iter().map(|b| b.base).fold(res, |mut acc, base| {
            if !acc.contains(&base) {
                acc.push(base);
            }
            acc
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        sequence::{ChunkRegion, Region, Segment},
        utils::default,
    };
    use insta::assert_debug_snapshot;
    use seqair_types::Strand;

    #[test]
    fn test_alleles_in_order() {
        let bases = SimpleReads(
            vec![
                SimpleRead { base: Base::A, strand: Strand::OT, ..default() },
                SimpleRead {
                    base: Base::C,
                    strand: Strand::OB,
                    reverse: true,
                    second: false,
                    ..default()
                },
                SimpleRead { base: Base::A, strand: Strand::OT, ..default() },
            ]
            .into(),
        );

        let segment = Segment {
            range: ChunkRegion {
                region: Region { contig: "chr19".into(), start: 1000, end: 1100 },
                last_position: 2000,
                overlap_start: 0,
                overlap_end: 0,
            },
            sequence: vec![],
            overlap_start: 0,
            overlap_end: 0,
        };
        let variant_candidate = Pileup {
            region: segment.range.clone(),
            context: SequenceContext::default(),
            pos: 1002, // Corresponds to index in the segment
            reads: bases,
            reference_base: Base::T, // Assume T is the reference base at this position
            indel_observations: default(),
            depth_offset: 0,
            homopolymer_run: 0,
            dinucleotide_run: 0,
            soft_clip_count: 0,
            noisy_ref_count: 0,
            indel_ref_window: default(),
            indel_ref_anchor: 0,
        };

        let alleles = variant_candidate.alleles();
        assert_eq!(alleles[0], Base::T, "Reference base should be first");

        assert_debug_snapshot!(alleles, @r"
        [
            T,
            A,
            C,
        ]
        ");
    }
}
