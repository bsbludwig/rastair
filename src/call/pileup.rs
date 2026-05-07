use crate::{sequence::ChunkRegion, utils::Base, vcf::SequenceContext};
use rastair_types::{SmallVec, SmolStr};

pub mod indels;
mod read;
pub use read::*;
pub(crate) mod from_hts;
pub(crate) mod overlapping_reads;

/// Rastair's representation of a pileup at a specific position in the genome
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
    #[serde(default)]
    pub indel_observations: SmallVec<indels::IndelObservation, 0>,
    /// Number of reference reads with problematic patterns (homopolymer, soft-clip)
    /// for indel depth adjustment.
    #[serde(default)]
    pub depth_offset: u32,
    #[serde(default)]
    pub homopolymer_run: u8,
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
    use rastair_types::Strand;

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
