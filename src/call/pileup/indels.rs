use crate::utils::{Base, Strand};
use rastair_types::SmallVec;

/// A specific indel allele observed in reads.
#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum IndelAllele {
    /// Bases inserted after the anchor position (not including the anchor).
    Insertion(SmallVec<Base, 4>),
    /// Reference bases deleted after the anchor position.
    Deletion(SmallVec<Base, 4>),
}

impl IndelAllele {
    pub fn len(&self) -> usize {
        match self {
            Self::Insertion(bases) | Self::Deletion(bases) => bases.len(),
        }
    }

    pub fn bases(&self) -> &[Base] {
        match self {
            Self::Insertion(bases) | Self::Deletion(bases) => bases,
        }
    }

    pub fn is_insertion(&self) -> bool {
        matches!(self, Self::Insertion(_))
    }

    pub fn is_deletion(&self) -> bool {
        matches!(self, Self::Deletion(_))
    }
}

/// A single read's indel observation at a pileup position.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IndelObservation {
    pub allele: IndelAllele,
    pub strand: Strand,
    pub reverse: bool,
    pub pos_in_read: u32,
    pub read_length: u32,
    pub mapq: u8,
    pub base_qual: u8,
    pub matching_bases: u32,
    pub num_indels_in_read: u32,
    pub insertion_base_quals: SmallVec<u8, 4>,
    pub post_del_base_qual: u8,
    pub has_repeat: bool,
}

/// Aggregated indel counts at a position, ready for calling.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct IndelCounts {
    /// Each unique indel allele with its forward and reverse strand counts.
    pub alleles: SmallVec<IndelAlleleCounts, 2>,
    /// Total reads at this position that do NOT have an indel (reference-supporting).
    pub ref_count: u32,
    /// Reads with problematic patterns (homopolymer, soft-clip) subtracted from depth
    /// for indel quality calculation.
    pub depth_offset: u32,
}

impl IndelCounts {
    pub fn total_indel_reads(&self) -> u32 {
        self.alleles.iter().map(|a| a.total()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.alleles.is_empty()
    }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IndelAlleleCounts {
    pub allele: IndelAllele,
    pub fwd: u32,
    pub rev: u32,
}

impl IndelAlleleCounts {
    pub fn total(&self) -> u32 {
        self.fwd + self.rev
    }

    pub fn on_both_strands(&self) -> bool {
        self.fwd > 0 && self.rev > 0
    }
}
