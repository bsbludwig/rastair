use crate::{call::variant_calling::indel_calling::IndelCall, utils::Base, utils::Strand};
use seqair_types::SmallVec;

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
    /// Terminal tandem repeat or soft-clip: this fragment's alignment is the kind
    /// that slips. Subtracted from the alternate side of the ratio.
    pub noisy: bool,
}

/// Aggregated indel counts at a position, ready for calling.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct IndelCounts {
    /// Each unique indel allele with its forward and reverse strand counts.
    pub alleles: SmallVec<IndelAlleleCounts, 2>,
    /// Total reads at this position that do NOT have an indel (reference-supporting).
    pub ref_count: u32,
    /// Reference-supporting fragments that look noisy. The reference-side
    /// counterpart of [`IndelAlleleCounts::noisy`]; both come off together so the
    /// noise exclusion cannot move the ratio on its own.
    ///
    /// A fragment carrying the indel is excluded even when one of its mates is
    /// soft-clipped over it: `ref_count` has already dropped that fragment via
    /// `total_indel_reads`, so counting it here too would subtract it twice.
    pub noisy_ref_count: u32,
}

impl IndelCounts {
    pub fn total_indel_reads(&self) -> u32 {
        self.alleles.iter().map(|a| a.total()).sum()
    }

    /// Depth with repeat-noisy reads removed from *both* sides of the ratio.
    ///
    /// `noisy_ref_count` covers only the reference side, so subtracting it alone
    /// would shrink the denominator and inflate the alternate fraction by the noise
    /// rate — worst exactly in the repeats where indels concentrate. Hence the
    /// alternate side comes in via `clean_total`, which drops its own noisy reads.
    pub fn clean_depth(&self) -> u32 {
        self.ref_count.saturating_sub(self.noisy_ref_count)
            + self.alleles.iter().map(|a| a.clean_total()).sum::<u32>()
    }

    pub fn is_empty(&self) -> bool {
        self.alleles.is_empty()
    }
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IndelAlleleCounts {
    pub allele: IndelAllele,
    /// Supporting reads by *bisulfite* strand, not by the alignment reverse flag.
    /// Both mates of a fragment share an OT/OB assignment but have opposite reverse
    /// flags, so a reverse-flag split is present on almost every real allele and
    /// carries no signal.
    pub ot: u32,
    pub ob: u32,
    pub unknown_strand: u32,
    /// Supporting reads carrying a terminal tandem repeat. The alternate-side
    /// counterpart of [`IndelCounts::noisy_ref_count`].
    pub noisy: u32,
}

impl IndelAlleleCounts {
    pub fn total(&self) -> u32 {
        self.ot + self.ob + self.unknown_strand
    }

    pub fn clean_total(&self) -> u32 {
        self.total().saturating_sub(self.noisy)
    }

    pub fn on_both_strands(&self) -> bool {
        self.ot > 0 && self.ob > 0
    }
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct IndelData {
    pub observations: SmallVec<IndelObservation, 0>,
    pub ref_window: SmallVec<Base, { super::INDEL_REF_WINDOW_LEN }>,
    pub ref_anchor: u8,
    pub homopolymer_run: u8,
    pub dinucleotide_run: u8,
    pub soft_clip_count: u32,
    pub counts: IndelCounts,
    pub calls: Vec<IndelCall>,
}
