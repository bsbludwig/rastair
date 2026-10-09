use crate::{
    call::{
        pair_hmm::HaplotypeScores,
        variant_calling::indel_calling::{IndelRecords, call_on_pair_hmm},
    },
    utils::Base,
};
use seqair_types::SmallVec;
use std::sync::Arc;

/// A specific indel allele observed in reads.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
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
}

/// The indel alleles a column's reads carry, one observation per fragment.
/// They decide which alleles the pair-HMM scores; the record's evidence is
/// the pair-HMM's.
#[derive(Debug, Clone, Default)]
pub struct IndelCounts {
    pub alleles: SmallVec<IndelAlleleCounts, 2>,
}

impl IndelCounts {
    pub fn total_indel_reads(&self) -> u32 {
        self.alleles.iter().map(|a| a.total).sum()
    }

    #[cfg(feature = "experimental-seqair")]
    pub(crate) fn aggregate(observations: &[IndelAllele]) -> Self {
        let mut alleles: SmallVec<IndelAlleleCounts, 2> = SmallVec::new();
        for observed in observations {
            if let Some(entry) = alleles.iter_mut().find(|e| e.allele == *observed) {
                entry.total += 1;
            } else {
                alleles.push(IndelAlleleCounts { allele: observed.clone(), total: 1 });
            }
        }
        Self { alleles }
    }
}

#[derive(Debug, Clone)]
pub struct IndelAlleleCounts {
    pub allele: IndelAllele,
    pub total: u32,
}

#[derive(Debug, Clone, Default)]
pub struct IndelData {
    pub counts: IndelCounts,
    /// `None` when no allele had the support to be scored, or no read to score it on.
    pub scored: Option<Scored>,
}

/// The pair-HMM's scores of a column's leading alleles, and the calls made on them.
#[derive(Debug, Clone)]
pub struct Scored {
    pub haplotypes: Arc<HaplotypeScores>,
    pub records: IndelRecords,
}

impl Scored {
    /// The column at `pos`'s records over `haplotypes`, which may have been
    /// scored together with neighbouring columns' alleles.
    pub fn new(haplotypes: Arc<HaplotypeScores>, pos: u32) -> Self {
        let records = call_on_pair_hmm(&haplotypes, pos);
        Self { haplotypes, records }
    }
}
