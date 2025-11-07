use std::fmt;

use crate::{metrics2::PileupMetrics, vcf};

/// Struct holding references to current and its surrounding items
pub struct Surrounding<'a, T> {
    pub before: Option<&'a T>,
    pub current: &'a mut T,
    pub after: Option<&'a T>,
}

impl<'a, T: fmt::Debug> fmt::Debug for Surrounding<'a, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Surrounding")
            .field("before", &self.before)
            .field("current", &self.current)
            .field("after", &self.after)
            .finish()
    }
}

/// Get the pileup metrics for the record at `index` and its direct neighbors, if they exist
pub fn surrounding_pileups(
    records: &mut [PileupMetrics],
    index: usize,
) -> Surrounding<'_, PileupMetrics> {
    // To appease the borrow checker and get a mutable reference to the current record,
    // we split the records into three parts.
    let (left, right) = records.split_at_mut(index);
    let (current_slice, next_slice) = right.split_at_mut(1);
    let current = &mut current_slice[0];

    let before = left.last();
    let after = next_slice.first();
    // we might not have the direct neighbors
    let before = before.filter(|p| {
        p.pileup.contig() == current.pileup.contig()
            && Some(p.pileup.pos) == current.pileup.pos.checked_sub(1)
    });
    let after = after.filter(|p| {
        p.pileup.contig() == current.pileup.contig()
            && Some(p.pileup.pos) == current.pileup.pos.checked_add(1)
    });

    Surrounding { before, current, after }
}

/// Get the surrounding records for a given index in the records slice.
pub fn surrounding_records(
    records: &mut [vcf::Record],
    index: usize,
) -> (Option<&vcf::Record>, &mut vcf::Record, Option<&vcf::Record>) {
    // To appease the borrow checker and get a mutable reference to the current record,
    // we split the records into three parts.
    let (left, right) = records.split_at_mut(index);
    let (current_slice, next_slice) = right.split_at_mut(1);
    let current = &mut current_slice[0];

    let before = left.last();
    let after = next_slice.first();
    // we might not have the direct neighbors
    let before = before.filter(|r| {
        r.main.chrom == current.main.chrom && Some(r.main.pos) == current.main.pos.checked_sub(1)
    });
    let after = after.filter(|r| {
        r.main.chrom == current.main.chrom && Some(r.main.pos) == current.main.pos.checked_add(1)
    });

    (before, current, after)
}
