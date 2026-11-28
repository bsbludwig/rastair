//! Haplotype tags: which haplotype of which phase set a read pair comes from,
//! written as `HP` (1 or 2) and `PS`, the tags `whatshap haplotag` writes.
//!
//! `call` phases heterozygous records into blocks (`GT` `0|1`, `PS`); a
//! [`PhasedSite`] is one such record. `HP` 1 carries the allele left of the
//! `|`, and `PS` is the record's.
//!
//! **Votes.** A read votes at each phased site it shows one haplotype's allele
//! at and not the other's, with the base quality that tells them apart.
//! Conversion is allowed for, as in the calls: on OT a `T` may be a methylated
//! `C`, on OB an `A` a methylated `G`, so a read showing the converted base at
//! a `C/T` (or `G/A`) site votes only where its own neighbouring base is not
//! the CpG partner. The neighbour is read off the read, not the reference,
//! because the CpG may be one only this haplotype has.
//! An indel site is compared as sequence: what follows the anchor base in the
//! read, soft clip included, against what each haplotype reads there (its
//! allele, then [`INDEL_FLANK`] reference bases). The read votes for the
//! haplotype it matches throughout when it does not match the other. That
//! does not depend on where the aligner put the indel, and a read ending
//! inside a repeat, which fits both, casts no vote. It is strict on purpose:
//! a sequencing error or another variant in the window costs the vote, not
//! the tag's correctness.
//!
//! **Fragments.** Votes are pooled over the alignments sharing a read name
//! whose starts lie within [`MATE_REACH`] of the read being tagged, so mates
//! get one tag and a read over no site takes its mate's. The reach is measured
//! between two alignments, not from a segment's edge, which is what makes a
//! tag independent of how the run is cut into segments; `rewrite_region`
//! fetches that much beyond its segment for it. The phase set with the widest
//! margin between its two haplotypes wins (the lowest `PS` on a tie); a read
//! pair whose votes cancel gets no tag, and any `HP`/`PS` it came with is
//! removed.

use super::base_modification::{remove_aux_if_present, replace_aux};
use color_eyre::eyre::{Context as _, ContextCompat as _, Result};
use rust_htslib::{
    bam::{Record, ext::BamRecordExtensions as _, record::Aux},
    bcf::{self, record::GenotypeAllele},
};
use rustc_hash::FxHashMap;
use seqair_types::{Base, SmallVec, Strand, StrandFromRecord};
use std::collections::BTreeMap;

/// How far apart two alignments of one read name may start and still vote
/// together. Covers the insert of a short-read library; a mate further away
/// is tagged on its own.
pub(super) const MATE_REACH: u64 = 1000;

/// Reference bases after an indel record's `REF` that its haplotypes are
/// compared over: enough to leave the short repeats most indels sit in.
pub(super) const INDEL_FLANK: u64 = 30;

/// One of a phase set's two haplotypes, in `GT` order: [`Self::One`] carries
/// the allele left of the `|`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Haplotype {
    One,
    Two,
}

impl Haplotype {
    const fn index(self) -> usize {
        match self {
            Self::One => 0,
            Self::Two => 1,
        }
    }

    /// The `HP` value.
    const fn number(self) -> u8 {
        match self {
            Self::One => 1,
            Self::Two => 2,
        }
    }
}

/// What the two haplotypes carry at a phased site, [`Haplotype::One`] first.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum SiteAlleles {
    Snv([Base; 2]),
    /// What each haplotype reads after the anchor base: the rest of its
    /// allele, then the reference that follows `REF`.
    Indel([Vec<u8>; 2]),
}

/// A phased heterozygous record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct PhasedSite {
    /// 0-based; an indel's anchor base.
    pub pos: u32,
    /// The record's `PS`.
    pub set: i32,
    pub alleles: SiteAlleles,
}

impl PhasedSite {
    /// The phased site a VCF record is, if its first sample is a phased
    /// heterozygote with a `PS` and alleles spelled in bases.
    ///
    /// `reference_after` is asked for the bases from a 0-based position on,
    /// [`INDEL_FLANK`] of them or as many as the contig has.
    pub(super) fn from_record(
        record: &bcf::Record,
        reference_after: impl FnOnce(u64) -> Result<Vec<u8>>,
    ) -> Result<Option<Self>> {
        let Some([first, second]) = phased_alleles(record) else { return Ok(None) };
        let Some(set) = phase_set(record) else { return Ok(None) };

        let alleles = record.alleles();
        let (Some(reference), Some(first), Some(second)) =
            (alleles.first(), alleles.get(first), alleles.get(second))
        else {
            return Ok(None);
        };
        let spelled = |allele: &[u8]| {
            !allele.is_empty() && allele.iter().all(|base| b"ACGTacgt".contains(base))
        };
        if !(spelled(first) && spelled(second)) {
            return Ok(None);
        }
        let pos = u32::try_from(record.pos()).wrap_err("VCF position does not fit in u32")?;

        let alleles = match (*reference, *first, *second) {
            ([_], &[first], &[second]) => SiteAlleles::Snv([Base::from(first), Base::from(second)]),
            // Past the anchor the alleles are compared as sequence, which
            // needs them to start on the same base.
            (reference, [anchor, first @ ..], [second_anchor, second @ ..])
                if anchor.eq_ignore_ascii_case(second_anchor) =>
            {
                let after = u64::from(pos)
                    .checked_add(reference.len() as u64)
                    .wrap_err("VCF record ends past u64")?;
                let flank = reference_after(after)?;
                let haplotype = |tail: &[u8]| -> Vec<u8> {
                    tail.iter().chain(&flank).map(u8::to_ascii_uppercase).collect()
                };
                SiteAlleles::Indel([haplotype(first), haplotype(second)])
            }
            _ => return Ok(None),
        };
        Ok(Some(Self { pos, set, alleles }))
    }
}

/// The allele indices of a phased heterozygous diploid `GT`, in `GT` order.
fn phased_alleles(record: &bcf::Record) -> Option<[usize; 2]> {
    let genotypes = record.genotypes().ok()?;
    let genotype = genotypes.get(0);
    // Up to VCF 4.3 the first allele carries no phase bit, so the second
    // allele's is the one that says `|`.
    let [
        GenotypeAllele::Phased(first) | GenotypeAllele::Unphased(first),
        GenotypeAllele::Phased(second),
    ] = genotype[..]
    else {
        return None;
    };
    (first != second).then_some([usize::try_from(first).ok()?, usize::try_from(second).ok()?])
}

/// The first sample's `PS`, where it has one.
fn phase_set(record: &bcf::Record) -> Option<i32> {
    use bcf::record::Numeric as _;
    let sets = record.format(b"PS").integer().ok()?;
    let set = *sets.first()?.first()?;
    (!set.is_missing()).then_some(set)
}

/// A read's say at one phased site.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct Vote {
    pub set: i32,
    pub haplotype: Haplotype,
    /// Base quality of the base that tells the haplotypes apart.
    pub weight: u8,
}

type Votes = SmallVec<Vote, 4>;

/// Whether a haplotype carrying `allele` can show `observed` in a read of
/// `strand` (either, if unknown), whose bases before and after are `prev`
/// and `next`.
fn reads_as(allele: Base, observed: Base, strand: Option<Strand>, prev: Base, next: Base) -> bool {
    if allele == observed {
        return true;
    }
    let on = |converting: Strand| strand.is_none_or(|strand| strand == converting);
    match (allele, observed) {
        (Base::C, Base::T) => on(Strand::OT) && next == Base::G,
        (Base::G, Base::A) => on(Strand::OB) && prev == Base::C,
        _ => false,
    }
}

/// The haplotype a read base supports at an SNV site: the one whose allele
/// can read as it, when the other's cannot.
///
/// A neighbour the read does not show (its end, an `N`) may or may not be the
/// CpG partner, and the read votes only where that makes no difference.
pub(super) fn snv_vote(
    alleles: [Base; 2],
    strand: Option<Strand>,
    observed: Base,
    prev: Option<Base>,
    next: Option<Base>,
) -> Option<Haplotype> {
    let shown = |neighbour: Option<Base>| neighbour.filter(|base| *base != Base::Unknown);
    let either = |neighbour: Option<Base>, partner: Base| -> SmallVec<Base, 2> {
        match shown(neighbour) {
            Some(base) => SmallVec::from_iter([base]),
            None => SmallVec::from_iter([partner, Base::Unknown]),
        }
    };
    let (prevs, nexts) = (either(prev, Base::C), either(next, Base::G));

    let mut votes = prevs.iter().flat_map(|&prev| {
        nexts.iter().map(move |&next| {
            let [one, two] = alleles.map(|allele| reads_as(allele, observed, strand, prev, next));
            match (one, two) {
                (true, false) => Some(Haplotype::One),
                (false, true) => Some(Haplotype::Two),
                _ => None,
            }
        })
    });
    let first = votes.next()??;
    votes.all(|vote| vote == Some(first)).then_some(first)
}

/// Whether a haplotype base can show as a read base on `strand`, in any
/// context.
fn base_reads_as(haplotype: u8, read: u8, strand: Option<Strand>) -> bool {
    let on = |converting: Strand| strand.is_none_or(|strand| strand == converting);
    haplotype == read
        || (haplotype == b'C' && read == b'T' && on(Strand::OT))
        || (haplotype == b'G' && read == b'A' && on(Strand::OB))
}

/// The haplotype the bases after an indel's anchor support, with the read
/// index (into `read`) of the first base at which the haplotypes differ.
///
/// Compared over what the read and both haplotypes cover. The read has to
/// match one haplotype there base for base, conversions aside, and not the
/// other.
pub(super) fn indel_vote(
    haplotypes: &[Vec<u8>; 2],
    strand: Option<Strand>,
    read: &[u8],
) -> Option<(Haplotype, usize)> {
    let [one, two] = haplotypes;
    let covered = read.len().min(one.len()).min(two.len());
    let decisive = one.iter().zip(two).take(covered).position(|(one, two)| one != two)?;
    let matches = |haplotype: &[u8]| {
        haplotype
            .iter()
            .zip(read)
            .take(covered)
            .all(|(&haplotype, &read)| base_reads_as(haplotype, read, strand))
    };
    match (matches(one), matches(two)) {
        (true, false) => Some((Haplotype::One, decisive)),
        (false, true) => Some((Haplotype::Two, decisive)),
        _ => None,
    }
}

/// A read's votes at the phased sites it covers. `sites` are sorted by
/// position.
pub(super) fn votes_of(record: &Record, sites: &[PhasedSite]) -> Result<Votes> {
    let mut votes = Votes::new();
    let start = record.pos();
    let end = record.reference_end();
    let first = sites.partition_point(|site| i64::from(site.pos) < start);
    let covered = sites.get(first..).unwrap_or_default();
    let covered = covered
        .get(..covered.partition_point(|site| i64::from(site.pos) < end))
        .unwrap_or_default();
    if covered.is_empty() || record.is_unmapped() {
        return Ok(votes);
    }

    let strand = StrandFromRecord::strand(record);
    let seq = record.seq().as_bytes();
    let qual = record.qual();
    let base_at = |index: Option<usize>| index.and_then(|index| seq.get(index)).map(Base::from);

    let mut pending = covered.iter().peekable();
    for [read_pos, ref_pos] in record.aligned_pairs() {
        // A site the read deletes or skips has no base to vote with.
        while pending.next_if(|site| i64::from(site.pos) < ref_pos).is_some() {}
        let index = usize::try_from(read_pos).wrap_err("read position does not fit in usize")?;
        while let Some(site) = pending.next_if(|site| i64::from(site.pos) == ref_pos) {
            let vote = match &site.alleles {
                SiteAlleles::Snv(alleles) => base_at(Some(index))
                    .and_then(|observed| {
                        let prev = base_at(index.checked_sub(1));
                        let next = base_at(index.checked_add(1));
                        snv_vote(*alleles, strand, observed, prev, next)
                    })
                    .map(|haplotype| (haplotype, index)),
                SiteAlleles::Indel(haplotypes) => {
                    let after = index.saturating_add(1);
                    indel_vote(haplotypes, strand, seq.get(after..).unwrap_or_default())
                        .map(|(haplotype, decisive)| (haplotype, after.saturating_add(decisive)))
                }
            };
            if let Some((haplotype, decisive)) = vote
                && let Some(&weight) = qual.get(decisive)
            {
                votes.push(Vote { set: site.set, haplotype, weight });
            }
        }
        if pending.peek().is_none() {
            break;
        }
    }
    Ok(votes)
}

/// The haplotype of a phase set a read pair is tagged with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct HaplotypeTag {
    pub set: i32,
    pub haplotype: Haplotype,
}

/// The tag a fragment's votes add up to: the phase set whose haplotypes are
/// furthest apart, the lowest such set on a tie, and its heavier haplotype.
pub(super) fn assign(votes: impl IntoIterator<Item = Vote>) -> Option<HaplotypeTag> {
    let mut weights: BTreeMap<i32, [u32; 2]> = BTreeMap::new();
    for Vote { set, haplotype, weight } in votes {
        if let Some(sum) = weights.entry(set).or_default().get_mut(haplotype.index()) {
            *sum = sum.saturating_add(u32::from(weight));
        }
    }
    let mut best: Option<(u32, HaplotypeTag)> = None;
    for (set, [one, two]) in weights {
        let margin = one.abs_diff(two);
        if margin > best.map_or(0, |(margin, _)| margin) {
            let haplotype = if one > two { Haplotype::One } else { Haplotype::Two };
            best = Some((margin, HaplotypeTag { set, haplotype }));
        }
    }
    best.map(|(_, tag)| tag)
}

/// The votes of one alignment, and where it starts.
#[derive(Debug)]
struct AlignmentVotes {
    pos: i64,
    votes: Votes,
}

/// The votes of every alignment fetched for a segment, by read name.
#[derive(Debug, Default)]
pub(super) struct Fragments {
    by_name: FxHashMap<Vec<u8>, SmallVec<AlignmentVotes, 2>>,
}

impl Fragments {
    /// Record the votes `record` casts at `sites` (sorted by position).
    pub(super) fn add(&mut self, record: &Record, sites: &[PhasedSite]) -> Result<()> {
        let votes = votes_of(record, sites)?;
        if !votes.is_empty() {
            self.by_name
                .entry(record.qname().to_vec())
                .or_default()
                .push(AlignmentVotes { pos: record.pos(), votes });
        }
        Ok(())
    }

    /// The tag of `record`'s fragment: its own votes and those of the
    /// alignments of its read name starting within [`MATE_REACH`] of it.
    pub(super) fn tag_of(&self, record: &Record) -> Option<HaplotypeTag> {
        let pos = record.pos();
        let alignments = self.by_name.get(record.qname())?;
        assign(
            alignments
                .iter()
                .filter(|alignment| alignment.pos.abs_diff(pos) <= MATE_REACH)
                .flat_map(|alignment| alignment.votes.iter().copied()),
        )
    }

    /// Write `record`'s `HP` and `PS`, or remove them where its fragment has
    /// no tag.
    pub(super) fn tag(&self, record: &mut Record) -> Result<()> {
        match self.tag_of(record) {
            Some(HaplotypeTag { set, haplotype }) => {
                replace_aux(record, b"HP", Aux::U8(haplotype.number()))?;
                replace_aux(record, b"PS", Aux::I32(set))
            }
            None => {
                remove_aux_if_present(record, b"HP")?;
                remove_aux_if_present(record, b"PS")
            }
        }
    }
}

#[cfg(test)]
mod tests;
