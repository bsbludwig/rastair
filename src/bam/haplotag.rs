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
//! between two alignments, not from a batch's edge, so a tag does not depend
//! on how the input is cut into batches. The phase set with the widest margin
//! between its two haplotypes wins (the lowest `PS` on a tie); a read pair
//! whose votes cancel gets no tag, and any `HP`/`PS` it came with is removed.
//!
//! **Two stages.** A worker gathers its batch's votes ([`BatchVotes`]) before
//! it rewrites the records, since `standard` mode rewrites the bases the votes
//! read. The mate may be in another batch, so the tags are set where the
//! batches come back in order: the [`Tagger`] holds a batch until the input
//! has moved more than [`MATE_REACH`] past it, then tags and hands it on. That
//! reads the input's order as coordinate order. On input in another order,
//! and past [`MAX_POOLED_RECORDS`] reads within the reach of one another
//! (about 7,500x of 150-base reads), mates pool only if they are in one batch
//! or in batches held together.
//! With `--region`, a mate outside the regions is not read and casts no vote.

use super::aligned_bases;
use color_eyre::eyre::{Context as _, ContextCompat as _, Result};
use rustc_hash::FxHashMap;
use seqair::{
    bam::{RawBamRecord, aligned_pairs::MatchPosition},
    reader::Tid,
};
use seqair_types::{Base, SmallVec, Strand};
use std::collections::{BTreeMap, VecDeque};

/// How far apart two alignments of one read name may start and still vote
/// together. Covers the insert of a short-read library; a mate further away
/// is tagged on its own.
pub(super) const MATE_REACH: u64 = 1000;

/// How many records the [`Tagger`] holds back, and how many written ones it
/// remembers the votes of, before it stops waiting for the input to move on.
const MAX_POOLED_RECORDS: usize = 50_000;

/// What a vote weighs in a read stored without base qualities, where every
/// base counts alike.
const NO_QUALITY_WEIGHT: u8 = 1;

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
    /// The site of a phased heterozygous record at the 0-based `pos` whose
    /// haplotypes carry `first` and `second`, if both are spelled in bases.
    ///
    /// `reference_after` is asked for the bases from a 0-based position on,
    /// [`INDEL_FLANK`] of them or as many as the contig has.
    pub(super) fn new(
        pos: u32,
        set: i32,
        reference: &[u8],
        [first, second]: [&[u8]; 2],
        reference_after: impl FnOnce(u64) -> Result<Vec<u8>>,
    ) -> Result<Option<Self>> {
        let spelled = |allele: &[u8]| {
            !allele.is_empty() && allele.iter().all(|base| b"ACGTacgt".contains(base))
        };
        if !(spelled(first) && spelled(second)) {
            return Ok(None);
        }

        let alleles = match (reference, first, second) {
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
pub(super) fn votes_of(record: &RawBamRecord, sites: &[PhasedSite]) -> Votes {
    let mut votes = Votes::new();
    let (Some(start), Some(last)) = (record.pos(), record.end_pos()) else { return votes };
    let flags = record.flags();
    let first = sites.partition_point(|site| site.pos < start.as_u32());
    let covered = sites.get(first..).unwrap_or_default();
    let covered = covered
        .get(..covered.partition_point(|site| site.pos <= last.as_u32()))
        .unwrap_or_default();
    if covered.is_empty() || flags.is_unmapped() {
        return votes;
    }

    let strand = Strand::from_flags(flags);
    let mut seq = Vec::with_capacity(record.seq_len());
    record.decode_seq_into(&mut seq);
    let qual = record.qual();
    let base_at = |index: Option<usize>| index.and_then(|index| seq.get(index)).map(Base::from);

    let mut pending = covered.iter().peekable();
    for MatchPosition { qpos, rpos, .. } in aligned_bases(record) {
        let (index, ref_pos) = (qpos.as_usize(), rpos.as_u32());
        // A site the read deletes or skips has no base to vote with.
        while pending.next_if(|site| site.pos < ref_pos).is_some() {}
        while let Some(site) = pending.next_if(|site| site.pos == ref_pos) {
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
            let Some((haplotype, decisive)) = vote else { continue };
            let weight = match qual {
                Some(qual) => qual.get(decisive).and_then(|quality| quality.get()),
                None => Some(NO_QUALITY_WEIGHT),
            };
            if let Some(weight) = weight {
                votes.push(Vote { set: site.set, haplotype, weight });
            }
        }
        if pending.peek().is_none() {
            break;
        }
    }
    votes
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

/// Where an alignment starts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Place {
    tid: Tid,
    pos: u64,
}

impl Place {
    fn of(record: &RawBamRecord) -> Option<Self> {
        Some(Self { tid: record.tid()?, pos: record.pos()?.as_u64() })
    }

    /// Whether alignments starting here and at `other` vote together.
    fn reaches(self, other: Self) -> bool {
        self.tid == other.tid && self.pos.abs_diff(other.pos) <= MATE_REACH
    }

    /// Whether a coordinate-sorted input that has come to `later` can still
    /// hold an alignment voting with one starting here.
    fn in_reach_of(self, later: Self) -> bool {
        self.tid == later.tid && later.pos <= self.pos.saturating_add(MATE_REACH)
    }
}

/// The votes of one alignment.
#[derive(Debug)]
struct AlignmentVotes {
    /// The batch it came in, by [`Tagger`]'s count.
    batch: u64,
    place: Place,
    votes: Votes,
}

/// The votes of a batch's records, by read name, gathered before the records
/// are rewritten.
#[derive(Debug, Default)]
pub(super) struct BatchVotes {
    alignments: Vec<(Box<[u8]>, Place, Votes)>,
}

impl BatchVotes {
    /// Record the votes `record` casts at `sites` (sorted by position).
    pub(super) fn add(&mut self, record: &RawBamRecord, sites: &[PhasedSite]) {
        let votes = votes_of(record, sites);
        if let (false, Some(place)) = (votes.is_empty(), Place::of(record)) {
            self.alignments.push((record.qname().into(), place, votes));
        }
    }
}

/// A batch the [`Tagger`] has seen.
#[derive(Debug)]
struct Seen<T> {
    number: u64,
    /// Until it is tagged and handed on.
    held: Option<(Vec<RawBamRecord>, T)>,
    records: usize,
    /// The read names it has votes under.
    names: Vec<Box<[u8]>>,
    /// Where its first and last placed records start.
    ends: Option<(Place, Place)>,
}

/// Tags batches arriving in input order, each once the batches after it have
/// brought the votes its reads' mates cast. `T` rides along with a batch.
#[derive(Debug)]
pub(super) struct Tagger<T> {
    /// By read name, the votes of the alignments of `window`.
    votes: FxHashMap<Box<[u8]>, SmallVec<AlignmentVotes, 2>>,
    /// In input order: batches whose votes are no longer needed leave at the
    /// front, and the held ones are at the back.
    window: VecDeque<Seen<T>>,
    batches: u64,
}

impl<T> Default for Tagger<T> {
    fn default() -> Self {
        Self { votes: FxHashMap::default(), window: VecDeque::new(), batches: 0 }
    }
}

impl<T> Tagger<T> {
    /// Take the next batch and return the batches that are now tagged, in
    /// input order.
    pub(super) fn push(
        &mut self,
        records: Vec<RawBamRecord>,
        votes: BatchVotes,
        with: T,
    ) -> Result<Vec<(Vec<RawBamRecord>, T)>> {
        self.batches += 1;
        let number = self.batches;
        let mut places = records.iter().filter_map(Place::of);
        let first = places.next();
        let ends = first.map(|first| (first, places.next_back().unwrap_or(first)));
        let mut names = Vec::with_capacity(votes.alignments.len());
        for (name, place, votes) in votes.alignments {
            self.votes.entry(name.clone()).or_default().push(AlignmentVotes {
                batch: number,
                place,
                votes,
            });
            names.push(name);
        }
        self.window.push_back(Seen {
            number,
            records: records.len(),
            held: Some((records, with)),
            names,
            ends,
        });

        let mut tagged = Vec::new();
        let mut held: usize =
            self.window.iter().filter(|seen| seen.held.is_some()).map(|seen| seen.records).sum();
        for index in 0..self.window.len() {
            let Some(seen) = self.window.get(index) else { break };
            if seen.held.is_none() {
                continue;
            }
            // The input has moved on when this batch's first record is out of
            // reach of that one's last.
            let waits = seen.ends.zip(first).is_some_and(|((_, last), now)| last.in_reach_of(now));
            if (waits || seen.number == number) && held <= MAX_POOLED_RECORDS {
                break;
            }
            held = held.saturating_sub(seen.records);
            tagged.extend(self.release(index)?);
        }
        self.forget();
        Ok(tagged)
    }

    /// The batches still held, tagged, once the input has ended.
    pub(super) fn finish(&mut self) -> Result<Vec<(Vec<RawBamRecord>, T)>> {
        let mut tagged = Vec::new();
        for index in 0..self.window.len() {
            tagged.extend(self.release(index)?);
        }
        self.window.clear();
        self.votes.clear();
        Ok(tagged)
    }

    /// Tag the batch at `index` of the window, if it is still held.
    fn release(&mut self, index: usize) -> Result<Option<(Vec<RawBamRecord>, T)>> {
        let Some((mut records, with)) =
            self.window.get_mut(index).and_then(|seen| seen.held.take())
        else {
            return Ok(None);
        };
        for record in &mut records {
            let tag = self.tag_of(record);
            write_tag(record, tag).wrap_err_with(|| {
                format!(
                    "failed to tag the haplotype of record {}",
                    String::from_utf8_lossy(record.qname())
                )
            })?;
        }
        Ok(Some((records, with)))
    }

    /// Drop the votes of batches at the front that are handed on and out of
    /// reach of everything still held or still to come.
    fn forget(&mut self) {
        let next = self
            .window
            .iter()
            .find(|seen| seen.held.is_some())
            .and_then(|seen| seen.ends)
            .map(|(first, _)| first)
            .or_else(|| self.window.back().and_then(|seen| seen.ends).map(|(_, last)| last));
        let mut remembered: usize =
            self.window.iter().filter(|seen| seen.held.is_none()).map(|seen| seen.records).sum();
        while let Some(front) = self.window.front() {
            let needed =
                front.ends.zip(next).is_some_and(|((_, last), next)| last.in_reach_of(next));
            if front.held.is_some() || (needed && remembered <= MAX_POOLED_RECORDS) {
                break;
            }
            remembered = remembered.saturating_sub(front.records);
            let Some(front) = self.window.pop_front() else { break };
            for name in front.names {
                if let Some(alignments) = self.votes.get_mut(&name) {
                    alignments.retain(|alignment| alignment.batch != front.number);
                    if alignments.is_empty() {
                        self.votes.remove(&name);
                    }
                }
            }
        }
    }

    /// The tag of `record`'s fragment: its own votes and those of the
    /// alignments of its read name starting within [`MATE_REACH`] of it.
    fn tag_of(&self, record: &RawBamRecord) -> Option<HaplotypeTag> {
        let place = Place::of(record)?;
        let alignments = self.votes.get(record.qname())?;
        assign(
            alignments
                .iter()
                .filter(|alignment| alignment.place.reaches(place))
                .flat_map(|alignment| alignment.votes.iter().copied()),
        )
    }
}

/// Write `record`'s `HP` and `PS`, or remove them where its fragment has no
/// tag.
fn write_tag(record: &mut RawBamRecord, tag: Option<HaplotypeTag>) -> Result<()> {
    let Some(HaplotypeTag { set, haplotype }) = tag else {
        // Most reads carry neither, and an edit copies the aux block; the
        // scan still fails on malformed tags, as an edit would.
        let mut tags = record.aux().iter_checked();
        let present = tags.try_fold(false, |present, tag| {
            tag.map(|(tag, _)| present || tag == *b"HP" || tag == *b"PS")
        });
        if !present.wrap_err("could not read the tags")? {
            return Ok(());
        }
        return record
            .edit_aux(|aux| {
                aux.remove(*b"HP");
                aux.remove(*b"PS");
                Ok(())
            })
            .wrap_err("could not remove the haplotype tags");
    };
    record
        .edit_aux(|aux| {
            aux.set_int(*b"HP", i64::from(haplotype.number()))?;
            aux.set_int(*b"PS", i64::from(set))
        })
        .wrap_err("could not write the haplotype tags")
}

#[cfg(test)]
mod tests;
