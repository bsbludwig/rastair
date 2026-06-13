#[cfg(not(feature = "experimental-seqair"))]
use super::{
    indels::{IndelAllele, IndelObservation},
    overlapping_reads::{NameCollector, resolve_pair},
    ref_features::{dinucleotide_run_at, homopolymer_run_at, indel_ref_window_at},
};
use crate::sequence::Segment;
#[cfg(not(feature = "experimental-seqair"))]
use crate::{
    call::{
        pileup::{Pileup, PositionInRead, SimpleRead, SimpleReads},
        process::PileupMappingParams,
    },
    utils::SequenceContext,
};
#[cfg(not(feature = "experimental-seqair"))]
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
#[cfg(not(feature = "experimental-seqair"))]
use rust_htslib::bam::{
    pileup::{Alignment, Indel, Pileup as HtsPileup},
    record::RecordView,
};
use rust_htslib::bam::{Record, ext::BamRecordExtensions as _};
#[cfg(not(feature = "experimental-seqair"))]
use rustc_hash::FxHashMap;
use rustc_hash::FxHasher;
#[cfg(not(feature = "experimental-seqair"))]
use seqair_types::SmallVec;
#[cfg(not(feature = "experimental-seqair"))]
use seqair_types::strand_from_flags;
use seqair_types::{Base, Strand};
use std::hash::{Hash, Hasher};
#[cfg(not(feature = "experimental-seqair"))]
use std::{mem, rc::Rc};
use tracing::trace;
#[cfg(not(feature = "experimental-seqair"))]
use tracing::{debug, instrument};

#[cfg(not(feature = "experimental-seqair"))]
#[derive(Default)]
pub(crate) struct ReadOrientationCache {
    strands: FxHashMap<ReadOrientationCacheKey, Strand>,
}

#[cfg(not(feature = "experimental-seqair"))]
impl ReadOrientationCache {
    fn strand_for_alignment(
        &mut self,
        alignment: &Alignment<'_>,
        segment: &Segment,
        params: &PileupMappingParams,
    ) -> Option<Strand> {
        let record = alignment.record_view();
        if !params.guess_read_orientation {
            return strand_from_flags(record.flags().into()).ok();
        }

        let key = ReadOrientationCacheKey::from_alignment(alignment);
        if let Some(&strand) = self.strands.get(&key) {
            return Some(strand);
        }

        let strand = infer_strand_from_mismatch_motifs(&alignment.record(), segment);
        self.strands.insert(key, strand);
        Some(strand)
    }
}

#[cfg(not(feature = "experimental-seqair"))]
#[derive(Default)]
pub(crate) struct ReadMismatchCache {
    counts: FxHashMap<ReadOrientationCacheKey, u32>,
}

#[cfg(not(feature = "experimental-seqair"))]
impl ReadMismatchCache {
    fn mismatch_count_for_alignment(
        &mut self,
        alignment: &Alignment<'_>,
        segment: &Segment,
        strand: Strand,
    ) -> u32 {
        let key = ReadOrientationCacheKey::from_alignment(alignment);
        if let Some(&count) = self.counts.get(&key) {
            return count;
        }
        let count = count_taps_aware_mismatches(&alignment.record(), segment, strand);
        self.counts.insert(key, count);
        count
    }
}

#[cfg(not(feature = "experimental-seqair"))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct ReadOrientationCacheKey {
    name_hash: u64,
    start: i64,
    flags: u16,
}

#[cfg(not(feature = "experimental-seqair"))]
impl ReadOrientationCacheKey {
    fn from_alignment(alignment: &Alignment<'_>) -> Self {
        let record = alignment.record_view();
        Self {
            name_hash: hash_bytes(record.qname()),
            start: alignment.record().pos(),
            flags: record.flags(),
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct ReadOrientationEvidence {
    tg: u32,
    ca: u32,
}

impl ReadOrientationEvidence {
    fn ot_likelihood(self) -> Option<f64> {
        let total = self.tg + self.ca;
        (total > 0).then(|| f64::from(self.tg) / f64::from(total))
    }
}

#[cfg(not(feature = "experimental-seqair"))]
/// Buffers and caches reused across every pileup column in a segment.
///
/// A whole-genome run calls [`Pileup::from_hts`] once per reference base, so
/// anything allocated per call is allocated billions of times; these clear and
/// refill instead.
pub(crate) struct PileupScratch {
    names: NameCollector,
    orientation: ReadOrientationCache,
    mismatches: ReadMismatchCache,
    fragments: FragmentVotes,
}

#[cfg(not(feature = "experimental-seqair"))]
impl PileupScratch {
    pub(crate) fn new(params: &PileupMappingParams) -> Self {
        Self {
            names: NameCollector::new(params),
            orientation: ReadOrientationCache::default(),
            mismatches: ReadMismatchCache::default(),
            fragments: FragmentVotes::default(),
        }
    }
}

#[cfg(not(feature = "experimental-seqair"))]
/// What each fragment contributes to one pileup column.
///
/// Entries are indexed by the slot of whichever mate the pileup reported first,
/// so both mates of an overlapping pair land on the same vote. Keeping separate
/// tallies instead lets a fragment whose mates disagree — one spanning the
/// indel, the other soft-clipped over it — be counted on the alternate side
/// *and* on the noisy-reference side, and [`IndelCounts::clean_depth`] then
/// subtracts it twice.
///
/// [`IndelCounts::clean_depth`]: super::indels::IndelCounts::clean_depth
#[derive(Default)]
struct FragmentVotes {
    votes: Vec<FragmentVote>,
}

#[cfg(not(feature = "experimental-seqair"))]
#[derive(Clone, Copy, Default)]
struct FragmentVote {
    /// Some mate of this fragment carries an indel here.
    indel: bool,
    /// Some mate supports the reference from an alignment of the kind that slips.
    noisy_ref: bool,
    soft_clipped: bool,
}

#[cfg(not(feature = "experimental-seqair"))]
#[derive(Default)]
struct FragmentTotals {
    noisy_ref_count: u32,
    soft_clip_count: u32,
}

#[cfg(not(feature = "experimental-seqair"))]
impl FragmentVotes {
    fn prepare(&mut self, capacity: usize) {
        self.votes.clear();
        self.votes.reserve(capacity);
    }

    fn slot(&mut self, fragment: usize) -> Option<&mut FragmentVote> {
        if fragment >= self.votes.len() {
            self.votes.resize(fragment + 1, FragmentVote::default());
        }
        self.votes.get_mut(fragment)
    }

    /// Record an indel for this fragment, reporting whether it is the first one:
    /// a pair whose mates both span the indel casts one vote, not two.
    fn saw_indel(&mut self, fragment: usize) -> bool {
        self.slot(fragment).is_some_and(|vote| !mem::replace(&mut vote.indel, true))
    }

    fn saw_reference(&mut self, fragment: usize, shape: AlignmentShape) {
        if let Some(vote) = self.slot(fragment) {
            vote.noisy_ref |= shape.noisy();
        }
    }

    fn saw_soft_clip(&mut self, fragment: usize) {
        if let Some(vote) = self.slot(fragment) {
            vote.soft_clipped = true;
        }
    }

    /// An indel-carrying fragment is already off the reference side, because
    /// `ref_count` is `reads.len() - total_indel_reads`, so it must not also
    /// appear in `noisy_ref_count`.
    fn totals(&self) -> FragmentTotals {
        let mut totals = FragmentTotals::default();
        for vote in &self.votes {
            if vote.soft_clipped {
                totals.soft_clip_count += 1;
            }
            if vote.indel {
                continue;
            }
            if vote.noisy_ref {
                totals.noisy_ref_count += 1;
            }
        }
        totals
    }
}

#[cfg(not(feature = "experimental-seqair"))]
/// The alignment shapes that make an indel call at this column unreliable.
#[derive(Clone, Copy)]
struct AlignmentShape {
    terminal_repeat: bool,
    soft_clipped: bool,
}

#[cfg(not(feature = "experimental-seqair"))]
impl AlignmentShape {
    fn of(record: &RecordView<'_>) -> Self {
        let (seq, _) = record.seq_and_qual();
        Self {
            terminal_repeat: has_repeat_seq(&seq, 1, HOMOPOLYMER_UNITS)
                || has_repeat_seq(&seq, 2, DINUCLEOTIDE_UNITS),
            soft_clipped: has_soft_clip(record.raw_cigar()),
        }
    }

    fn noisy(self) -> bool {
        self.terminal_repeat || self.soft_clipped
    }
}

#[cfg(not(feature = "experimental-seqair"))]
impl Pileup {
    /// Convert a pileup from htslib into our internal Pileup representation.
    #[instrument(level = "trace", skip_all)]
    pub(crate) fn from_hts(
        pile: &HtsPileup,
        segment: Rc<Segment>,
        params: &PileupMappingParams,
        scratch: &mut PileupScratch,
    ) -> Result<Pileup> {
        let pos = pile.pos();
        let idx = segment.pos_to_idx(pos)?;
        let depth = pile.depth();
        let max_reads = depth.min(params.max_coverage);
        if depth > max_reads {
            debug!(pos, depth, "Capping number of reads in pileup to {max_reads}");
        }
        let max_reads = usize::try_from(max_reads).wrap_err("max_reads exceeds usize")?;

        let mut raw_reads = Vec::with_capacity(max_reads);
        let mut to_remove = SmallVec::<usize, 16>::new();
        let mut indel_observations = SmallVec::new();

        scratch.fragments.prepare(max_reads);
        if let NameCollector::Collect(buf) = &mut scratch.names {
            buf.prepare(max_reads);
        }

        // Both sides of `VAF = alt / depth` have to be drawn from the same reads,
        // so indel observations are collected in the same pass, behind the same
        // filters, as the reads that make up the depth. Collecting them separately
        // let an alignment that `alignment_to_read` rejects — MAPQ 0, a base below
        // `--min-baseq`, the interior of a deletion, anything past `--max-coverage`
        // — still cast an indel vote while contributing nothing to the depth, which
        // floors `ref_count` to zero and reads VAF as 1.0.
        //
        // NOTE: flag and read-group filtering already happened in the pileup-level
        // filter, so only read masking and quality are applied here.
        for a in pile.alignments() {
            if raw_reads.len() >= max_reads {
                break;
            }
            let Some((name, read)) =
                alignment_to_read(&a, segment.as_ref(), params, &mut scratch.orientation)
            else {
                continue;
            };
            if !params.read_masking.filter(&read) || !params.quality.filter(&read) {
                continue;
            }

            // Indel accounting is per *fragment*, matching the deduplication that
            // `reads` gets below: counting observations per alignment while the
            // depth is per fragment double-counts overlapping mates. With
            // `--keep-overlapping-reads` nothing is deduplicated and both sides
            // stay at the alignment level, which is equally consistent.
            let slot = raw_reads.len();
            let duplicate = match &mut scratch.names {
                NameCollector::Skip => None,
                NameCollector::Collect(buf) => buf.see(name, slot),
            };
            let fragment = duplicate.unwrap_or(slot);

            let shape = AlignmentShape::of(&a.record_view());
            if shape.soft_clipped {
                scratch.fragments.saw_soft_clip(fragment);
            }
            match a.indel() {
                Indel::None => scratch.fragments.saw_reference(fragment, shape),
                _ => {
                    let observation = indel_observation(
                        &a,
                        &read,
                        pos,
                        segment.as_ref(),
                        params,
                        &mut scratch.mismatches,
                        shape,
                    );
                    if let Some(observation) = observation
                        && scratch.fragments.saw_indel(fragment)
                    {
                        indel_observations.push(observation);
                    }
                }
            }

            raw_reads.push(read);
            if let Some(first) = duplicate {
                resolve_pair(&raw_reads, slot, first, &mut to_remove);
            }
        }

        // Three-plus alignments sharing a qname can nominate the same slot twice,
        // and swap_remove of an already-removed slot deletes an unrelated read.
        to_remove.sort_unstable();
        to_remove.dedup();
        for &slot in to_remove.iter().rev() {
            raw_reads.swap_remove(slot);
        }
        let reads = SimpleReads(raw_reads.into());

        let FragmentTotals { noisy_ref_count, soft_clip_count } = scratch.fragments.totals();

        let reference_base: Base =
            segment.sequence.get(idx).wrap_err("failed to get reference base")?.into();

        let context =
            SequenceContext::new(idx, &segment).wrap_err("failed to get sequence context")?;

        // The reference window is only needed for indel slippage features, so
        // skip the allocation at the (vast majority of) positions without indels.
        let (indel_ref_window, indel_ref_anchor) = if indel_observations.is_empty() {
            (SmallVec::new(), 0)
        } else {
            indel_ref_window_at(idx, &segment)
        };

        let segment_start = segment.range.region.start as usize;

        Ok(Pileup {
            region: segment.range.clone(),
            context,
            pos,
            reads,
            reference_base,
            indel_observations,
            noisy_ref_count,
            // `pos` is the anchor, i.e. the base *before* the indel; for a
            // left-aligned indel that is the base before the repeat. Measured at the
            // anchor these read ~1 exactly where the tract is longest.
            homopolymer_run: homopolymer_run_at(pos as usize + 1, &segment, segment_start),
            dinucleotide_run: dinucleotide_run_at(pos as usize + 1, &segment, segment_start),
            soft_clip_count,
            indel_ref_window,
            indel_ref_anchor,
        })
    }
}

#[cfg(not(feature = "experimental-seqair"))]
/// The indel this alignment shows at `pos`, if it survives the indel-specific
/// filters.
///
/// Everything both an observation and a [`SimpleRead`] need is taken from `read`
/// rather than recomputed, so the CIGAR is walked once per alignment.
fn indel_observation(
    a: &Alignment<'_>,
    read: &SimpleRead,
    pos: u32,
    segment: &Segment,
    params: &PileupMappingParams,
    mismatch_cache: &mut ReadMismatchCache,
    shape: AlignmentShape,
) -> Option<IndelObservation> {
    let qpos = read.position.pos as usize;
    let read_len = read.position.read_length as usize;
    let cutoff = params.indel_end_of_read_cutoff;
    // End-of-read filter, stricter than the one SNVs get.
    if qpos < cutoff || qpos >= read_len.saturating_sub(cutoff) {
        return None;
    }

    let mismatches = mismatch_cache.mismatch_count_for_alignment(a, segment, read.strand);
    if mismatches > params.indel_max_mismatches {
        trace!(
            mismatches,
            max = params.indel_max_mismatches,
            "Indel skipped: too many non-TAPS mismatches"
        );
        return None;
    }

    let record = a.record_view();
    let (seq, qual) = record.seq_and_qual();
    let allele = match a.indel() {
        Indel::Ins(len) => {
            let start = qpos + 1;
            let end = start + len as usize;
            let bases: SmallVec<Base, 4> =
                (start..end).filter_map(|i| seq.get(i).map(Base::from)).collect();
            IndelAllele::Insertion(bases)
        }
        Indel::Del(len) => {
            let ref_start = (pos as usize + 1).saturating_sub(segment.range.region.start as usize);
            let ref_end = ref_start + len as usize;
            let bases: SmallVec<Base, 4> = segment
                .sequence
                .get(ref_start..ref_end)
                .map(|slice| slice.iter().copied().map(Base::from).collect())
                .unwrap_or_default();
            IndelAllele::Deletion(bases)
        }
        Indel::None => return None,
    };
    if allele.bases().is_empty() {
        return None;
    }

    let insertion_base_quals = match &allele {
        IndelAllele::Insertion(bases) => {
            let start = qpos + 1;
            (start..start + bases.len()).filter_map(|i| qual.get(i).copied()).collect()
        }
        IndelAllele::Deletion(_) => SmallVec::new(),
    };
    let post_del_base_qual = match &allele {
        IndelAllele::Deletion(_) => qual.get(qpos + 1).copied().unwrap_or(0),
        IndelAllele::Insertion(_) => 0,
    };

    Some(IndelObservation {
        allele,
        strand: read.strand,
        reverse: read.reverse,
        pos_in_read: read.position.pos,
        read_length: read.position.read_length,
        mapq: read.mapq,
        base_qual: read.qual,
        matching_bases: read.matching_bases,
        num_indels_in_read: read.indels,
        insertion_base_quals,
        post_del_base_qual,
        has_repeat: shape.terminal_repeat,
        noisy: shape.noisy(),
    })
}

#[cfg(not(feature = "experimental-seqair"))]
fn alignment_to_read<'a>(
    a: &Alignment<'a>,
    segment: &Segment,
    params: &PileupMappingParams,
    orientation_cache: &mut ReadOrientationCache,
) -> Option<(&'a [u8], SimpleRead)> {
    let pos = a.qpos()?;
    let record = a.record_view();
    let flags = record.flags();
    let strand = orientation_cache.strand_for_alignment(a, segment, params)?;
    let (matches, indels) = calc_cigar_data(record.raw_cigar());
    let (seq, qual) = record.seq_and_qual();

    let before_base = pos.checked_sub(1).map(|i| Base::from(seq[i]));
    let after_base = match a.indel() {
        Indel::None if pos + 1 < seq.len() => Some(Base::from(seq[pos + 1])),
        _ => None,
    };

    Some((
        record.qname(),
        SimpleRead {
            base: seq[pos].into(),
            qual: *qual.get(pos)?,
            mapq: record.mapq(),
            strand,
            reverse: flags & 0x10 != 0,
            second: flags & 0x80 != 0,
            position: PositionInRead {
                pos: u32::try_from(pos).expect("position fits in u32"),
                read_length: u32::try_from(seq.len()).expect("read length fits in u32"),
            },
            matching_bases: matches,
            indels,
            before_base,
            after_base,
        },
    ))
}

pub(crate) fn infer_strand_from_mismatch_motifs(record: &Record, segment: &Segment) -> Strand {
    let evidence = collect_mismatch_motif_evidence(record, segment);
    let likelihood = evidence.ot_likelihood();
    let strand = match evidence.tg.cmp(&evidence.ca) {
        std::cmp::Ordering::Greater => Strand::OT,
        std::cmp::Ordering::Less => Strand::OB,
        std::cmp::Ordering::Equal => pseudo_random_strand(record),
    };

    trace!(
        read_id = %String::from_utf8_lossy(record.qname()),
        tg_evidence = evidence.tg,
        ca_evidence = evidence.ca,
        ot_likelihood = likelihood,
        assigned = %strand,
        "Assigned read orientation from mismatch motifs",
    );

    strand
}

fn collect_mismatch_motif_evidence(record: &Record, segment: &Segment) -> ReadOrientationEvidence {
    let seq = record.seq().as_bytes();
    let mut evidence = ReadOrientationEvidence::default();

    for [pos_in_read, pos_in_ref] in record.aligned_pairs_full() {
        let (Some(pos_in_read), Some(pos_in_ref)) = (pos_in_read, pos_in_ref) else {
            continue;
        };
        let Some(read_idx) = usize::try_from(pos_in_read).ok() else {
            continue;
        };
        let Some(observed) = seq.get(read_idx).copied().map(Base::from) else {
            continue;
        };
        let Some(reference) = reference_base(segment, pos_in_ref) else {
            continue;
        };

        if observed == Base::Unknown || reference == Base::Unknown || observed == reference {
            continue;
        }

        add_motif_evidence(&seq, read_idx, &mut evidence);
    }

    evidence
}

fn add_motif_evidence(seq: &[u8], read_idx: usize, evidence: &mut ReadOrientationEvidence) {
    if let Some(motif) = dinucleotide(seq, read_idx, read_idx.checked_add(1)) {
        count_motif(motif, evidence);
    }
    if let Some(previous_idx) = read_idx.checked_sub(1)
        && let Some(motif) = dinucleotide(seq, previous_idx, Some(read_idx))
    {
        count_motif(motif, evidence);
    }
}

fn dinucleotide(seq: &[u8], first: usize, second: Option<usize>) -> Option<[Base; 2]> {
    let second = second?;
    Some([Base::from(*seq.get(first)?), Base::from(*seq.get(second)?)])
}

fn count_motif(motif: [Base; 2], evidence: &mut ReadOrientationEvidence) {
    match motif {
        [Base::T, Base::G] => evidence.tg += 1,
        [Base::C, Base::A] => evidence.ca += 1,
        _ => {}
    }
}

#[cfg(not(feature = "experimental-seqair"))]
fn count_taps_aware_mismatches(record: &Record, segment: &Segment, strand: Strand) -> u32 {
    let seq = record.seq().as_bytes();
    let mut count = 0u32;

    for [pos_in_read, pos_in_ref] in record.aligned_pairs_full() {
        let (Some(pos_in_read), Some(pos_in_ref)) = (pos_in_read, pos_in_ref) else {
            continue; // insertion or deletion gap — not a mismatch
        };
        let Some(read_idx) = usize::try_from(pos_in_read).ok() else { continue };
        let Some(observed) = seq.get(read_idx).copied().map(Base::from) else { continue };
        let Some(reference) = reference_base(segment, pos_in_ref) else { continue };

        if observed == Base::Unknown || reference == Base::Unknown || observed == reference {
            continue;
        }

        // Skip mismatches that are expected TAPS signal, not sequencing errors.
        // For Unknown strand, conservatively skip both patterns to avoid penalising
        // methylated reads whose strand couldn't be determined.
        let is_taps_signal = match strand {
            Strand::OT => observed == Base::T && reference == Base::C,
            Strand::OB => observed == Base::A && reference == Base::G,
            Strand::Unknown => {
                (observed == Base::T && reference == Base::C)
                    || (observed == Base::A && reference == Base::G)
            }
        };

        if !is_taps_signal {
            count += 1;
        }
    }

    count
}

fn reference_base(segment: &Segment, pos_in_ref: i64) -> Option<Base> {
    let pos_in_ref = u32::try_from(pos_in_ref).ok()?;
    let idx = segment.pos_to_idx(pos_in_ref).ok()?;
    segment.sequence.get(idx).copied().map(Base::from)
}

fn pseudo_random_strand(record: &Record) -> Strand {
    let mut hasher = FxHasher::default();
    record.qname().hash(&mut hasher);
    record.pos().hash(&mut hasher);
    record.flags().hash(&mut hasher);
    if hasher.finish() & 1 == 0 { Strand::OT } else { Strand::OB }
}

#[cfg(not(feature = "experimental-seqair"))]
fn hash_bytes(bytes: &[u8]) -> u64 {
    let mut hasher = FxHasher::default();
    bytes.hash(&mut hasher);
    hasher.finish()
}

#[cfg(not(feature = "experimental-seqair"))]
/// Calculate matches and indels from a packed CIGAR array.
///
/// Lower 4 bits encode the operation; upper 28 bits encode the length.
fn calc_cigar_data(cigar: &[u32]) -> (u32, u32) {
    let mut matches = 0;
    let mut indels = 0;
    for c in cigar {
        let len = c >> 4;
        match c & 0b1111 {
            0 => matches += len,
            1 | 2 => indels += len,
            _ => {}
        }
    }
    (matches, indels)
}

#[cfg(not(feature = "experimental-seqair"))]
/// Check if a CIGAR array contains a soft-clip operation (op 4).
fn has_soft_clip(cigar: &[u32]) -> bool {
    cigar.iter().any(|&c| c & 0xF == 4)
}

#[cfg(not(feature = "experimental-seqair"))]
/// Check if first or last `cutoff` bases of a read form a repeating pattern of length `n`.
/// Whether either read terminus is a tandem repeat of period `n`, measured in whole
/// repeat units.
///
/// Units rather than a shared base window: with a 3 bp window the period-2 arm
/// reduces to `seq[0] == seq[2] || seq[len-3] == seq[len-1]`, true for 43.75% of
/// random reads, which makes the noise flag fire on a typical read rather than an
/// unusual one. At 4 units a terminal homopolymer occurs ~3% of the time and a
/// 3-unit dinucleotide repeat ~0.8%.
/// A flagged read should be unusual, not typical: see [`has_repeat_seq`].
const HOMOPOLYMER_UNITS: usize = 4;
const DINUCLEOTIDE_UNITS: usize = 3;

fn has_repeat_seq(seq: &rust_htslib::bam::record::Seq<'_>, n: usize, units: usize) -> bool {
    let len = seq.len();
    let Some(window) = n.checked_mul(units).filter(|w| *w <= len) else {
        return false;
    };
    if n == 0 || units < 2 {
        return false;
    }

    let periodic = |start: usize| {
        (start..start + window - n).all(|i| match (seq.get(i), seq.get(i + n)) {
            (Some(a), Some(b)) => a == b,
            _ => false,
        })
    };

    periodic(0) || periodic(len - window)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(not(feature = "experimental-seqair"))]
    use crate::call::{process::PileupMappingParams, variant_calling::VariantCallingParams};
    use crate::sequence::{ChunkRegion, Region, Segment};
    #[cfg(not(feature = "experimental-seqair"))]
    use crate::utils::default;
    use rust_htslib::bam::record::{Cigar, CigarString};

    /// A run counter is a `u8` but a reference tract can be longer than 255. The
    /// wrapped value would land back in the short-tract rows of
    /// `ALLELE_FRACTION_BY_TRACT`, applying the simple-sequence prior inside the
    /// longest repeat in the genome — and panic outright in a debug build.
    #[test]
    fn a_reference_run_longer_than_the_counter_saturates() {
        let homopolymer = test_segment(&[b'A'; 400]);
        assert_eq!(homopolymer_run_at(300, &homopolymer, 100), u8::MAX);

        let dinucleotide = test_segment(&b"AT".repeat(200));
        assert!(dinucleotide_run_at(300, &dinucleotide, 100) >= u8::MAX - 1);
    }

    fn test_segment(sequence: &[u8]) -> Segment {
        let start = 100u64;
        let end = start + u64::try_from(sequence.len()).expect("sequence length fits") - 1;
        Segment {
            range: ChunkRegion {
                region: Region { contig: "chrTest".into(), start, end },
                last_position: end,
                overlap_start: 0,
                overlap_end: 0,
            },
            sequence: sequence.to_vec(),
            overlap_start: 0,
            overlap_end: 0,
        }
    }

    fn test_record(qname: &[u8], flags: u16, start: i64, sequence: &[u8]) -> Record {
        let mut record = Record::new();
        let cigar = CigarString(
            vec![Cigar::Match(u32::try_from(sequence.len()).expect("sequence length fits in u32"))]
                .into(),
        );
        record.set(qname, Some(&cigar), sequence, &vec![40; sequence.len()]);
        record.set_flags(flags);
        record.set_pos(start);
        record
    }

    #[cfg(not(feature = "experimental-seqair"))]
    #[test]
    fn name_collector_skip_when_keep_overlapping() {
        let params = PileupMappingParams {
            variant_calling: VariantCallingParams { keep_overlapping_reads: true, ..default() },
            ..default()
        };
        assert!(matches!(NameCollector::new(&params), NameCollector::Skip));
    }

    #[cfg(not(feature = "experimental-seqair"))]
    #[test]
    fn name_collector_collect_by_default() {
        let params = PileupMappingParams::default();
        assert!(matches!(NameCollector::new(&params), NameCollector::Collect(_)));
    }

    #[test]
    fn mismatch_motifs_forward_tg_is_ot() {
        let segment = test_segment(b"ACGT");
        let record = test_record(b"forward-tg", 99, 100, b"ATGT");

        assert_eq!(infer_strand_from_mismatch_motifs(&record, &segment), Strand::OT);
    }

    #[test]
    fn mismatch_motifs_forward_ca_is_ob() {
        let segment = test_segment(b"CCGT");
        let record = test_record(b"forward-ca", 99, 100, b"CAGT");

        assert_eq!(infer_strand_from_mismatch_motifs(&record, &segment), Strand::OB);
    }

    #[test]
    fn mismatch_motifs_reverse_ca_is_ob() {
        let segment = test_segment(b"CCGT");
        let record = test_record(b"reverse-ca", 83, 100, b"CAGT");

        assert_eq!(infer_strand_from_mismatch_motifs(&record, &segment), Strand::OB);
    }

    #[test]
    fn mismatch_motifs_reverse_tg_is_ot() {
        let segment = test_segment(b"AAGA");
        let record = test_record(b"reverse-tg", 83, 100, b"ATGA");

        assert_eq!(infer_strand_from_mismatch_motifs(&record, &segment), Strand::OT);
    }

    #[test]
    fn mismatch_motifs_use_next_base_window_for_denovo_like_tg() {
        let segment = test_segment(b"AAGA");
        let record = test_record(b"denovo-next", 99, 100, b"ATGA");

        assert_eq!(infer_strand_from_mismatch_motifs(&record, &segment), Strand::OT);
    }

    #[test]
    fn mismatch_motifs_use_previous_base_window_for_denovo_like_ca() {
        let segment = test_segment(b"CTGT");
        let record = test_record(b"denovo-prev", 99, 100, b"CAGT");

        assert_eq!(infer_strand_from_mismatch_motifs(&record, &segment), Strand::OB);
    }

    #[test]
    fn mismatch_motif_tie_break_is_reproducible() {
        let segment = test_segment(b"AAAA");
        let record = test_record(b"no-evidence", 99, 100, b"CCCC");

        let first = infer_strand_from_mismatch_motifs(&record, &segment);
        let second = infer_strand_from_mismatch_motifs(&record, &segment);

        assert_eq!(first, second);
        assert!(matches!(first, Strand::OT | Strand::OB));
    }

    #[cfg(not(feature = "experimental-seqair"))]
    #[test]
    fn taps_aware_ot_ct_not_counted() {
        // C→T on OT is TAPS methylation signal — must not count
        let segment = test_segment(b"ACGT");
        let record = test_record(b"ot-ct", 99, 100, b"ATGT"); // pos 1: C→T
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::OT), 0);
    }

    #[cfg(not(feature = "experimental-seqair"))]
    #[test]
    fn taps_aware_ob_ga_not_counted() {
        // G→A on OB is TAPS methylation signal — must not count
        let segment = test_segment(b"ACGT");
        let record = test_record(b"ob-ga", 83, 100, b"ACAT"); // pos 2: G→A
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::OB), 0);
    }

    #[cfg(not(feature = "experimental-seqair"))]
    #[test]
    fn taps_aware_ot_other_mismatch_counted() {
        // A→G on OT is a real sequencing error — must count
        let segment = test_segment(b"ACGT");
        let record = test_record(b"ot-ag", 99, 100, b"GCGT"); // pos 0: A→G
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::OT), 1);
    }

    #[cfg(not(feature = "experimental-seqair"))]
    #[test]
    fn taps_aware_ob_ct_counted() {
        // C→T on OB is not a TAPS signal (wrong strand) — must count
        let segment = test_segment(b"ACGT");
        let record = test_record(b"ob-ct", 83, 100, b"ATGT"); // pos 1: C→T on OB
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::OB), 1);
    }

    #[cfg(not(feature = "experimental-seqair"))]
    #[test]
    fn taps_aware_unknown_strand_excludes_both_patterns() {
        // Unknown strand: conservatively skip both C→T and G→A to avoid
        // penalising methylated reads whose strand couldn't be determined
        let segment = test_segment(b"CGCG");
        let record = test_record(b"unknown", 99, 100, b"TATA"); // C→T and G→A alternating
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::Unknown), 0);
    }

    #[cfg(not(feature = "experimental-seqair"))]
    #[test]
    fn taps_aware_no_mismatches_returns_zero() {
        let segment = test_segment(b"ACGT");
        let record = test_record(b"perfect", 99, 100, b"ACGT");
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::OT), 0);
    }

    #[cfg(not(feature = "experimental-seqair"))]
    #[test]
    fn taps_aware_mixed_counts_only_non_taps() {
        // A→G (real error, counted) + C→T (TAPS on OT, excluded) = 1
        let segment = test_segment(b"ACGT");
        let record = test_record(b"mixed", 99, 100, b"GTGT"); // pos 0: A→G, pos 1: C→T
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::OT), 1);
    }

    #[cfg(not(feature = "experimental-seqair"))]
    #[test]
    fn taps_aware_multiple_real_mismatches_all_counted() {
        // Two non-TAPS mismatches on OT — both must count
        let segment = test_segment(b"ACGT");
        let record = test_record(b"two-errors", 99, 100, b"GCGG"); // A→G and T→G
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::OT), 2);
    }

    const CLEAN: AlignmentShape = AlignmentShape { terminal_repeat: false, soft_clipped: false };
    const CLIPPED: AlignmentShape = AlignmentShape { terminal_repeat: false, soft_clipped: true };
    const REPEAT: AlignmentShape = AlignmentShape { terminal_repeat: true, soft_clipped: false };

    #[test]
    fn both_mates_spanning_an_indel_cast_one_vote() {
        let mut votes = FragmentVotes::default();
        votes.prepare(2);
        assert!(votes.saw_indel(0));
        // Second mate of the same fragment: `duplicate` resolves it to slot 0.
        assert!(!votes.saw_indel(0));
    }

    #[test]
    fn a_fragment_is_alternate_or_noisy_reference_but_never_both() {
        // The mate configuration that used to be subtracted twice: one mate
        // soft-clipped over the indel, the other spanning it. Whichever the
        // pileup reports first, the fragment leaves `ref_count` via
        // `total_indel_reads` and must not leave it again via `noisy_ref_count`.
        for indel_first in [true, false] {
            let mut votes = FragmentVotes::default();
            votes.prepare(2);
            if indel_first {
                assert!(votes.saw_indel(0));
                votes.saw_reference(0, CLIPPED);
            } else {
                votes.saw_reference(0, CLIPPED);
                assert!(votes.saw_indel(0));
            }
            assert_eq!(votes.totals().noisy_ref_count, 0);
        }
    }

    #[test]
    fn noisy_reference_mates_are_counted_once_per_fragment() {
        let mut votes = FragmentVotes::default();
        votes.prepare(4);
        // One fragment, both mates in a terminal repeat.
        votes.saw_reference(0, REPEAT);
        votes.saw_reference(0, REPEAT);
        // A second, clean fragment.
        votes.saw_reference(2, CLEAN);

        let totals = votes.totals();
        assert_eq!(totals.noisy_ref_count, 1);
    }

    #[test]
    fn soft_clips_count_per_fragment_including_indel_carriers() {
        // `soft_clip_rate` divides by `reads.len()`, which is per fragment, so
        // the numerator has to be too or the rate can exceed 1.
        let mut votes = FragmentVotes::default();
        votes.prepare(3);
        votes.saw_soft_clip(0);
        votes.saw_soft_clip(0);
        assert!(votes.saw_indel(0));
        votes.saw_soft_clip(2);

        assert_eq!(votes.totals().soft_clip_count, 2);
    }

    #[test]
    fn votes_do_not_leak_between_positions() {
        let mut votes = FragmentVotes::default();
        votes.prepare(1);
        votes.saw_reference(0, REPEAT);
        assert_eq!(votes.totals().noisy_ref_count, 1);

        votes.prepare(1);
        assert_eq!(votes.totals().noisy_ref_count, 0);
    }
}
