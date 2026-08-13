use super::{
    INDEL_REF_WINDOW_DOWN, INDEL_REF_WINDOW_LEN, INDEL_REF_WINDOW_UP,
    indels::{IndelAllele, IndelObservation, TerminalRepeatLimits},
    overlapping_reads::{NameCollector, resolve_pair},
};
use crate::{
    call::{
        pileup::{Pileup, PositionInRead, SimpleRead, SimpleReads},
        process::PileupMappingParams,
    },
    sequence::Segment,
    utils::SequenceContext,
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use rust_htslib::bam::{
    Record,
    ext::BamRecordExtensions as _,
    pileup::{Alignment, Indel, Pileup as HtsPileup},
};
use rustc_hash::{FxHashMap, FxHasher};
use seqair_types::{Base, SmallVec, Strand, strand_from_flags};
use std::{
    hash::{Hash, Hasher},
    rc::Rc,
};
use tracing::{debug, instrument, trace};

#[derive(Default)]
pub(crate) struct ReadOrientationCache {
    strands: FxHashMap<ReadOrientationCacheKey, Strand>,
}

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

#[derive(Default)]
pub(crate) struct ReadMismatchCache {
    counts: FxHashMap<ReadOrientationCacheKey, u32>,
}

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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct ReadOrientationCacheKey {
    name_hash: u64,
    start: i64,
    flags: u16,
}

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

impl Pileup {
    /// Convert a pileup from htslib into our internal Pileup representation.
    #[instrument(level = "trace", skip_all)]
    pub(crate) fn from_hts(
        pile: &HtsPileup,
        segment: Rc<Segment>,
        params: &PileupMappingParams,
        collector: &mut NameCollector,
        orientation_cache: &mut ReadOrientationCache,
        mismatch_cache: &mut ReadMismatchCache,
    ) -> Result<Pileup> {
        let pos = pile.pos();
        let idx = segment.pos_to_idx(pos)?;
        let depth = pile.depth();
        let max_reads = depth.min(params.max_coverage);
        if depth > max_reads {
            debug!(pos, depth, "Capping number of reads in pileup to {max_reads}");
        }
        let max_reads = usize::try_from(max_reads).wrap_err("max_reads exceeds usize")?;

        let reference_base: Base =
            segment.sequence.get(idx).wrap_err("failed to get reference base")?.into();

        let segment_start = segment.range.region.start as usize;

        let indel_ctx = params.collect_indels.then_some(IndelContext {
            segment_start,
            pos,
            end_of_read_cutoff: params.indel_end_of_read_cutoff,
            repeat_limits: params.indel_repeat_limits,
            max_mismatches: params.indel_max_mismatches,
        });

        let mut raw_reads = Vec::with_capacity(max_reads);
        // Indel bookkeeping is built in lock-step with `raw_reads` so that the read
        // filters, the `max_coverage` cap and the overlap deduplication below apply
        // to both identically. Indel depth arithmetic (`ref_count = reads.len() -
        // total_indel_reads`) is only coherent if the two describe the same set of
        // fragments; computing them over separate passes over `pile.alignments()`
        // silently mixes granularities and biases VAF. Left empty when indels are off.
        let mut indel_data: Vec<IndelReadData> =
            if indel_ctx.is_some() { Vec::with_capacity(max_reads) } else { Vec::new() };

        // NOTE: The pileup might have already had some reads filtered out by
        // the pileup-level filter, so we don't need to worry about flag and
        // read-group filtering here. We do still apply read masking and quality
        // filtering, however.
        let alignments = pile
            .alignments()
            .filter_map(|alignment| {
                alignment_to_read(
                    alignment,
                    segment.as_ref(),
                    params,
                    orientation_cache,
                    mismatch_cache,
                    indel_ctx.as_ref(),
                )
            })
            .filter(|(_, seen_base, _)| params.read_masking.filter(seen_base))
            .filter(|(_, seen_base, _)| params.quality.filter(seen_base))
            .take(max_reads);

        let reads = match collector {
            NameCollector::Skip => {
                for (_, read, indel) in alignments {
                    raw_reads.push(read);
                    if params.collect_indels {
                        indel_data.push(indel);
                    }
                }
                SimpleReads(raw_reads.into())
            }
            NameCollector::Collect(buf) => {
                buf.prepare(max_reads);
                let mut to_remove = SmallVec::<usize, 16>::new();
                for (name, read, indel) in alignments {
                    let this_idx = raw_reads.len();
                    raw_reads.push(read);
                    if params.collect_indels {
                        indel_data.push(indel);
                    }
                    if let Some(other_idx) = buf.see(name, this_idx) {
                        resolve_pair(&raw_reads, this_idx, other_idx, &mut to_remove);
                    }
                }
                to_remove.sort_unstable();
                for &idx in to_remove.iter().rev() {
                    raw_reads.swap_remove(idx);
                    if params.collect_indels {
                        indel_data.swap_remove(idx);
                    }
                }
                SimpleReads(raw_reads.into())
            }
        };

        let context =
            SequenceContext::new(idx, &segment).wrap_err("failed to get sequence context")?;

        // Fold the per-alignment indel bookkeeping into position-level counts. The
        // vector has already had the overlap dedup applied, so every entry is one
        // fragment and each fragment is either indel-supporting or reference-supporting
        // — never both.
        let mut indel_observations = SmallVec::new();
        let mut depth_offset: u32 = 0;
        let mut soft_clip_count: u32 = 0;
        let mut noisy_ref_count: u32 = 0;
        for data in indel_data {
            if data.soft_clipped {
                soft_clip_count += 1;
            }
            match data.observation {
                // Supporting fragments carry their own noise flag on the
                // observation, so `aggregate_indels` can subtract it from the
                // matching allele's count.
                Some(observation) => indel_observations.push(*observation),
                None => {
                    if data.terminal_repeat {
                        depth_offset += 1;
                    }
                    if data.soft_clipped || data.terminal_repeat {
                        noisy_ref_count += 1;
                    }
                }
            }
        }

        // The reference window is only needed for indel slippage features, so
        // skip the allocation at the (vast majority of) positions without indels.
        let (indel_ref_window, indel_ref_anchor) = if indel_observations.is_empty() {
            (SmallVec::new(), 0)
        } else {
            indel_ref_window_at(idx, &segment)
        };

        Ok(Pileup {
            region: segment.range.clone(),
            context,
            pos: pile.pos(),
            reads,
            reference_base,
            indel_observations,
            depth_offset,
            homopolymer_run: homopolymer_run_at(pos as usize, &segment, segment_start),
            dinucleotide_run: dinucleotide_run_at(pos as usize, &segment, segment_start),
            soft_clip_count,
            noisy_ref_count,
            indel_ref_window,
            indel_ref_anchor,
        })
    }
}

fn alignment_to_read<'a>(
    a: Alignment<'a>,
    segment: &Segment,
    params: &PileupMappingParams,
    orientation_cache: &mut ReadOrientationCache,
    mismatch_cache: &mut ReadMismatchCache,
    indel_ctx: Option<&IndelContext>,
) -> Option<(&'a [u8], SimpleRead, IndelReadData)> {
    let pos = a.qpos()?;
    let record = a.record_view();
    let flags = record.flags();
    let strand = orientation_cache.strand_for_alignment(&a, segment, params)?;
    let (matches, indels) = calc_cigar_data(record.raw_cigar());
    let (seq, qual) = record.seq_and_qual();

    let before_base = pos.checked_sub(1).map(|i| Base::from(seq[i]));
    let after_base = match a.indel() {
        Indel::None if pos + 1 < seq.len() => Some(Base::from(seq[pos + 1])),
        _ => None,
    };

    let observed_base: Base = seq[pos].into();
    let indel_data = match indel_ctx {
        Some(ctx) => ctx.read_data(
            &a,
            segment,
            mismatch_cache,
            IndelReadInput {
                strand,
                qpos: pos,
                seq: &seq,
                qual,
                matches,
                indels_in_read: indels,
                mapq: record.mapq(),
                reverse: flags & 0x10 != 0,
                soft_clipped: has_soft_clip(record.raw_cigar()),
            },
        ),
        None => IndelReadData::default(),
    };

    Some((
        record.qname(),
        SimpleRead {
            base: observed_base,
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
        indel_data,
    ))
}

/// Position-level inputs the indel bookkeeping needs, hoisted out of the
/// per-alignment loop.
struct IndelContext {
    segment_start: usize,
    pos: u32,
    end_of_read_cutoff: usize,
    repeat_limits: TerminalRepeatLimits,
    max_mismatches: u32,
}

/// Per-alignment values already computed by [`alignment_to_read`].
struct IndelReadInput<'a> {
    strand: Strand,
    qpos: usize,
    seq: &'a rust_htslib::bam::record::Seq<'a>,
    qual: &'a [u8],
    matches: u32,
    indels_in_read: u32,
    mapq: u8,
    reverse: bool,
    soft_clipped: bool,
}

/// Indel bookkeeping for one alignment, carried alongside its [`SimpleRead`] so
/// that overlap deduplication removes both together.
///
/// An entry either carries an [`IndelObservation`] or contributes to the depth
/// penalties — the two are mutually exclusive, which is what makes
/// `ref_count = reads.len() - total_indel_reads` a real partition.
#[derive(Default)]
struct IndelReadData {
    /// Boxed: present on a small minority of reads, and this struct is stored
    /// once per read in a pileup-sized vector.
    observation: Option<Box<IndelObservation>>,
    soft_clipped: bool,
    /// Terminal homopolymer or dinucleotide repeat; drives the ML-facing
    /// `depth_offset`.
    terminal_repeat: bool,
}

impl IndelContext {
    fn read_data(
        &self,
        a: &Alignment<'_>,
        segment: &Segment,
        mismatch_cache: &mut ReadMismatchCache,
        read: IndelReadInput<'_>,
    ) -> IndelReadData {
        let terminal_repeat = has_repeat_seq(read.seq, 1, self.repeat_limits.homopolymer_units)
            || has_repeat_seq(read.seq, 2, self.repeat_limits.dinucleotide_units);
        let noisy = read.soft_clipped || terminal_repeat;
        let mut data =
            IndelReadData { observation: None, soft_clipped: read.soft_clipped, terminal_repeat };

        let Some(observation) =
            self.observation(a, segment, mismatch_cache, &read, terminal_repeat, noisy)
        else {
            return data;
        };

        data.observation = Some(Box::new(observation));
        data
    }

    /// The indel this alignment supports at the anchor, if any survives the
    /// read-level indel filters.
    fn observation(
        &self,
        a: &Alignment<'_>,
        segment: &Segment,
        mismatch_cache: &mut ReadMismatchCache,
        read: &IndelReadInput<'_>,
        has_repeat: bool,
        noisy: bool,
    ) -> Option<IndelObservation> {
        let indel = match a.indel() {
            Indel::None => return None,
            indel => indel,
        };

        let read_len = read.seq.len();
        // End-of-read filter (stricter than SNVs)
        if read.qpos < self.end_of_read_cutoff
            || read.qpos >= read_len.saturating_sub(self.end_of_read_cutoff)
        {
            return None;
        }

        let mismatches = mismatch_cache.mismatch_count_for_alignment(a, segment, read.strand);
        if mismatches > self.max_mismatches {
            trace!(
                mismatches,
                max = self.max_mismatches,
                "Indel skipped: too many non-TAPS mismatches"
            );
            return None;
        }

        let allele = match indel {
            Indel::Ins(len) => {
                let start = read.qpos + 1;
                let end = start + len as usize;
                let bases: SmallVec<Base, 4> =
                    (start..end).filter_map(|i| read.seq.get(i).map(Base::from)).collect();
                if bases.is_empty() {
                    return None;
                }
                IndelAllele::Insertion(bases)
            }
            Indel::Del(len) => {
                let ref_start = (self.pos as usize + 1).saturating_sub(self.segment_start);
                let ref_end = ref_start + len as usize;
                let bases: SmallVec<Base, 4> = segment
                    .sequence
                    .get(ref_start..ref_end)
                    .map(|slice| slice.iter().copied().map(Base::from).collect())
                    .unwrap_or_default();
                if bases.is_empty() {
                    return None;
                }
                IndelAllele::Deletion(bases)
            }
            Indel::None => return None,
        };

        let insertion_base_quals = match &allele {
            IndelAllele::Insertion(bases) => {
                let start = read.qpos + 1;
                let end = start + bases.len();
                (start..end).filter_map(|i| read.qual.get(i).copied()).collect()
            }
            IndelAllele::Deletion(_) => SmallVec::new(),
        };
        let post_del_base_qual = match &allele {
            IndelAllele::Deletion(_) => read.qual.get(read.qpos + 1).copied().unwrap_or(0),
            IndelAllele::Insertion(_) => 0,
        };

        Some(IndelObservation {
            allele,
            strand: read.strand,
            reverse: read.reverse,
            pos_in_read: read.qpos as u32,
            read_length: read_len as u32,
            mapq: read.mapq,
            base_qual: read.qual.get(read.qpos).copied().unwrap_or(0),
            matching_bases: read.matches,
            num_indels_in_read: read.indels_in_read,
            insertion_base_quals,
            post_del_base_qual,
            has_repeat,
            noisy,
        })
    }
}

/// Whether an observed/reference base pair is expected TAPS conversion signal
/// rather than a sequencing error or a variant.
///
/// For [`Strand::Unknown`] both patterns are accepted, so that reads whose
/// orientation could not be determined are never penalised for methylation.
fn is_taps_signal(observed: Base, reference: Base, strand: Strand) -> bool {
    match strand {
        Strand::OT => observed == Base::T && reference == Base::C,
        Strand::OB => observed == Base::A && reference == Base::G,
        Strand::Unknown => {
            (observed == Base::T && reference == Base::C)
                || (observed == Base::A && reference == Base::G)
        }
    }
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

        // Mismatches that are expected TAPS signal are not sequencing errors.
        if !is_taps_signal(observed, reference, strand) {
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

fn hash_bytes(bytes: &[u8]) -> u64 {
    let mut hasher = FxHasher::default();
    bytes.hash(&mut hasher);
    hasher.finish()
}

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

/// Check if a CIGAR array contains a soft-clip operation (op 4).
fn has_soft_clip(cigar: &[u32]) -> bool {
    cigar.iter().any(|&c| c & 0xF == 4)
}

/// Reference bases around the anchor at segment index `idx`, plus the anchor's
/// index within the returned window. Clamped at segment boundaries.
fn indel_ref_window_at(
    idx: usize,
    segment: &Segment,
) -> (SmallVec<Base, INDEL_REF_WINDOW_LEN>, u8) {
    let seq = &segment.sequence;
    let start = idx.saturating_sub(INDEL_REF_WINDOW_UP);
    let end = (idx + INDEL_REF_WINDOW_DOWN + 1).min(seq.len());
    let window = seq.get(start..end).unwrap_or(&[]).iter().map(|&b| Base::from(b)).collect();
    let anchor = u8::try_from(idx - start).unwrap_or(0);
    (window, anchor)
}

fn homopolymer_run_at(pos: usize, segment: &Segment, segment_start: usize) -> u8 {
    let seq = &segment.sequence;
    let idx = pos.saturating_sub(segment_start);
    let Some(&center) = seq.get(idx) else { return 0 };
    let mut run = 1u8;
    let mut i = idx;
    while i > 0 {
        i -= 1;
        if seq.get(i) == Some(&center) {
            run += 1;
        } else {
            break;
        }
    }
    i = idx;
    while i + 1 < seq.len() {
        i += 1;
        if seq.get(i) == Some(&center) {
            run += 1;
        } else {
            break;
        }
    }
    run
}

fn dinucleotide_run_at(pos: usize, segment: &Segment, segment_start: usize) -> u8 {
    let seq = &segment.sequence;
    let idx = pos.saturating_sub(segment_start);
    let try_phase = |start: usize| -> u8 {
        if start + 1 >= seq.len() {
            return 0;
        }
        let p0 = seq[start];
        let p1 = seq[start + 1];
        if p0 == p1 {
            return 0;
        }
        let mut run = 2u8;
        let mut i = start;
        while i >= 2 {
            if seq.get(i - 2) == Some(&p0) && seq.get(i - 1) == Some(&p1) {
                run += 2;
                i -= 2;
            } else {
                break;
            }
        }
        i = start + 2;
        while i + 1 < seq.len() {
            if seq.get(i) == Some(&p0) && seq.get(i + 1) == Some(&p1) {
                run += 2;
                i += 2;
            } else {
                break;
            }
        }
        run
    };
    try_phase(idx).max(try_phase(idx.saturating_sub(1)))
}

/// Check if the first or last `cutoff` bases of a read form a repeating pattern of
/// period `n`.
///
/// A window is period-`n` exactly when every base equals the one `n` positions
/// along, which is what both branches test. Phrasing the tail check against a
/// captured trailing pattern instead is easy to get wrong: the pattern has to be
/// indexed from the *end* of the window, so anchoring it at the window start only
/// agrees when `n` divides `cutoff` — with the default `cutoff = 3`, `n = 2` then
/// silently degrades into a 2 bp homopolymer test.
/// Whether either terminus of `seq` is a tandem repeat of period `n` spanning at
/// least `units` whole repeat units.
///
/// `units` counts repeat units, not bases, so one threshold means the same thing
/// for every period. Fewer than two units is rejected outright: the window would
/// then be shorter than a single comparison and `all()` over an empty range is
/// vacuously true, i.e. every read would look like a repeat.
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
    use crate::{
        call::{process::PileupMappingParams, variant_calling::VariantCallingParams},
        sequence::{ChunkRegion, Region, Segment},
        utils::default,
    };
    use rust_htslib::bam::record::{Cigar, CigarString};

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

    fn repeats(sequence: &[u8], n: usize, units: usize) -> bool {
        let record = test_record(b"repeat", 0, 100, sequence);
        has_repeat_seq(&record.seq(), n, units)
    }

    #[test]
    fn homopolymer_repeats_detected_at_either_terminus() {
        assert!(repeats(b"AAAACGTCGT", 1, 4), "leading AAAA");
        assert!(repeats(b"ACGTCGTTTTT", 1, 4), "trailing TTTT");
        assert!(!repeats(b"ACGTCGTACG", 1, 4));
        assert!(!repeats(b"AAACGTCGTA", 1, 4), "three is one unit short of the threshold");
    }

    /// The tail window has to be indexed from the end, not the start: anchoring at
    /// the window start compares adjacent bases instead of bases 2 apart, which
    /// silently turns the dinucleotide check into a 2 bp homopolymer check and
    /// misses real trailing dinucleotide repeats.
    #[test]
    fn dinucleotide_repeats_detected_at_either_terminus() {
        assert!(repeats(b"ACGTATATAT", 2, 3), "trailing ATATAT is three period-2 units");
        assert!(repeats(b"ATATATCCGG", 2, 3), "leading ATATAT is three period-2 units");
        assert!(!repeats(b"ACGTATATAG", 2, 3));
        assert!(!repeats(b"ACGTCGATAT", 2, 3), "two units is one short of the threshold");
        // Not a homopolymer at either end, so period 1 must not fire.
        assert!(!repeats(b"ACGTATATAT", 1, 4));
    }

    #[test]
    fn repeat_check_needs_a_window() {
        assert!(!repeats(b"AA", 1, 4), "sequence shorter than the window");
        assert!(!repeats(b"AAAAAA", 0, 4), "period 0 is not a repeat");
    }

    /// A window of fewer than two units has nothing to compare, so the periodicity
    /// check runs over an empty range and `all()` is vacuously true. Guarding on
    /// `units` rather than on the derived window length is what keeps that from
    /// flagging every read: the old shared-base-window form let `n == cutoff`
    /// through and returned `true` unconditionally.
    #[test]
    fn a_window_below_two_units_is_never_a_repeat() {
        for n in 1..=4 {
            for units in 0..=1 {
                assert!(!repeats(b"ACGTACGTAC", n, units), "period {n}, {units} unit(s)");
                assert!(!repeats(b"AAAAAAAAAA", n, units), "period {n}, {units} unit(s)");
            }
        }
    }

    /// The flag is meant to mark a read as unusual. A shared 3 bp window made the
    /// period-2 arm reduce to `seq[0] == seq[2] || seq[len - 3] == seq[len - 1]`,
    /// which fires on 43.75% of random reads; the unit-based thresholds have to
    /// stay far below that or whatever consumes the flag is mostly noise.
    #[test]
    fn default_limits_flag_only_a_small_share_of_random_reads() {
        let limits = TerminalRepeatLimits::default();
        // Deterministic xorshift over 2-bit base codes — no dev-dependency needed.
        let mut state: u64 = 0x2545_F491_4F6C_DD1D;
        let mut flagged = 0usize;
        const TRIALS: usize = 20_000;
        for _ in 0..TRIALS {
            let seq: Vec<u8> = (0..60)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    b"ACGT"[(state >> 33) as usize % 4]
                })
                .collect();
            if repeats(&seq, 1, limits.homopolymer_units)
                || repeats(&seq, 2, limits.dinucleotide_units)
            {
                flagged += 1;
            }
        }
        let rate = flagged as f64 / TRIALS as f64;
        assert!(rate < 0.10, "default limits flagged {:.1}% of random reads", rate * 100.0);
    }

    #[test]
    fn name_collector_skip_when_keep_overlapping() {
        let params = PileupMappingParams {
            variant_calling: VariantCallingParams { keep_overlapping_reads: true, ..default() },
            ..default()
        };
        assert!(matches!(NameCollector::new(&params), NameCollector::Skip));
    }

    #[test]
    fn name_collector_collect_by_default() {
        let params = PileupMappingParams::default();
        assert!(matches!(NameCollector::new(&params), NameCollector::Collect(_)));
    }

    /// Regression guard for the indel overlap double-counting fix: an indel supported
    /// by overlapping mate-pairs must contribute one vote per *fragment* (the depth
    /// denominator's granularity), not one per alignment. Drives the real
    /// `get_pileups` → `Pileup::from_hts` path over a synthetic BAM.
    #[test]
    fn indel_observations_are_deduped_per_fragment() -> color_eyre::Result<()> {
        use crate::{
            call::process::get_pileups,
            sequence::ReaderParams,
            utils::{CliRegionInput, RegionString},
        };
        use rust_htslib::bam::{self, Header, Record, header::HeaderRecord};
        use seqair_types::Pos1;

        // Non-repetitive reference (no homopolymer/dinucleotide runs) so the ref-noise
        // and depth offsets stay zero and we isolate the alt-count dedup.
        let refseq = "ACGT".repeat(15);
        let refbytes = refseq.as_bytes();
        let dir = tempfile::TempDir::new()?;
        let fasta = dir.path().join("ref.fasta");
        std::fs::write(&fasta, format!(">chrT\n{refseq}\n"))?;
        std::fs::write(
            dir.path().join("ref.fasta.fai"),
            format!("chrT\t{0}\t6\t{0}\t{1}\n", refseq.len(), refseq.len() + 1),
        )?;

        // A 1 bp insertion at the CIGAR anchor (ref index 29), supported by six
        // fragments, each written as two fully-overlapping mates with the same qname.
        // Three fragments are OT (flags 99/147) and three OB (flags 83/163). CIGAR
        // 10M1I10M from 0-based pos 20; the M bases match the reference.
        let start = 20usize;
        let seq: Vec<u8> = refbytes
            .iter()
            .skip(start)
            .take(10)
            .chain(b"A")
            .chain(refbytes.iter().skip(start + 10).take(10))
            .copied()
            .collect();
        let cigar = CigarString(vec![Cigar::Match(10), Cigar::Ins(1), Cigar::Match(10)].into());
        let quals = vec![40u8; seq.len()];
        let start_pos = i64::try_from(start)?;

        let bam_path = dir.path().join("reads.bam");
        {
            let mut header = Header::new();
            header.push_record(
                HeaderRecord::new(b"SQ").push_tag(b"SN", "chrT").push_tag(b"LN", refseq.len()),
            );
            let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam)?;
            // 99/147 is an F1R2 pair (OT); 83/163 is R1F2 (OB). Both mates of a
            // fragment share the OT/OB assignment but have opposite reverse flags.
            let fragments: [(&str, [u16; 2]); 2] = [("ot", [99, 147]), ("ob", [83, 163])];
            for (label, mate_flags) in fragments {
                for i in 0..3 {
                    let qname = format!("{label}frag{i}");
                    for &flags in &mate_flags {
                        let mut rec = Record::new();
                        rec.set(qname.as_bytes(), Some(&cigar), &seq, &quals);
                        rec.set_tid(0);
                        rec.set_pos(start_pos);
                        rec.set_mtid(0);
                        rec.set_mpos(start_pos);
                        rec.set_mapq(60);
                        rec.set_flags(flags);
                        writer.write(&rec)?;
                    }
                }
            }
        }
        bam::index::build(&bam_path, None, bam::index::Type::Bai, 1)?;

        // (total observations, OT observations, OB observations)
        let observation_strands =
            |keep_overlapping: bool| -> color_eyre::Result<(usize, usize, usize)> {
                let mut params = ReaderParams::test_with(
                    bam_path.to_str().expect("utf8 bam path"),
                    fasta.to_str().expect("utf8 fasta path"),
                );
                params.regions = Some(CliRegionInput::from_region(RegionString {
                    chromosome: "chrT".into(),
                    start: Some(Pos1::new(1).expect("valid pos")),
                    end: Some(Pos1::new(60).expect("valid pos")),
                }));
                let mut readers = params.readers()?;
                let chunk = readers
                    .segments(1000, 0)?
                    .next()
                    .ok_or_else(|| color_eyre::eyre::eyre!("no segment fetched"))?;
                let pileup_params = PileupMappingParams {
                    variant_calling: VariantCallingParams {
                        keep_overlapping_reads: keep_overlapping,
                        ..default()
                    },
                    collect_indels: true,
                    ..default()
                };
                let (_segment, pileups) = get_pileups(&mut readers, &chunk, &pileup_params)?;
                let observations: Vec<Strand> = pileups
                    .flat_map(|p| p.indel_observations.iter().map(|o| o.strand).collect::<Vec<_>>())
                    .collect();
                let ot = observations.iter().filter(|s| **s == Strand::OT).count();
                let ob = observations.iter().filter(|s| **s == Strand::OB).count();
                Ok((observations.len(), ot, ob))
            };

        // Default: overlapping mates collapse to one vote per fragment (6 fragments).
        // Both strands must survive: OT/OB is a property of the fragment, so dropping
        // the second mate must not make the allele look single-stranded — that would
        // make the hard-filter strand-bias test see spurious skew for every indel
        // inside a mate overlap.
        assert_eq!(
            observation_strands(false)?,
            (6, 3, 3),
            "expected fragment-level indel counting with both strands represented"
        );
        // `--keep-overlapping-reads`: every alignment votes (2 mates x 6 fragments).
        assert_eq!(
            observation_strands(true)?,
            (12, 6, 6),
            "expected alignment-level counting with overlaps kept"
        );

        Ok(())
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

    #[test]
    fn taps_aware_ot_ct_not_counted() {
        // C→T on OT is TAPS methylation signal — must not count
        let segment = test_segment(b"ACGT");
        let record = test_record(b"ot-ct", 99, 100, b"ATGT"); // pos 1: C→T
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::OT), 0);
    }

    #[test]
    fn taps_aware_ob_ga_not_counted() {
        // G→A on OB is TAPS methylation signal — must not count
        let segment = test_segment(b"ACGT");
        let record = test_record(b"ob-ga", 83, 100, b"ACAT"); // pos 2: G→A
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::OB), 0);
    }

    #[test]
    fn taps_aware_ot_other_mismatch_counted() {
        // A→G on OT is a real sequencing error — must count
        let segment = test_segment(b"ACGT");
        let record = test_record(b"ot-ag", 99, 100, b"GCGT"); // pos 0: A→G
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::OT), 1);
    }

    #[test]
    fn taps_aware_ob_ct_counted() {
        // C→T on OB is not a TAPS signal (wrong strand) — must count
        let segment = test_segment(b"ACGT");
        let record = test_record(b"ob-ct", 83, 100, b"ATGT"); // pos 1: C→T on OB
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::OB), 1);
    }

    #[test]
    fn taps_aware_unknown_strand_excludes_both_patterns() {
        // Unknown strand: conservatively skip both C→T and G→A to avoid
        // penalising methylated reads whose strand couldn't be determined
        let segment = test_segment(b"CGCG");
        let record = test_record(b"unknown", 99, 100, b"TATA"); // C→T and G→A alternating
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::Unknown), 0);
    }

    #[test]
    fn taps_aware_no_mismatches_returns_zero() {
        let segment = test_segment(b"ACGT");
        let record = test_record(b"perfect", 99, 100, b"ACGT");
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::OT), 0);
    }

    #[test]
    fn taps_aware_mixed_counts_only_non_taps() {
        // A→G (real error, counted) + C→T (TAPS on OT, excluded) = 1
        let segment = test_segment(b"ACGT");
        let record = test_record(b"mixed", 99, 100, b"GTGT"); // pos 0: A→G, pos 1: C→T
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::OT), 1);
    }

    #[test]
    fn taps_aware_multiple_real_mismatches_all_counted() {
        // Two non-TAPS mismatches on OT — both must count
        let segment = test_segment(b"ACGT");
        let record = test_record(b"two-errors", 99, 100, b"GCGG"); // A→G and T→G
        assert_eq!(count_taps_aware_mismatches(&record, &segment, Strand::OT), 2);
    }
}
