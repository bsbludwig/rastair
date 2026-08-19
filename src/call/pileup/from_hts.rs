use super::{
    INDEL_REF_WINDOW_DOWN, INDEL_REF_WINDOW_LEN, INDEL_REF_WINDOW_UP,
    indels::{IndelAllele, IndelObservation},
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

        let mut raw_reads = Vec::with_capacity(max_reads);

        // NOTE: The pileup might have already had some reads filtered out by
        // the pileup-level filter, so we don't need to worry about flag and
        // read-group filtering here. We do still apply read masking and quality
        // filtering, however.
        let alignments = pile
            .alignments()
            .filter_map(|alignment| {
                alignment_to_read(alignment, segment.as_ref(), params, orientation_cache)
            })
            .filter(|(_, seen_base)| params.read_masking.filter(seen_base))
            .filter(|(_, seen_base)| params.quality.filter(seen_base))
            .take(max_reads);

        let reads = match collector {
            NameCollector::Skip => {
                for (_, read) in alignments {
                    raw_reads.push(read);
                }
                SimpleReads(raw_reads.into())
            }
            NameCollector::Collect(buf) => {
                buf.prepare(max_reads);
                let mut to_remove = SmallVec::<usize, 16>::new();
                for (name, read) in alignments {
                    let this_idx = raw_reads.len();
                    raw_reads.push(read);
                    if let Some(other_idx) = buf.see(name, this_idx) {
                        resolve_pair(&raw_reads, this_idx, other_idx, &mut to_remove);
                    }
                }
                to_remove.sort_unstable();
                for &idx in to_remove.iter().rev() {
                    raw_reads.swap_remove(idx);
                }
                SimpleReads(raw_reads.into())
            }
        };

        let reference_base: Base =
            segment.sequence.get(idx).wrap_err("failed to get reference base")?.into();

        let context =
            SequenceContext::new(idx, &segment).wrap_err("failed to get sequence context")?;

        // Second pass: collect indel observations and compute depth_offset.
        // FIXME: Do in same pass
        let segment_start = segment.range.region.start as usize;
        let indel_cutoff = params.indel_end_of_read_cutoff;
        let mut indel_observations = SmallVec::new();
        let mut depth_offset: u32 = 0;
        let mut soft_clip_count: u32 = 0;

        // Indel accounting is per *fragment*, matching the deduplication already
        // applied to `reads` above. Both sides of `VAF = alt / depth` must share a
        // granularity: counting observations per alignment while the depth is per
        // fragment double-counts overlapping mates, and `ref_count` (computed as
        // `reads.len() - total_indel_reads`) can then floor to zero and read VAF as
        // 1.0. With `--keep-overlapping-reads` nothing is deduplicated and both
        // sides stay at the alignment level, which is equally consistent.
        let dedup_indels = !params.keep_overlapping_reads;
        let mut seen_indel_fragments: SmallVec<u64, 8> = SmallVec::new();
        let mut seen_depth_offset: SmallVec<u64, 8> = SmallVec::new();
        let mut noisy_ref_count: u32 = 0;

        for a in pile.alignments() {
            let record = a.record_view();
            let flags = record.flags();
            let (seq, qual) = record.seq_and_qual();
            let read_len = seq.len();
            let (matches, indels_in_read) = calc_cigar_data(record.raw_cigar());
            let fragment_key = dedup_indels.then(|| hash_bytes(record.qname()));

            let soft_clipped = has_soft_clip(record.raw_cigar());
            if soft_clipped {
                soft_clip_count += 1;
            }

            match a.indel() {
                Indel::None => {
                    let terminal_repeat = has_repeat_seq(&seq, 1, HOMOPOLYMER_UNITS)
                        || has_repeat_seq(&seq, 2, DINUCLEOTIDE_UNITS);
                    if (terminal_repeat || soft_clipped)
                        && first_seen(&mut seen_depth_offset, fragment_key)
                    {
                        if terminal_repeat {
                            depth_offset += 1;
                        }
                        noisy_ref_count += 1;
                    }
                }
                indel => {
                    let Some(qpos) = a.qpos() else { continue };

                    // End-of-read filter (stricter than SNVs)
                    if qpos < indel_cutoff || qpos >= read_len.saturating_sub(indel_cutoff) {
                        continue;
                    }

                    let strand = orientation_cache
                        .strand_for_alignment(&a, &segment, params)
                        .unwrap_or(Strand::Unknown);

                    let mismatches =
                        mismatch_cache.mismatch_count_for_alignment(&a, &segment, strand);
                    if mismatches > params.indel_max_mismatches {
                        trace!(
                            mismatches,
                            max = params.indel_max_mismatches,
                            "Indel skipped: too many non-TAPS mismatches"
                        );
                        continue;
                    }

                    let allele = match indel {
                        Indel::Ins(len) => {
                            let start = qpos + 1;
                            let end = start + len as usize;
                            let bases: SmallVec<Base, 4> =
                                (start..end).filter_map(|i| seq.get(i).map(Base::from)).collect();
                            if bases.is_empty() {
                                continue;
                            }
                            IndelAllele::Insertion(bases)
                        }
                        Indel::Del(len) => {
                            let ref_start = (pos as usize + 1).saturating_sub(segment_start);
                            let ref_end = ref_start + len as usize;
                            let bases: SmallVec<Base, 4> = segment
                                .sequence
                                .get(ref_start..ref_end)
                                .map(|slice| slice.iter().copied().map(Base::from).collect())
                                .unwrap_or_default();
                            if bases.is_empty() {
                                continue;
                            }
                            IndelAllele::Deletion(bases)
                        }
                        Indel::None => unreachable!(),
                    };

                    let base_qual = qual.get(qpos).copied().unwrap_or(0);
                    let mapq = record.mapq();

                    let insertion_base_quals = match &allele {
                        IndelAllele::Insertion(bases) => {
                            let start = qpos + 1;
                            let end = start + bases.len();
                            (start..end).filter_map(|i| qual.get(i).copied()).collect()
                        }
                        IndelAllele::Deletion(_) => SmallVec::new(),
                    };

                    // One vote per fragment: the second overlapping mate of a
                    // fragment that already contributed is skipped.
                    if !first_seen(&mut seen_indel_fragments, fragment_key) {
                        continue;
                    }

                    let has_repeat = has_repeat_seq(&seq, 1, HOMOPOLYMER_UNITS)
                        || has_repeat_seq(&seq, 2, DINUCLEOTIDE_UNITS);
                    let post_del_base_qual = match &allele {
                        IndelAllele::Deletion(_) => qual.get(qpos + 1).copied().unwrap_or(0),
                        IndelAllele::Insertion(_) => 0,
                    };

                    indel_observations.push(IndelObservation {
                        allele,
                        strand,
                        reverse: flags & 0x10 != 0,
                        pos_in_read: qpos as u32,
                        read_length: read_len as u32,
                        mapq,
                        base_qual,
                        matching_bases: matches,
                        num_indels_in_read: indels_in_read,
                        insertion_base_quals,
                        post_del_base_qual,
                        has_repeat,
                        noisy: has_repeat || soft_clipped,
                    });
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

fn alignment_to_read<'a>(
    a: Alignment<'a>,
    segment: &Segment,
    params: &PileupMappingParams,
    orientation_cache: &mut ReadOrientationCache,
) -> Option<(&'a [u8], SimpleRead)> {
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

/// Records a per-fragment contribution once. `None` (dedup disabled) always counts.
fn first_seen(seen: &mut SmallVec<u64, 8>, key: Option<u64>) -> bool {
    match key {
        None => true,
        Some(k) if seen.contains(&k) => false,
        Some(k) => {
            seen.push(k);
            true
        }
    }
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

    #[test]
    fn first_seen_dedups_by_key() {
        // Per-fragment dedup: a repeated key (second overlapping mate) is not counted.
        let mut seen: SmallVec<u64, 8> = SmallVec::new();
        assert!(first_seen(&mut seen, Some(7)));
        assert!(!first_seen(&mut seen, Some(7)));
        assert!(first_seen(&mut seen, Some(9)));
        assert!(!first_seen(&mut seen, Some(9)));
    }

    #[test]
    fn first_seen_none_is_always_true() {
        // dedup disabled (`--keep-overlapping-reads`): every alignment counts, nothing tracked.
        let mut seen: SmallVec<u64, 8> = SmallVec::new();
        assert!(first_seen(&mut seen, None));
        assert!(first_seen(&mut seen, None));
        assert!(seen.is_empty());
    }
}
