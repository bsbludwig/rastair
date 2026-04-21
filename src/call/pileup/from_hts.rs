use super::{
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
use rastair_types::{Base, SmallVec, Strand, strand_from_flags};
use rust_htslib::bam::{
    Record,
    ext::BamRecordExtensions as _,
    pileup::{Alignment, Indel, Pileup as HtsPileup},
};
use rustc_hash::{FxHashMap, FxHasher};
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
            return strand_from_flags(record.flags()).ok();
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
        let repeat_limit = params.indel_repeat_limit;
        let mut indel_observations = SmallVec::new();
        let mut depth_offset: u32 = 0;

        for a in pile.alignments() {
            let record = a.record_view();
            let flags = record.flags();
            let (seq, _qual) = record.seq_and_qual();
            let read_len = seq.len();

            match a.indel() {
                Indel::None => {
                    if has_soft_clip(record.raw_cigar())
                        || has_repeat_seq(&seq, 1, repeat_limit)
                        || has_repeat_seq(&seq, 2, repeat_limit)
                    {
                        depth_offset += 1;
                    }
                }
                indel => {
                    let Some(qpos) = a.qpos() else { continue };

                    // End-of-read filter (stricter than SNVs)
                    if qpos < indel_cutoff || qpos >= read_len.saturating_sub(indel_cutoff) {
                        continue;
                    }

                    let strand = strand_from_flags(flags);

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

                    indel_observations.push(IndelObservation {
                        allele,
                        strand,
                        reverse: flags & 0x10 != 0,
                        pos_in_read: qpos as u32,
                        read_length: read_len as u32,
                    });
                }
            }
        }

        Ok(Pileup {
            region: segment.range.clone(),
            context,
            pos: pile.pos(),
            reads,
            reference_base,
            indel_observations,
            depth_offset,
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

/// Check if first or last `cutoff` bases of a read form a repeating pattern of length `n`.
fn has_repeat_seq(seq: &rust_htslib::bam::record::Seq<'_>, n: usize, cutoff: usize) -> bool {
    let len = seq.len();
    if len < cutoff || n == 0 || cutoff < n {
        return false;
    }

    // Check start: do first `cutoff` bases repeat a pattern of length `n`?
    let start_pattern: SmallVec<u8, 4> = (0..n).filter_map(|i| seq.get(i)).collect();
    if start_pattern.len() == n {
        let start_repeat = (n..cutoff).all(|i| {
            seq.get(i).map_or(false, |b| start_pattern.get(i % n).map_or(false, |&p| b == p))
        });
        if start_repeat {
            return true;
        }
    }

    // Check end: do last `cutoff` bases repeat a pattern of length `n`?
    let end_start = len.saturating_sub(n);
    let end_pattern: SmallVec<u8, 4> = (end_start..len).filter_map(|i| seq.get(i)).collect();
    if end_pattern.len() == n {
        let check_start = len.saturating_sub(cutoff);
        (check_start..end_start).all(|i| {
            seq.get(i).map_or(false, |b| {
                let offset = (i - check_start) % n;
                // Align pattern index: the tail pattern repeats back from end
                end_pattern.get(offset % n).map_or(false, |&p| b == p)
            })
        })
    } else {
        false
    }
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
}
