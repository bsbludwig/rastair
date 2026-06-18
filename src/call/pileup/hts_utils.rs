use crate::{call::process::PileupMappingParams, sequence::Segment};
use rust_htslib::bam::{Record, ext::BamRecordExtensions as _, pileup::Alignment};
use rustc_hash::{FxHashMap, FxHasher};
use seqair_types::{Base, Strand, strand_from_flags};
use std::hash::{Hash as _, Hasher as _};
use tracing::trace;

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

fn pseudo_random_strand(record: &Record) -> Strand {
    let mut hasher = FxHasher::default();
    record.qname().hash(&mut hasher);
    record.pos().hash(&mut hasher);
    record.flags().hash(&mut hasher);
    if hasher.finish() & 1 == 0 { Strand::OT } else { Strand::OB }
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

/// Calculate matches and indels from a packed CIGAR array.
///
/// Lower 4 bits encode the operation; upper 28 bits encode the length.
pub fn calc_cigar_data(cigar: &[u32]) -> (u32, u32) {
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

#[derive(Default)]
pub(crate) struct ReadOrientationCache {
    strands: FxHashMap<ReadOrientationCacheKey, Strand>,
}

impl ReadOrientationCache {
    pub fn strand_for_alignment(
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
pub struct ReadMismatchCache {
    counts: FxHashMap<ReadOrientationCacheKey, u32>,
}

impl ReadMismatchCache {
    pub fn mismatch_count_for_alignment(
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

fn hash_bytes(bytes: &[u8]) -> u64 {
    let mut hasher = FxHasher::default();
    bytes.hash(&mut hasher);
    hasher.finish()
}

#[cfg(test)]
mod tests {
    use rust_htslib::bam::record::{Cigar, CigarString};

    use super::*;
    use crate::sequence::{ChunkRegion, Region};

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
