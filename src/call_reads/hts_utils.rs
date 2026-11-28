use crate::sequence::Segment;
use rustc_hash::FxHasher;
use seqair::bam::{RawBamRecord, aligned_pairs::MatchPosition};
use seqair_types::{Base, Pos0, Strand};
use std::hash::{Hash as _, Hasher as _};
use tracing::trace;

pub(crate) fn infer_strand_from_mismatch_motifs(
    record: &RawBamRecord,
    segment: &Segment,
) -> Strand {
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

fn collect_mismatch_motif_evidence(
    record: &RawBamRecord,
    segment: &Segment,
) -> ReadOrientationEvidence {
    let mut seq = Vec::with_capacity(record.seq_len());
    record.decode_seq_into(&mut seq);
    let mut evidence = ReadOrientationEvidence::default();

    for MatchPosition { qpos, rpos, .. } in crate::bam::aligned_bases(record) {
        let read_idx = qpos.as_usize();
        let Some(observed) = seq.get(read_idx).copied().map(Base::from) else {
            continue;
        };
        let Some(reference) = reference_base(segment, rpos) else {
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

fn pseudo_random_strand(record: &RawBamRecord) -> Strand {
    let mut hasher = FxHasher::default();
    record.qname().hash(&mut hasher);
    // As 64-bit position and 16-bit flags: the hash decides ties, so its
    // input keeps the shape it always had.
    record.pos().map_or(-1, |pos| i64::from(pos.as_u32())).hash(&mut hasher);
    record.flags().raw().hash(&mut hasher);
    if hasher.finish() & 1 == 0 { Strand::OT } else { Strand::OB }
}

fn reference_base(segment: &Segment, pos_in_ref: Pos0) -> Option<Base> {
    let idx = segment.pos_to_idx(pos_in_ref).ok()?;
    segment.sequence.get(idx).copied().map(Base::from)
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

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::sequence::{ChunkRegion, Region};
    use seqair::reader::AlignmentStream;

    pub(crate) fn test_segment(sequence: &[u8]) -> Segment {
        let start = 100u64;
        let end = start + u64::try_from(sequence.len()).expect("sequence length fits") - 1;
        Segment {
            range: std::sync::Arc::new(ChunkRegion {
                region: Region { contig: "chrTest".into(), start, end },
                last_position: end,
                overlap_start: 0,
                overlap_end: 0,
            }),
            sequence: sequence.to_vec(),
        }
    }

    /// A read of `sequence` aligned without gaps at the 0-based `start`, its
    /// mate at the contig's first base; read through the stream reader, since
    /// a record has no public constructor.
    pub(crate) fn test_record(
        qname: &[u8],
        flags: u16,
        start: i64,
        sequence: &[u8],
    ) -> RawBamRecord {
        let sam = format!(
            "@HD\tVN:1.6\n@SQ\tSN:chrTest\tLN:100000\n{}\t{flags}\tchrTest\t{}\t60\t{}M\t=\t1\t0\t{}\t{}\n",
            String::from_utf8_lossy(qname),
            start + 1,
            sequence.len(),
            String::from_utf8_lossy(sequence),
            "I".repeat(sequence.len()),
        );
        let mut stream = AlignmentStream::from_reader(std::io::Cursor::new(sam.into_bytes()))
            .expect("the SAM text opens");
        let mut record = RawBamRecord::default();
        assert!(stream.read_record(&mut record).expect("the SAM text parses"), "no record");
        record
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
