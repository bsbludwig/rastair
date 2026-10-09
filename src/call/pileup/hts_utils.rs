#![cfg_attr(feature = "experimental-seqair", allow(unused, reason = "seqair integration WIP"))]

use crate::{call::process::PileupMappingParams, sequence::Segment};
use rust_htslib::bam::{Record, ext::BamRecordExtensions as _, pileup::Alignment};
use rustc_hash::{FxHashMap, FxHasher};
use seqair_types::{BamFlags, Base, Strand};
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
            let flags = BamFlags::from(record.flags());
            return Strand::from(flags).ok();
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

fn hash_bytes(bytes: &[u8]) -> u64 {
    let mut hasher = FxHasher::default();
    bytes.hash(&mut hasher);
    hasher.finish()
}
