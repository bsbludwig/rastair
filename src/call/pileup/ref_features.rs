use super::{INDEL_REF_WINDOW_DOWN, INDEL_REF_WINDOW_LEN, INDEL_REF_WINDOW_UP};
use crate::sequence::Segment;
use seqair_types::{Base, SmallVec};

/// Reference bases around the anchor at segment index `idx`, plus the anchor's
/// index within the returned window. Clamped at segment boundaries.
pub(crate) fn indel_ref_window_at(
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

pub(crate) fn homopolymer_run_at(pos: usize, segment: &Segment, segment_start: usize) -> u8 {
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

pub(crate) fn dinucleotide_run_at(pos: usize, segment: &Segment, segment_start: usize) -> u8 {
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
