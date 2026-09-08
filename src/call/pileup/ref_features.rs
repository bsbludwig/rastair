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

/// The reference tract lengths an indel anchored at `anchor` sits in.
///
/// `anchor` is the pileup column: the base *before* the indel. A left-aligned
/// indel therefore starts one base later, and that next base is where the tract
/// has to be measured — at the anchor these read ~1 exactly where the tract is
/// longest, which inverts the feature.
///
/// The `+ 1` and the segment-relative index live here, and the run functions
/// are private, so the two backends cannot disagree about the convention. They
/// did: `from_seqair` measured at the anchor while `from_hts` measured at the
/// tract, which cost ~1 point of indel F1 on chr12.
pub(crate) fn indel_tract_runs_at(anchor: u64, segment: &Segment) -> TractRuns {
    let idx = (anchor + 1).saturating_sub(segment.range.region.start) as usize;
    TractRuns {
        homopolymer: homopolymer_run_at(idx, segment),
        dinucleotide: dinucleotide_run_at(idx, segment),
    }
}

/// Lengths of the reference repeat tracts an indel sits in, in bases.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct TractRuns {
    pub(crate) homopolymer: u8,
    pub(crate) dinucleotide: u8,
}

fn homopolymer_run_at(idx: usize, segment: &Segment) -> u8 {
    let seq = &segment.sequence;
    let Some(&center) = seq.get(idx) else { return 0 };
    let mut run = 1u8;
    let mut i = idx;
    while i > 0 {
        i -= 1;
        if seq.get(i) == Some(&center) {
            run = run.saturating_add(1);
        } else {
            break;
        }
    }
    i = idx;
    while i + 1 < seq.len() {
        i += 1;
        if seq.get(i) == Some(&center) {
            run = run.saturating_add(1);
        } else {
            break;
        }
    }
    run
}

fn dinucleotide_run_at(idx: usize, segment: &Segment) -> u8 {
    let seq = &segment.sequence;
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
                run = run.saturating_add(2);
                i -= 2;
            } else {
                break;
            }
        }
        i = start + 2;
        while i + 1 < seq.len() {
            if seq.get(i) == Some(&p0) && seq.get(i + 1) == Some(&p1) {
                run = run.saturating_add(2);
                i += 2;
            } else {
                break;
            }
        }
        run
    };
    try_phase(idx).max(try_phase(idx.saturating_sub(1)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::{ChunkRegion, Region};

    /// A run counter is a `u8` but a reference tract can be longer than 255. The
    /// wrapped value would land back in the short-tract rows of
    /// `ALLELE_FRACTION_BY_TRACT`, applying the simple-sequence prior inside the
    /// longest repeat in the genome — and panic outright in a debug build.
    #[test]
    fn a_reference_run_longer_than_the_counter_saturates() {
        let homopolymer = test_segment(&[b'A'; 400]);
        assert_eq!(indel_tract_runs_at(299, &homopolymer).homopolymer, u8::MAX);

        let dinucleotide = test_segment(&b"AT".repeat(200));
        assert!(indel_tract_runs_at(299, &dinucleotide).dinucleotide >= u8::MAX - 1);
    }

    /// The runs are measured one base past the anchor, because that is where a
    /// left-aligned indel starts. Anchored on the `C` before an 8-base `A`
    /// tract, the homopolymer feature must read the tract, not the anchor's own
    /// run of one.
    #[test]
    fn runs_are_measured_at_the_tract_not_the_anchor() {
        //                        idx: 0123456789...
        let segment = test_segment(b"GGGCAAAAAAAATTT");
        // Anchor 103 is segment index 3, the `C`.
        let runs = indel_tract_runs_at(103, &segment);
        assert_eq!(runs.homopolymer, 8, "must measure the A tract the indel sits in");
        assert_eq!(indel_tract_runs_at(102, &segment).homopolymer, 1, "the C is a run of one");
    }

    fn test_segment(sequence: &[u8]) -> Segment {
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
            overlap_start: 0,
            overlap_end: 0,
        }
    }
}
