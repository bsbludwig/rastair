#![cfg(not(feature = "experimental-seqair"))]

use super::overlapping_reads::{DedupInfo, NameCollector, resolve_pair};
use crate::{call::pileup::hts_utils::*, sequence::Segment};
use crate::{
    call::{
        pileup::{Pileup, PositionInRead, SimpleRead, SimpleReads},
        process::PileupMappingParams,
    },
    utils::SequenceContext,
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use rust_htslib::bam::pileup::{Alignment, Indel, Pileup as HtsPileup};
use seqair_types::Base;
use seqair_types::SmallVec;
use std::rc::Rc;
use tracing::{debug, instrument};

/// Buffers and caches reused across every pileup column in a segment.
///
/// A whole-genome run calls [`Pileup::from_hts`] once per reference base, so
/// anything allocated per call is allocated billions of times; these clear and
/// refill instead.
pub(crate) struct PileupScratch {
    names: NameCollector,
    orientation: ReadOrientationCache,
}

#[cfg(not(feature = "experimental-seqair"))]
impl PileupScratch {
    pub(crate) fn new(params: &PileupMappingParams) -> Self {
        Self { names: NameCollector::new(params), orientation: ReadOrientationCache::default() }
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
        let depth = usize::try_from(depth).wrap_err("pileup depth exceeds usize")?;
        // The cap counts reads that pass the filters below, so a MAPQ 0 pileup
        // cannot crowd out the usable reads; `clamp` only sizes the buffers.
        let capacity = params.max_coverage.clamp(depth);
        let max_reads = params.max_coverage.kept();
        if depth > max_reads {
            debug!(pos, depth, "Capping number of counted reads in pileup to {max_reads}");
        }

        let mut raw_reads = Vec::with_capacity(capacity);
        let mut to_remove = SmallVec::<usize, 16>::new();

        scratch.names.prepare(capacity);

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

            let info = DedupInfo { idx: raw_reads.len(), base: read.base, second: read.second };
            let duplicate = scratch.names.see(name, info);

            raw_reads.push(read);
            if let Some(other) = duplicate {
                resolve_pair(&info, info.idx, &other, other.idx, &mut to_remove);
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

        let reference_base: Base =
            segment.sequence.get(idx).wrap_err("failed to get reference base")?.into();

        let context =
            SequenceContext::new(idx, &segment).wrap_err("failed to get sequence context")?;

        Ok(Pileup { region: segment.range.clone(), context, pos, reads, reference_base })
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::{process::PileupMappingParams, variant_calling::VariantCallingParams};
    use crate::sequence::{ChunkRegion, Region, Segment};
    use crate::utils::default;
    use rust_htslib::bam::Record;
    use rust_htslib::bam::record::{Cigar, CigarString};
    use seqair_types::Strand;

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
