use crate::{
    call::{
        pileup::{Pileup, PositionInRead, SimpleRead, SimpleReads},
        process::PileupMappingParams,
    },
    sequence::Segment,
    utils::{SequenceContext, StrandFromRecord, logging::ThisIsABug},
    vcf::InCpG,
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use rust_htslib::bam::pileup::{Alignment, Pileup as HtsPileup};
use smallvec::SmallVec;
use std::sync::Arc;
use tracing::instrument;

impl Pileup {
    #[instrument(level = "trace", skip_all)]
    pub fn from_hts(
        pile: &HtsPileup,
        segment: Arc<Segment>,
        params: &PileupMappingParams,
    ) -> Result<Pileup> {
        let pos = pile.pos();
        let segment_start_pos = usize::try_from(segment.range.start)
            .wrap_err("segment range fits in usize")
            .this_is_a_bug()?;
        let idx = usize::try_from(pos)
            .wrap_err("position fits in usize")
            .this_is_a_bug()?
            .checked_sub(segment_start_pos)
            .wrap_err_with(|| {
                format!("pile position {pos} is not in segment {}", segment.region)
            })?;

        let seen_bases = pile
            .alignments()
            .filter_map(|pile| alignment_to_read(params, pile))
            .filter(|seen_base| params.read_masking.filter(seen_base))
            .filter(|seen_base| params.quality.filter(seen_base))
            .collect();

        let mut reads = SimpleReads(seen_bases);

        if !params.keep_overlapping_reads {
            reads.remove_overlapping_pairs();
        }

        let reference_base =
            segment.sequence.get(idx).wrap_err("failed to get reference base")?.into();

        let context = SequenceContext::new(reference_base, idx, &segment)
            .wrap_err("failed to get sequence context")?;
        let is_cpg = *InCpG::new(reference_base, context.before_1, context.after_1);

        Ok(Pileup {
            region: segment.range.clone(),
            context,
            pos: pile.pos(),
            reads,
            reference_base,
            is_cpg,
        })
    }
}

/// Collect info from a pileup alignment
fn alignment_to_read(params: &PileupMappingParams, a: Alignment<'_>) -> Option<SimpleRead> {
    let pos = a.qpos()?;
    let record = a.record();
    let cigar = record.raw_cigar();

    if !params.read_flags.filter(&record) {
        return None;
    }
    let (matches, indels) = calc_cigar_data(cigar);

    Some(SimpleRead {
        qname: SmallVec::from(record.qname()),
        base: record.seq()[pos].into(),
        qual: *record.qual().get(pos)?,
        mapq: record.mapq(),
        // Strand of the read, derived from the record. Early return if strand cannot be determined.
        // TODO: handle "lenient mode"
        strand: StrandFromRecord::strand(&record).ok()?,
        reverse: record.is_reverse(),
        second: record.is_last_in_template(),
        position: PositionInRead {
            pos: u32::try_from(pos).expect("position fits in u32"),
            read_length: u32::try_from(record.seq_len()).expect("read length fits in u32"),
        },
        matching_bases: matches,
        indels,
    })
}

/// Calculate the number of matches and indels from a packed CIGAR array.
///
/// Packed CIGAR data is encoded as follows:
/// - lower 4 bits for the operation
/// - upper 28 bits for the length
fn calc_cigar_data(cigar: &[u32]) -> (u32, u32) {
    let mut matches = 0;
    let mut indels = 0;
    for c in cigar {
        let len = c >> 4;
        match c & 0b1111 {
            // Match
            0 => matches += len,
            // Insertion or deletion
            1 | 2 => indels += len,
            _ => {
                // Other operations (like soft clipping, padding, etc.) are ignored
                // for the purpose of counting matches and indels.
            }
        }
    }
    (matches, indels)
}
