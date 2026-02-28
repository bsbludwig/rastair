use crate::{
    call::{
        pileup::{Pileup, PositionInRead, SimpleRead, SimpleReads},
        process::PileupMappingParams,
    },
    sequence::Segment,
    utils::SequenceContext,
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use rastair_types::strand_from_flags;
use rust_htslib::bam::pileup::{Alignment, Pileup as HtsPileup};
use std::rc::Rc;
use tracing::{debug, instrument};

impl Pileup {
    #[instrument(level = "trace", skip_all)]
    pub fn from_hts(
        pile: &HtsPileup,
        segment: Rc<Segment>,
        params: &PileupMappingParams,
    ) -> Result<Pileup> {
        let pos = pile.pos();
        let idx = segment.pos_to_idx(pos)?;
        let depth = pile.depth();
        let max_reads = depth.min(params.max_coverage);
        if depth > max_reads {
            debug!(pos, depth, "Capping number of reads in pileup to {max_reads}");
        }
        let max_reads = usize::try_from(max_reads).wrap_err("max_reads exceeds usize")?;

        let mut reads = Vec::with_capacity(max_reads);
        let mut names: Vec<&[u8]> = Vec::with_capacity(max_reads);

        pile.alignments()
            .filter_map(|a| alignment_to_read(params, a))
            .filter(|(_, seen_base)| params.read_masking.filter(seen_base))
            .filter(|(_, seen_base)| params.quality.filter(seen_base))
            .take(max_reads)
            .for_each(|(name, seen_base)| {
                names.push(name);
                reads.push(seen_base);
            });

        let mut reads = SimpleReads(reads);

        if should_deduplicate_overlapping_pairs(params) {
            reads.remove_overlapping_pairs(&names);
        }

        let reference_base =
            segment.sequence.get(idx).wrap_err("failed to get reference base")?.into();

        let context =
            SequenceContext::new(idx, &segment).wrap_err("failed to get sequence context")?;

        Ok(Pileup {
            region: segment.range.clone(),
            context,
            pos: pile.pos(),
            reads,
            reference_base,
        })
    }
}

fn should_deduplicate_overlapping_pairs(params: &PileupMappingParams) -> bool {
    !params.unpaired && !params.keep_overlapping_reads
}

/// Collect info from a pileup alignment
fn alignment_to_read<'a>(
    params: &PileupMappingParams,
    a: Alignment<'a>,
) -> Option<(&'a [u8], SimpleRead)> {
    let pos = a.qpos()?;
    let record = a.record_view();
    let flags = record.flags();

    if !params.read_flags.filter_flags(flags, params.unpaired) {
        return None;
    }
    let (matches, indels) = calc_cigar_data(record.raw_cigar());

    Some((
        record.qname(),
        SimpleRead {
            base: record.seq()[pos].into(),
            qual: *record.qual().get(pos)?,
            mapq: record.mapq(),
            // Strand of the read, derived from flags. Early return if strand cannot be determined.
            // TODO: handle "lenient mode"
            strand: strand_from_flags(flags).ok()?,
            reverse: record.is_reverse(),
            second: record.is_last_in_template(),
            position: PositionInRead {
                pos: u32::try_from(pos).expect("position fits in u32"),
                read_length: u32::try_from(record.seq_len()).expect("read length fits in u32"),
            },
            matching_bases: matches,
            indels,
        },
    ))
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

#[cfg(test)]
mod tests {
    use super::should_deduplicate_overlapping_pairs;
    use crate::call::{process::PileupMappingParams, variant_calling::VariantCallingParams};

    #[test]
    fn deduplicates_overlapping_pairs_by_default() {
        let params = PileupMappingParams { variant_calling: VariantCallingParams::default() };
        assert!(should_deduplicate_overlapping_pairs(&params));
    }

    #[test]
    fn does_not_deduplicate_when_keep_overlapping_reads_is_enabled() {
        let mut variant_calling = VariantCallingParams::default();
        variant_calling.keep_overlapping_reads = true;
        let params = PileupMappingParams { variant_calling };
        assert!(!should_deduplicate_overlapping_pairs(&params));
    }

    #[test]
    fn does_not_deduplicate_in_unpaired_mode() {
        let mut variant_calling = VariantCallingParams::default();
        variant_calling.unpaired = true;
        let params = PileupMappingParams { variant_calling };
        assert!(!should_deduplicate_overlapping_pairs(&params));
    }

    #[test]
    fn unpaired_mode_takes_precedence_over_dedup_setting() {
        let mut variant_calling = VariantCallingParams::default();
        variant_calling.unpaired = true;
        variant_calling.keep_overlapping_reads = false;
        let params = PileupMappingParams { variant_calling };
        assert!(!should_deduplicate_overlapping_pairs(&params));
    }
}
