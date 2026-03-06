use super::overlapping_reads::{NameCollector, resolve_pair};
use crate::{
    call::{
        pileup::{Pileup, PositionInRead, SimpleRead, SimpleReads},
        process::PileupMappingParams,
    },
    sequence::Segment,
    utils::SequenceContext,
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use rastair_types::{SmallVec, strand_from_flags};
use rust_htslib::bam::pileup::{Alignment, Pileup as HtsPileup};
use std::rc::Rc;
use tracing::{debug, instrument};

impl Pileup {
    /// Convert a pileup from htslib into our internal Pileup representation.
    #[instrument(level = "trace", skip_all)]
    pub(crate) fn from_hts(
        pile: &HtsPileup,
        segment: Rc<Segment>,
        params: &PileupMappingParams,
        collector: &mut NameCollector,
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
            .filter_map(alignment_to_read)
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

fn alignment_to_read<'a>(a: Alignment<'a>) -> Option<(&'a [u8], SimpleRead)> {
    let pos = a.qpos()?;
    let record = a.record_view();
    let flags = record.flags();
    let (matches, indels) = calc_cigar_data(record.raw_cigar());
    let (seq, qual) = record.seq_and_qual();

    Some((
        record.qname(),
        SimpleRead {
            base: seq[pos].into(),
            qual: *qual.get(pos)?,
            mapq: record.mapq(),
            strand: strand_from_flags(flags).ok()?,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::process::PileupMappingParams;
    use crate::call::variant_calling::VariantCallingParams;
    use crate::utils::default;

    #[test]
    fn name_collector_skip_when_keep_overlapping() {
        let params = PileupMappingParams {
            variant_calling: VariantCallingParams { keep_overlapping_reads: true, ..default() },
            require_tags: default(),
        };
        assert!(matches!(NameCollector::new(&params), NameCollector::Skip));
    }

    #[test]
    fn name_collector_collect_by_default() {
        let params = PileupMappingParams { variant_calling: default(), require_tags: default() };
        assert!(matches!(NameCollector::new(&params), NameCollector::Collect(_)));
    }
}
