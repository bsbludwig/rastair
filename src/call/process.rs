use crate::{
    call::{
        pileup::{Pileup, PositionInRead, SimpleRead, SimpleReads},
        variant_calling::VariantCallingParams,
    },
    sequence::{ChunkRegion, Readers, Segment},
    utils::{SequenceContext, StrandFromRecord},
    vcf::{self, Filters, InCpG},
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use rust_htslib::bam::{
    FetchDefinition, Read as _,
    pileup::{Alignment, Pileup as HtsPileup},
};
use smallvec::SmallVec;
use std::{ops::Deref, sync::Arc};
use tracing::{Level, debug, instrument, trace, warn};

#[derive(Debug, Clone)]
pub struct PileupMappingParams {
    pub include_cpgs: IncludeAllCpGs,
    pub variant_calling: VariantCallingParams,
}

impl Deref for PileupMappingParams {
    type Target = VariantCallingParams;

    fn deref(&self) -> &Self::Target {
        &self.variant_calling
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IncludeAllCpGs {
    Yes,
    No,
}

impl Deref for IncludeAllCpGs {
    type Target = bool;

    fn deref(&self) -> &Self::Target {
        match self {
            IncludeAllCpGs::Yes => &true,
            IncludeAllCpGs::No => &false,
        }
    }
}

impl ChunkRegion {
    /// Process the chunk region to collect pileups
    ///
    /// # Returns
    /// - The segment corresponding to the chunk region
    /// - An iterator over pileups in the region
    #[instrument(level = "info", skip_all)]
    pub fn process(
        &self,
        readers: &mut Readers,
        params: &PileupMappingParams,
    ) -> Result<(Arc<Segment>, impl Iterator<Item = Pileup>)> {
        let segment = readers.segment(self, 2).wrap_err("failed to fetch segment")?;
        trace!(len = segment.sequence.len(), "Processing region");

        // Fetch the pileups for the segment
        FetchDefinition::try_from(&segment.region)
            .wrap_err("Could not convert region string")
            .and_then(|r| readers.bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))
            .wrap_err_with(|| format!("Could not fetch region `{}` from BAM file", self.region))?;

        let segment = Arc::new(segment);
        let segment_clone = segment.clone();

        // Go over each column in the pileup and collect variant candidates
        let piles = readers
            .bam
            .pileup()
            .filter_map(|p| {
                if tracing::enabled!(Level::TRACE) {
                    match p {
                        Ok(p) => Some(p),
                        Err(e) => {
                            trace!(%e, "Failed to read pileup, skipping");
                            None
                        }
                    }
                } else {
                    p.ok()
                }
            })
            .filter(|p| {
                // Filter out pileups that are not in the region of interest
                self.contains(u64::from(p.pos()))
            })
            .flat_map(move |pile| {
                collect_candidate(&pile, segment.clone(), params)
                    .wrap_err_with(|| {
                        format!("Failed to get candidate from pileup at position {}", pile.pos())
                    })
                    .transpose()
            })
            .filter_map(|res| match res {
                Ok(x) => Some(x),
                Err(error) => {
                    warn!(error = format!("{error:#}"), "Failed to get pileup, skipping");
                    None
                }
            });
        Ok((segment_clone, piles))
    }
}

/// Is this pileup a candidate for a variant?
#[instrument(level = "trace", skip_all)]
fn collect_candidate(
    pile: &HtsPileup,
    segment: Arc<Segment>,
    params: &PileupMappingParams,
) -> Result<Option<Pileup>> {
    let segment_start_pos =
        usize::try_from(segment.range.start).expect("segment range fits in usize");
    let idx = usize::try_from(pile.pos())
        .expect("position fits in usize")
        .checked_sub(segment_start_pos)
        .wrap_err_with(|| {
            format!(
                "pile position {} is not in segment range {}..{}",
                pile.pos(),
                segment.range.start,
                segment.range.end
            )
        })?;

    let seen_bases = pile
        .alignments()
        .filter_map(|pile| pileup_mapper(params, pile))
        .filter(|seen_base| params.read_masking.filter(seen_base))
        .filter(|seen_base| params.quality.filter(seen_base))
        .collect();

    let mut reads = SimpleReads(seen_bases);

    if !params.keep_overlapping_reads {
        reads.remove_overlapping_pairs();
    }

    let reference_base = segment.sequence.get(idx).wrap_err("failed to get reference base")?.into();

    let context = SequenceContext::new(reference_base, idx, &segment)
        .wrap_err("failed to get sequence context")?;
    let is_cpg = *InCpG::new(reference_base, context.before_1, context.after_1);
    let has_alts = !reads.matches(reference_base);

    let res = Pileup {
        region: segment.range.clone(),
        context,
        pos: pile.pos(),
        reads,
        reference_base,
        is_cpg,
    };

    if has_alts || (*params.include_cpgs && res.is_cpg) {
        Ok(Some(res))
    } else {
        // Matches reference base, boring.
        // trace!(?bases, pos = pile.pos(), ?reference_base, ?next_base, "pile matches reference");
        Ok(None)
    }
}

/// Collect info from a pileup alignment
pub(crate) fn pileup_mapper(params: &PileupMappingParams, a: Alignment<'_>) -> Option<SimpleRead> {
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

impl Pileup {
    /// Collect metrics
    #[instrument(level = "trace", skip_all)]
    #[deprecated = "Use `PileupMetrics` instead"]
    pub fn variant_metrics(&self, params: &VariantCallingParams) -> Result<vcf::Record> {
        let metrics = self.metrics().wrap_err("Failed to calculate metrics")?;
        let calling_metrics = self
            .calling_metrics(params.error_model)
            .wrap_err("Failed to calculate calling metrics")?;

        Ok(vcf::Record {
            main: self.fixed_fields(),
            filters: Filters::new(),
            info: metrics,
            samples: smallvec::smallvec![calling_metrics],
        })
    }
}

impl SimpleReads {
    /// Remove overlapping reads from the same fragment.
    pub fn remove_overlapping_pairs(&mut self) {
        // For each read, check if we already saw one with the same name.
        //
        // If the bases agree, keep only the first one. If they disagree, keep none.
        //
        // But this is rust -- so we can't just remove elements while iterating.
        // Instead, we keep a little list of indices to remove, and then remove them afterwards.
        // This should be fine since the amount of items to remove is typically small.
        let mut to_remove = SmallVec::<usize, 16>::new();
        for i in 0..self.0.len() {
            let base_i = &self.0[i];
            for j in (i + 1)..self.0.len() {
                let base_j = &self.0[j];
                if base_i.qname == base_j.qname {
                    // Same read name
                    if base_i.base == base_j.base {
                        // Same base, keep only the first one
                        to_remove.push(j);
                    } else {
                        // Different bases, ignore the second in pair
                        // NOTE: This is different from rastair1
                        if base_i.second {
                            to_remove.push(i);
                        } else {
                            to_remove.push(j);
                        }
                    }
                    // No need to check further
                    break;
                }
            }
        }
        // Remove duplicates
        to_remove.sort_unstable();
        for &idx in to_remove.iter().rev() {
            self.0.swap_remove(idx);
        }
    }
}
