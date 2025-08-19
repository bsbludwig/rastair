use crate::{
    call::{
        variant_calling::{ReadFlags, ReadMaskParams, VariantCallingParams},
        variants::{PositionInRead, SeenBase, SeenBases, VariantCandidatePileup},
    },
    sequence::{ChunkRegion, Readers, Segment},
    utils::{Base, ReadDeduplicator, StrandFromRecord},
    vcf::{self, Filters, InCpG},
};
use color_eyre::eyre::{ContextCompat as _, Result, WrapErr};
use rust_htslib::bam::{
    FetchDefinition, Read as _,
    pileup::{Alignment, Pileup},
};
use std::{ops::Deref, rc::Rc};
use tracing::{Level, debug, instrument, trace, warn};

#[derive(Debug, Clone)]
pub struct PileupMappingParams {
    pub include_cpgs: IncludeAllCpGs,
    pub keep_overlapping_reads: bool,
    pub read_masking: ReadMaskParams,
    pub read_flags: ReadFlags,
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
    #[instrument(level = "info", skip_all)]
    pub fn process(
        &self,
        readers: &mut Readers,
        params: &PileupMappingParams,
    ) -> Result<Vec<VariantCandidatePileup>> {
        let segment = readers.segment(self, 2).wrap_err("failed to fetch segment")?;
        trace!(len = segment.sequence.len(), "Processing region");

        // Fetch the pileups for the segment
        FetchDefinition::try_from(&segment.region)
            .wrap_err("Could not convert region string")
            .and_then(|r| readers.bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))
            .wrap_err_with(|| format!("Could not fetch region `{}` from BAM file", self.region))?;

        let segment = Rc::new(segment);

        // Allocate a set to track seen read names with enough capacity to avoid reallocation.
        let mut resusable_seen_names_set = ReadDeduplicator::with_capacity(64);

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
            .flat_map(|pile| {
                collect_candidate(&pile, segment.clone(), params, &mut resusable_seen_names_set)
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
            })
            .collect::<Vec<_>>();

        if tracing::enabled!(Level::DEBUG) {
            if piles.is_empty() {
                trace!("No candidate piles found in region, skipping");
                return Ok(piles);
            } else {
                let count = readable::num::Unsigned::from(piles.len());
                let bytes = readable::byte::Byte::from(
                    piles.len() * std::mem::size_of::<VariantCandidatePileup>(),
                );
                debug!(%count, %bytes, "Collected candidates");
            }
        }

        Ok(piles)
    }
}

/// Is this pileup a candidate for a variant?
#[instrument(level = "trace", skip_all)]
fn collect_candidate(
    pile: &Pileup,
    segment: Rc<Segment>,
    params: &PileupMappingParams,
    // This set is used to track seen read names. It is reused across calls to
    // avoid reallocation.
    resusable_seen_names_set: &mut ReadDeduplicator,
) -> Result<Option<VariantCandidatePileup>> {
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

    let seen_names = {
        resusable_seen_names_set.clear();
        resusable_seen_names_set
    };

    let seen_bases = pile
        .alignments()
        .filter(|pile| {
            params.keep_overlapping_reads || !seen_names.is_duplicate(pile.record().qname())
        })
        .filter_map(|pile| pileup_mapper(params, pile))
        .filter(|seen_base| params.read_masking.filter(seen_base))
        .collect();

    let bases = SeenBases(seen_bases);
    let reference_base = segment.sequence.get(idx).wrap_err("failed to get reference base")?.into();
    let has_alts = !bases.matches(reference_base);

    let before = idx.checked_sub(1).and_then(|idx| segment.sequence.get(idx)).map(Base::from);
    let after = idx.checked_add(1).and_then(|idx| segment.sequence.get(idx)).map(Base::from);
    let res = VariantCandidatePileup {
        segment: segment.clone(),
        pos: pile.pos(),
        bases,
        reference_base,
        is_cpg: *InCpG::new(reference_base, before, after),
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
pub(crate) fn pileup_mapper(params: &PileupMappingParams, a: Alignment<'_>) -> Option<SeenBase> {
    let pos = a.qpos()?;
    let record = a.record();
    let cigar = record.raw_cigar();

    if !params.read_flags.filter(&record) {
        return None;
    }
    let (matches, indels) = calc_cigar_data(cigar);

    // if !record.is_proper_pair() {
    //     // fixme: maybe be more lenient here
    //     return None;
    // }
    // if record.is_quality_check_failed() {
    //     return None;
    // }
    // fixme: understand this better:
    // if record.cigar().iter().any(|c| matches!(c, Cigar::SoftClip(_))) {
    //     return None;
    // }

    Some(SeenBase {
        // qname: SmallVec::from(record.qname()),
        // fixme: handle error or at least check usual error modes
        base: record.seq()[pos].into(),
        qual: *record.qual().get(pos)?,
        mapq: record.mapq(),
        // Strand of the read, derived from the record. Early return if strand cannot be determined.
        strand: StrandFromRecord::strand(&record).ok()?,
        reverse: record.is_reverse(),
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

impl VariantCandidatePileup {
    /// Collect metrics
    #[instrument(level = "trace", skip_all)]
    pub fn variant_metrics(&self, params: &VariantCallingParams) -> Result<vcf::Record> {
        let metrics = self.metrics().wrap_err("Failed to calculate metrics")?;
        let calling_metrics =
            self.calling_metrics(params).wrap_err("Failed to calculate calling metrics")?;

        Ok(vcf::Record {
            main: self.fixed_fields(),
            filters: Filters::new(),
            info: metrics,
            samples: smallvec::smallvec![calling_metrics],
        })
    }
}
