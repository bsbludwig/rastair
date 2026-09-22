use crate::{
    call::{RecordFilters, require_tags::TagRequirement, variant_calling::VariantCallingParams},
    sequence::{ChunkRegion, Segment},
};
use color_eyre::eyre::Result;
use std::{ops::Deref, rc::Rc};

#[cfg(not(feature = "experimental-seqair"))]
use crate::{
    call::pileup::{Pileup, from_hts::PileupScratch},
    sequence::Readers,
    utils::logging::ThisIsABug,
};
#[cfg(not(feature = "experimental-seqair"))]
use color_eyre::eyre::WrapErr;
#[cfg(not(feature = "experimental-seqair"))]
use rust_htslib::bam::{FetchDefinition, Read as _};
#[cfg(not(feature = "experimental-seqair"))]
use tracing::{Level, debug, instrument, trace, warn};

#[cfg(feature = "experimental-seqair")]
use crate::{
    call::pileup::from_seqair::{ColumnDraft, ColumnScratch, DenovoNeighbour},
    metrics::{PileupMetrics, entropy::SlidingEntropy},
    sequence::{PileupReaders, ReferenceWindow},
};
#[cfg(feature = "experimental-seqair")]
use color_eyre::eyre::WrapErr as _;
#[cfg(feature = "experimental-seqair")]
use seqair::reader::{DepthLimit, SegmentOptions};
#[cfg(feature = "experimental-seqair")]
use seqair_types::{Base, Pos0};
#[cfg(feature = "experimental-seqair")]
use std::{
    num::{NonZeroU32, NonZeroU64},
    sync::Arc,
};
#[cfg(feature = "experimental-seqair")]
use tracing::{debug, instrument, warn};

#[derive(better_default::Default)]
pub struct PileupMappingParams {
    pub variant_calling: VariantCallingParams,
    pub require_tags: TagRequirement,
    pub call_indels: bool,
    /// Ignore indels within this distance from read ends.
    #[default(0)]
    pub indel_end_of_read_cutoff: usize,
    /// Maximum number of non-TAPS mismatches allowed on a read supporting an indel.
    #[default(5)]
    pub indel_max_mismatches: u32,
    /// Per-segment compressed-byte budget for the seqair backend; segments
    /// estimated above this are subdivided. `0` disables the budget.
    #[default(32 * 1024 * 1024)]
    pub segment_max_bytes: u64,
    /// Rescue the soft-clipped CpG-partner base adjacent to each alignment
    /// (seqair backend only). See `MethylationCallingParams::rescue_soft_clip_cpg`.
    pub rescue_soft_clip_cpg: bool,
    /// Which columns are worth finishing (seqair backend only).
    ///
    /// `None` finishes every column the reader emits, for a caller that wants
    /// the raw pileup. Production passes the record filters, so the seven
    /// columns in eight they drop never become a `PileupMetrics` at all.
    pub early_reject: Option<RecordFilters>,
}

impl Deref for PileupMappingParams {
    type Target = VariantCallingParams;

    fn deref(&self) -> &Self::Target {
        &self.variant_calling
    }
}

#[cfg(not(feature = "experimental-seqair"))]
#[instrument(level = "info", skip_all)]
pub fn get_pileups(
    readers: &mut Readers,
    region: &ChunkRegion,
    params: &PileupMappingParams,
) -> Result<(Rc<Segment>, impl Iterator<Item = Pileup>)> {
    let segment = readers.segment(region, 2).wrap_err("Failed to fetch segment")?;
    debug!(len = segment.sequence.len(), "Processing region");

    // Fetch the pileups for the segment
    FetchDefinition::try_from(&segment.region)
        .wrap_err("Could not convert region string")
        .this_is_a_bug()
        .and_then(|r| readers.bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))
        .wrap_err_with(|| format!("Could not fetch region `{}` from BAM file", region.region))?;

    // Install a pileup-level filter so that reads can be rejected once at read
    // time instead of once per pileup column they span.
    {
        let read_flags = params.read_flags.clone();
        let unpaired = params.unpaired;
        let tag_filter = params.require_tags.clone();
        readers.bam.set_pileup_filter(move |record| {
            read_flags.filter_flags(record.flags(), unpaired) && tag_filter.allows(&record)
        });
    }

    let segment = Rc::new(segment);
    let segment_clone = segment.clone();

    // Go over each column in the pileup from htslib and build our own pileup
    let mut pileup = readers.bam.pileup();
    // htslib's own cap only bounds memory (the real cap is applied to filtered
    // reads per column, see `MaxCoverage::load_ceiling`); `0` there means "use
    // the default", not "unlimited", so an unlimited run asks for the largest
    // cap htslib can take instead. `set_max_depth` panics above `i32::MAX`
    // (it hands the value to `bam_plp_set_maxcnt`, which takes a C `int`), so
    // both the unlimited sentinel and a huge `--max-coverage` have to be
    // clamped down to it rather than passed through.
    pileup.set_max_depth(
        params
            .max_coverage
            .load_ceiling()
            .map_or(i32::MAX as u32, std::num::NonZeroU32::get)
            .min(i32::MAX as u32),
    );
    let mut scratch = PileupScratch::new(params);
    let piles = pileup
        .filter_map(|p| match p {
            Ok(p) => Some(p),
            Err(e) => {
                if tracing::enabled!(Level::TRACE) {
                    trace!(%e, "Failed to read pileup, skipping");
                }
                None
            }
        })
        .filter(|p| {
            // We might get pileups from htslib that are not in the region of
            // interest (but actually before it). Since our segments only cover
            // the specified region, we can just skip these (we won't have the
            // reference sequence for this anyway). They'll be part of the
            // next/previous segment anyway.
            region.contains(u64::from(p.pos()))
        })
        .map(move |pile| {
            Pileup::from_hts(&pile, segment.clone(), params, &mut scratch).wrap_err_with(|| {
                format!("Failed to get candidate from pileup at position {}", pile.pos())
            })
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

#[cfg(feature = "experimental-seqair")]
#[instrument(level = "info", skip_all)]
pub fn get_pileups(
    readers: &mut PileupReaders,
    region: &ChunkRegion,
    params: &PileupMappingParams,
) -> Result<(Rc<Segment>, impl Iterator<Item = PileupMetrics>)> {
    // Build the rastair Segment (FASTA fetch only) and extract ref window.
    let segment = readers.segment(region, 2).wrap_err("Failed to fetch segment")?;
    debug!(len = segment.sequence.len(), "Processing region (seqair)");

    let ref_window = {
        let bases: Arc<[Base]> =
            segment.sequence.iter().map(|&b| Base::from(b)).collect::<Vec<_>>().into();
        let start = Pos0::try_from(region.start).wrap_err("region start out of Pos0 range")?;
        ReferenceWindow { bases, start }
    };

    // Install the reference window so compute() can do mismatch inference.
    readers.inner_mut().customize_mut().reference = Some(ref_window);
    readers.inner_mut().customize_mut().read_flags = params.read_flags.clone();
    readers.inner_mut().customize_mut().unpaired = params.unpaired;
    readers.inner_mut().customize_mut().tag_requirement = params.require_tags.clone();
    readers.inner_mut().customize_mut().guess_orientation = params.guess_read_orientation;
    readers.inner_mut().customize_mut().read_masking = params.read_masking.clone();

    // Loading and the engine are bounded by the memory ceiling only; the
    // user's cap counts filtered reads, inside `ColumnDraft::accumulate`.
    let depth_limit = match params.max_coverage.load_ceiling() {
        Some(ceiling) => DepthLimit::PerColumn(ceiling),
        None => DepthLimit::Unlimited,
    };

    // Build the seqair Segment covering [region.start .. region.end + overfetch].
    // seqair's segment `end` is an *inclusive* 0-based position and must be
    // `<= contig_last_pos` (= contig length − 1). `region.last_position` is the
    // contig length, so clamp to `last_position - 1`; htslib tolerates ends one
    // past the contig, seqair rejects them.
    let overfetch: u64 = 2;
    let contig_last_pos = region.last_position.saturating_sub(1);
    let last = region.end.saturating_add(overfetch).min(contig_last_pos);
    let start = Pos0::try_from(region.start).wrap_err("region start out of Pos0 range")?;
    let end = Pos0::try_from(last).wrap_err("region end out of Pos0 range")?;
    let len_u32 =
        u32::try_from(last.saturating_sub(region.start).saturating_add(1)).unwrap_or(u32::MAX);
    const ONE: NonZeroU32 = NonZeroU32::MIN;
    let max_len = NonZeroU32::new(len_u32.max(1)).unwrap_or(ONE);
    // Bound peak memory per worker: split this region into sub-segments whose
    // estimated compressed load stays within the byte budget (`0` disables it).
    let opts = match NonZeroU64::new(params.segment_max_bytes) {
        Some(budget) => SegmentOptions::new(max_len).with_max_bytes(budget),
        None => SegmentOptions::new(max_len).without_byte_budget(),
    };
    // Collect the (possibly several) sub-segments before piling up: `segments()`
    // borrows the reader immutably while `pileup()` needs it mutably. The
    // sub-segments are disjoint (overlap 0), but seqair fetches every read
    // overlapping each one, so no boundary reads are lost and no column is
    // emitted twice.
    let seqair_segments: Vec<_> = readers
        .inner_mut()
        .segments((region.contig.as_str(), (start..=end).into()), opts)
        .wrap_err("Failed to plan seqair segments")?
        .collect();

    let segment = Rc::new(segment);
    // `early_reject` drops seven covered columns in eight before anything
    // downstream sees them, so this holds the survivors and not one entry per
    // covered position. A quarter of the region is a deliberate
    // over-reservation — a `PileupMetrics` is 928 bytes and doubling into it
    // both copies and overshoots — against a chr12 keep rate of ~12 % without
    // `--cpgs-only` and far below that with it. A caller that asked for every
    // column doubles up into the rest, which costs it two copies of a vector
    // that was going to be that size anyway.
    //
    // `region.len()` has to be the *inclusive* count for this to land where it
    // is meant to.
    let mut pileup_metrics: Vec<PileupMetrics> =
        Vec::with_capacity(usize::try_from(region.len() / 4).unwrap_or(0));
    // Reused across every column of every sub-segment; see `ColumnDraft::accumulate`.
    let mut scratch = ColumnScratch::default();

    // The one-column delay that lets a column be rejected before it is
    // finished. `set_denovo_adj` lets the *next* emitted column rescue a
    // reference C, so such a draft waits here until that column is known;
    // every other verdict is settled by what precedes it. The draft carries
    // the `before` it was judged against, because by the time its successor
    // arrives `previous` has moved on to the draft itself.
    //
    // Both live outside the sub-segment loop because the sequence they
    // describe — this region's emitted columns, in ascending position — spans
    // sub-segments, just as the vector this replaces did.
    let mut deferred: Option<(ColumnDraft, Option<DenovoNeighbour>)> = None;
    let mut previous: Option<DenovoNeighbour> = None;
    let keeps = |inputs| params.early_reject.as_ref().is_none_or(|f| f.keeps(inputs));
    let mut sliding_entropy = SlidingEntropy::new(&segment);

    for seqair_seg in &seqair_segments {
        // Fetch BAM records + FASTA into PileupEngine (compute() runs here).
        // The store is reused (cleared) per sub-segment, so peak memory is one
        // sub-segment's worth, not the whole region's.
        let mut guard = readers
            .inner_mut()
            .pileup(seqair_seg, depth_limit)
            .run()
            .wrap_err("Failed to start seqair pileup")?;
        if let Some(ceiling) = params.max_coverage.load_ceiling() {
            guard.set_max_depth(ceiling);
        }
        if params.rescue_soft_clip_cpg {
            // Recover exactly the single CpG-partner base the aligner clipped.
            guard.set_soft_clip_overhang(1);
        }

        while let Some(col) = guard.pileups() {
            let pos = col.pos().as_u64();
            if !region.contains(pos) {
                continue;
            }
            let draft = match ColumnDraft::accumulate(&col, segment.clone(), params, &mut scratch) {
                Ok(draft) => draft,
                Err(error) => {
                    warn!(error = format!("{error:#}"), pos, "Failed to get pileup, skipping");
                    continue;
                }
            };
            let neighbour = draft.denovo_neighbour();

            // The deferred column comes first: output stays in ascending
            // position order.
            if let Some((waiting, before)) = deferred.take()
                && keeps(waiting.pre_filter_inputs(before, Some(neighbour)))
            {
                keep(&mut pileup_metrics, &mut sliding_entropy, waiting);
            }

            if keeps(draft.pre_filter_inputs(previous, None)) {
                keep(&mut pileup_metrics, &mut sliding_entropy, draft);
            } else if draft.awaits_successor() {
                deferred = Some((draft, previous));
            }
            previous = Some(neighbour);
        }
    }

    // Nothing follows the last column, so a still-deferred draft was only ever
    // waiting for a rescue that cannot come.
    if let Some((waiting, before)) = deferred.take()
        && keeps(waiting.pre_filter_inputs(before, None))
    {
        keep(&mut pileup_metrics, &mut sliding_entropy, waiting);
    }

    Ok((segment, pileup_metrics.into_iter()))
}

/// Finish a surviving column and give it its region entropy.
///
/// The entropy is only read by the ML feature extractor, so it is computed
/// here and not for every covered position; the sliding window is still fed
/// ascending indices and its counts are integers, so the values are the ones
/// the whole-region pass produced.
#[cfg(feature = "experimental-seqair")]
fn keep(
    pileup_metrics: &mut Vec<PileupMetrics>,
    sliding_entropy: &mut SlidingEntropy<'_>,
    draft: ColumnDraft,
) {
    let idx = draft.idx();
    match draft.finish() {
        Ok(mut metrics) => {
            metrics.pos_metrics.extended.region_entropy = sliding_entropy.entropy_at(idx);
            pileup_metrics.push(metrics);
        }
        Err(error) => {
            warn!(error = format!("{error:#}"), "Failed to finish pileup, skipping");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::ReaderParams;
    use color_eyre::eyre::Result;

    #[test]
    fn test_reading_bounds() -> Result<()> {
        // check that we can read exactly the right positions
        let params = ReaderParams {
            regions: Some("chr19:6105700-6105800".parse().unwrap()),
            ..ReaderParams::test_data()
        };
        let mut readers = params.pileup_readers()?;
        let segments: Vec<_> = readers.segments(10_000, 100)?.collect();
        readers.segment(&segments[0], 0)?;

        let pileup_mapping_params = PileupMappingParams::default();
        let (_segment, pileups) = get_pileups(&mut readers, &segments[0], &pileup_mapping_params)?;
        let pileups: Vec<_> = pileups.collect();

        assert!(!pileups.is_empty());
        assert_eq!(pileups.first().unwrap().pos, 6_105_700);
        assert_eq!(pileups.last().unwrap().pos, 6_105_800);

        Ok(())
    }

    /// The early rejection is exactly [`RecordFilters::pre_filter`] moved in
    /// front of the work it saves. Both halves of that claim are checked: the
    /// surviving positions are the ones the downstream `retain` would have
    /// kept, in the same order, and their region entropy is what the
    /// whole-region sliding pass produced — the window only ever sees the kept
    /// columns now.
    ///
    /// `set_denovo_adj` is the reason this is not obvious: a column with no
    /// evidence of its own is kept when its genomic neighbour carries an alt
    /// that would create the other half of a CpG, so the reference run has to
    /// go through `map_surrounding` before it filters.
    #[cfg(feature = "experimental-seqair")]
    #[test]
    fn early_rejection_keeps_exactly_what_the_pre_filter_would() -> Result<()> {
        use crate::{call::RecordFilters, utils::map_surrounding};

        let params = ReaderParams {
            regions: Some("chr19:6105700-6106500".parse()?),
            ..ReaderParams::test_data()
        };
        let mut readers = params.pileup_readers()?;
        let segments: Vec<_> = readers.segments(10_000, 100)?.collect();

        let mut compared = 0usize;
        for filters in [
            RecordFilters { vcf_all: false, cpgs_only: false },
            RecordFilters { vcf_all: true, cpgs_only: false },
            RecordFilters { vcf_all: false, cpgs_only: true },
        ] {
            let unfiltered = PileupMappingParams::default();
            for chunk in &segments {
                let (_s, everything) = get_pileups(&mut readers, chunk, &unfiltered)?;
                let mut everything: Vec<_> = everything.collect();
                map_surrounding(&mut everything, super::super::set_denovo_adj, "test mapper");
                everything.retain(|p| filters.pre_filter(p));
                let expected: Vec<_> = everything
                    .iter()
                    .map(|p| (p.pos, p.pos_metrics.extended.region_entropy))
                    .collect();

                let early = PileupMappingParams {
                    early_reject: Some(filters.clone()),
                    ..Default::default()
                };
                let (_s, kept) = get_pileups(&mut readers, chunk, &early)?;
                let actual: Vec<_> =
                    kept.map(|p| (p.pos, p.pos_metrics.extended.region_entropy)).collect();

                assert_eq!(expected, actual, "early rejection diverged for {filters:?}");
                compared += expected.len();
            }
        }
        assert!(compared > 0, "the fixture produced nothing to compare");
        Ok(())
    }

    /// Byte-aware sub-segmentation must be transparent: a tiny budget (which
    /// forces the region to be split into many sub-segments) yields exactly the
    /// same pileups as an effectively-unlimited budget. Only meaningful on the
    /// seqair backend; on htslib the budget is ignored and both runs are
    /// trivially identical.
    #[cfg(feature = "experimental-seqair")]
    #[test]
    fn byte_budget_subdivision_is_transparent() -> Result<()> {
        let params = ReaderParams {
            regions: Some("chr19:6105700-6105900".parse().unwrap()),
            ..ReaderParams::test_data()
        };
        let mut readers = params.pileup_readers()?;
        let segments: Vec<_> = readers.segments(10_000, 100)?.collect();

        // Everything a boundary read could move: depth, the ALT list *in order*
        // (it becomes the VCF ALT column), per-allele depth and strand split,
        // and the indel observation count. A read fetched twice, or dropped at
        // a sub-segment edge, shows up in one of these. Types stay inferred so
        // this compiles on both the htslib and seqair backends.
        #[allow(clippy::type_complexity)]
        fn fingerprint(
            pileups: impl Iterator<Item = PileupMetrics>,
        ) -> Vec<(u32, Base, usize, Vec<(Base, u32, u32, u32)>, usize)> {
            pileups
                .map(|p| {
                    let alleles = std::iter::once(&p.ref_metrics)
                        .chain(p.alts.iter().map(|alt| &alt.metrics))
                        .map(|m| (m.base, m.depth, m.strand_count.ot, m.strand_count.ob))
                        .collect();
                    let indels = p.indel_data.as_ref().map_or(0, |data| data.observations.len());
                    (p.pos, p.reference_base, p.pos_metrics.depth as usize, alleles, indels)
                })
                .collect()
        }

        let huge = PileupMappingParams { segment_max_bytes: u64::MAX, ..Default::default() };
        let (_s1, p1) = get_pileups(&mut readers, &segments[0], &huge)?;
        let baseline = fingerprint(p1);

        // A 1-byte budget forces the region to split into many sub-segments.
        let tiny = PileupMappingParams { segment_max_bytes: 1, ..Default::default() };
        let (_s2, p2) = get_pileups(&mut readers, &segments[0], &tiny)?;
        let subdivided = fingerprint(p2);

        assert!(!baseline.is_empty());
        assert!(
            baseline.iter().any(|(_, _, _, alleles, _)| alleles.len() > 1),
            "the region must contain alts, or the ALT comparison proves nothing"
        );
        assert_eq!(baseline, subdivided, "subdivision changed the pileups");
        Ok(())
    }
}
