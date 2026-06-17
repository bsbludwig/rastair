use crate::{
    sequence::{
        ChunkRegion, Region, SelectedRegion, chunked::ChunkedRegions, segementation::Segment,
    },
    utils::{
        CliRegionInput, RegionString, cli,
        file_helpers::{FastaReader, open_fasta},
    },
};
use clap::value_parser;
use clio::ClioPath;
use color_eyre::{
    Result, Section,
    eyre::{Context, ContextCompat, ensure},
};
use rust_htslib::bam::{self, HeaderView, Read as _};
use seqair_types::SmolStr;
use tracing::{debug, instrument, trace};

// ── seqair path ───────────────────────────────────────────────────────────────

#[cfg(feature = "experimental-seqair")]
pub use seqair_readers::{RastairReadExtras, ReferenceWindow, SeqairReaders};

#[cfg(feature = "experimental-seqair")]
mod seqair_readers {
    use crate::{
        call::{require_tags::TagRequirement, variant_calling::ReadFlags},
        sequence::{
            ChunkRegion, Region, SelectedRegion, chunked::ChunkedRegions, segementation::Segment,
        },
        utils::{CliRegionInput, RegionString},
    };
    use color_eyre::eyre::{ContextCompat as _, Result, WrapErr as _, ensure};
    use seqair::bam::cigar::CigarOpType;
    use seqair::bam::record_store::{
        CustomizeRecordStore, FilterRawFields, RecordStore, SlimRecord,
    };
    use seqair_types::{Base, Pos0, SmallVec, Strand, strand_from_flags};
    use std::sync::Arc;
    use tracing::{debug, instrument};

    /// Reference sequence window installed by the pileup driver before each
    /// `readers.pileup()` call so that `RastairRecordFilter::compute` can do
    /// mismatch-motif inference and deletion reference-base lookup.
    #[derive(Clone)]
    pub struct ReferenceWindow {
        pub bases: Arc<[Base]>,
        pub start: Pos0,
    }

    impl Default for ReferenceWindow {
        fn default() -> Self {
            Self { bases: Arc::from([] as [Base; 0]), start: Pos0::ZERO }
        }
    }

    impl ReferenceWindow {
        pub fn base_at(&self, pos: u32) -> Option<Base> {
            let offset = usize::try_from(u64::from(pos).checked_sub(self.start.as_u64())?).ok()?;
            self.bases.get(offset).copied()
        }

        pub fn range(&self, pos: u32, len: u32) -> Option<&[Base]> {
            let offset = usize::try_from(u64::from(pos).checked_sub(self.start.as_u64())?).ok()?;
            let end = offset.checked_add(usize::try_from(len).ok()?)?;
            self.bases.get(offset..end)
        }
    }

    /// Per-record data computed once at fetch time.
    #[derive(Clone)]
    pub struct RastairReadExtras {
        pub strand: Strand,
        pub has_soft_clip: bool,
        /// True if first or last `repeat_limit` bases repeat a mono- or di-nucleotide.
        pub has_repeat: bool,
        /// TAPS-aware mismatch count. Only computed when `indel_bases > 0`
        /// (expensive walk); zero otherwise.
        pub taps_aware_mismatches: u32,
    }

    /// Push-time record filter + per-record extras provider for seqair.
    #[derive(Clone, Default)]
    pub struct RastairRecordFilter {
        pub read_flags: ReadFlags,
        pub unpaired: bool,
        pub tag_requirement: TagRequirement,
        pub guess_orientation: bool,
        pub reference: Option<ReferenceWindow>,
        /// Repeat limit for `has_repeat` (matches `PileupMappingParams::indel_repeat_limit`).
        pub repeat_limit: usize,
    }

    impl CustomizeRecordStore for RastairRecordFilter {
        type Extra = RastairReadExtras;

        fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
            self.read_flags.filter_flags(fields.flags.raw(), self.unpaired)
        }

        fn filter(&mut self, rec: &SlimRecord, store: &RecordStore<Self::Extra>) -> bool {
            match &self.tag_requirement {
                TagRequirement::All => true,
                TagRequirement::AllOf(filters) => {
                    // Parity with the htslib path: a record whose aux block can't
                    // be read fails a tag *requirement* (fail closed, not open).
                    let Ok(aux) = rec.aux(store) else { return false };
                    filters
                        .iter()
                        .all(|f| aux.get::<&[u8]>(f.tag()).is_ok_and(|v| v == f.value_bytes()))
                }
            }
        }

        fn compute(
            &mut self,
            rec: &SlimRecord,
            store: &RecordStore<RastairReadExtras>,
        ) -> RastairReadExtras {
            let from_flags = strand_from_flags(rec.flags);

            let strand = if self.guess_orientation {
                if let Some(reference) = self.reference.as_ref() {
                    infer_strand_from_motifs(rec, store, reference).unwrap_or(from_flags)
                } else {
                    from_flags
                }
            } else {
                from_flags
            };

            let has_soft_clip = rec
                .cigar(store)
                .map(|ops| ops.iter().any(|op| op.op_type() == CigarOpType::SoftClip))
                .unwrap_or(false);

            let repeat_limit = self.repeat_limit;
            let has_repeat = if repeat_limit > 0 {
                rec.seq(store)
                    .map(|seq| {
                        has_repeat_seq(seq, 1, repeat_limit) || has_repeat_seq(seq, 2, repeat_limit)
                    })
                    .unwrap_or(false)
            } else {
                false
            };

            let has_indels = rec.indel_bases > 0;
            let taps_aware_mismatches = if has_indels {
                if let Some(reference) = self.reference.as_ref() {
                    count_taps_aware_mismatches(rec, store, reference, strand)
                } else {
                    0
                }
            } else {
                0
            };

            RastairReadExtras { strand, has_soft_clip, has_repeat, taps_aware_mismatches }
        }
    }

    /// Walk aligned pairs counting non-TAPS mismatches (same logic as `from_hts.rs`).
    fn count_taps_aware_mismatches<U>(
        rec: &SlimRecord,
        store: &RecordStore<U>,
        reference: &ReferenceWindow,
        strand: Strand,
    ) -> u32 {
        let Ok(pairs) = rec.aligned_pairs_with_read(store) else { return 0 };
        let mut count = 0u32;

        for matched in pairs.matches_only() {
            let rpos_u32 = u32::try_from(matched.rpos.as_u64()).unwrap_or(u32::MAX);
            let Some(ref_base) = reference.base_at(rpos_u32) else { continue };
            let observed = matched.query;

            if observed == Base::Unknown || ref_base == Base::Unknown || observed == ref_base {
                continue;
            }

            let is_taps_signal = match strand {
                Strand::OT => observed == Base::T && ref_base == Base::C,
                Strand::OB => observed == Base::A && ref_base == Base::G,
                Strand::Unknown => {
                    (observed == Base::T && ref_base == Base::C)
                        || (observed == Base::A && ref_base == Base::G)
                }
            };

            if !is_taps_signal {
                count += 1;
            }
        }
        count
    }

    /// Infer OT/OB strand from mismatch-motif evidence (same logic as `from_hts.rs`).
    fn infer_strand_from_motifs<U>(
        rec: &SlimRecord,
        store: &RecordStore<U>,
        reference: &ReferenceWindow,
    ) -> Option<Strand> {
        let pairs = rec.aligned_pairs_with_read(store).ok()?;
        let mut tg: u32 = 0;
        let mut ca: u32 = 0;

        for matched in pairs.matches_only() {
            let rpos_u32 = u32::try_from(matched.rpos.as_u64()).ok()?;
            let ref_base = reference.base_at(rpos_u32)?;
            let observed = matched.query;

            if ref_base == Base::Unknown || observed == Base::Unknown || ref_base == observed {
                continue;
            }

            count_motifs_at_with_store(rec, store, matched.qpos as usize, &mut tg, &mut ca);
        }

        match tg.cmp(&ca) {
            std::cmp::Ordering::Greater => Some(Strand::OT),
            std::cmp::Ordering::Less => Some(Strand::OB),
            std::cmp::Ordering::Equal => {
                // Tie: pseudo-random from qname+pos+flags (same approach as from_hts.rs)
                let name = rec.qname(store).unwrap_or(&[]);
                Some(pseudo_random_strand_from_bytes(name, rec.pos.as_u64(), rec.flags.raw()))
            }
        }
    }

    fn count_motifs_at_with_store<U>(
        rec: &SlimRecord,
        store: &RecordStore<U>,
        idx: usize,
        tg: &mut u32,
        ca: &mut u32,
    ) {
        let Ok(seq) = rec.seq(store) else { return };
        let pair_at = |first: usize, second: usize| -> Option<[Base; 2]> {
            Some([*seq.get(first)?, *seq.get(second)?])
        };
        let mut bump = |motif: [Base; 2]| match motif {
            [Base::T, Base::G] => *tg = tg.saturating_add(1),
            [Base::C, Base::A] => *ca = ca.saturating_add(1),
            _ => {}
        };
        if let Some(next) = idx.checked_add(1)
            && let Some(motif) = pair_at(idx, next)
        {
            bump(motif);
        }
        if let Some(prev) = idx.checked_sub(1)
            && let Some(motif) = pair_at(prev, idx)
        {
            bump(motif);
        }
    }

    fn pseudo_random_strand_from_bytes(name: &[u8], pos: u64, flags: u16) -> Strand {
        use rustc_hash::FxHasher;
        use std::hash::{Hash, Hasher};
        let mut hasher = FxHasher::default();
        name.hash(&mut hasher);
        pos.hash(&mut hasher);
        flags.hash(&mut hasher);
        if hasher.finish() & 1 == 0 { Strand::OT } else { Strand::OB }
    }

    /// Check if first or last `cutoff` bases repeat a pattern of length `n`.
    pub(crate) fn has_repeat_seq(seq: &[Base], n: usize, cutoff: usize) -> bool {
        let len = seq.len();
        if len < cutoff || n == 0 || cutoff < n {
            return false;
        }
        let start_pattern: SmallVec<Base, 4> =
            seq.get(..n).map(|s| s.iter().copied().collect()).unwrap_or_default();
        if start_pattern.len() == n {
            let start_repeat = (n..cutoff).all(|i| {
                seq.get(i).map_or(false, |&b| start_pattern.get(i % n).map_or(false, |&p| b == p))
            });
            if start_repeat {
                return true;
            }
        }
        let end_start = len.saturating_sub(n);
        let end_pattern: SmallVec<Base, 4> =
            seq.get(end_start..len).map(|s| s.iter().copied().collect()).unwrap_or_default();
        if end_pattern.len() == n {
            let check_start = len.saturating_sub(cutoff);
            (check_start..end_start).all(|i| {
                seq.get(i).map_or(false, |&b| {
                    let offset = (i - check_start) % n;
                    end_pattern.get(offset % n).map_or(false, |&p| b == p)
                })
            })
        } else {
            false
        }
    }

    /// Newtype wrapping `seqair::Readers<RastairRecordFilter>` that exposes
    /// the same `segments()` / `segment()` interface as the htslib `Readers`,
    /// so the `call` path can switch between backends with minimal call-site
    /// changes.
    pub struct SeqairReaders {
        pub(super) inner: seqair::Readers<RastairRecordFilter>,
        regions: Option<CliRegionInput>,
    }

    impl SeqairReaders {
        pub fn new(
            inner: seqair::Readers<RastairRecordFilter>,
            regions: Option<CliRegionInput>,
        ) -> Self {
            Self { inner, regions }
        }

        /// Produces the same `ChunkRegion` iterator as the htslib `Readers::segments()`.
        #[instrument(level = "debug", skip_all)]
        pub fn segments(
            &self,
            segment_max_length: u64,
            segment_overlap: u64,
        ) -> Result<impl Iterator<Item = ChunkRegion> + use<>> {
            let header = self.inner.header();
            let mut full_regions = if let Some(input) = &self.regions {
                input
                    .regions()
                    .iter()
                    .map(|r| get_selected_region(r, header))
                    .collect::<Result<Vec<_>>>()?
            } else {
                debug!("fetching all regions");
                get_full_regions(header)
            };
            ensure!(!full_regions.is_empty(), "No regions found");

            // Emit records in coordinate order regardless of the order regions were
            // given on the CLI: the VCF index builder requires tids (and positions
            // within a tid) to be monotonically non-decreasing. `--region chr7 chr1`
            // is reordered to chr1, chr7 here so the output stays sorted and indexable.
            full_regions.sort_by_key(|region| (header.tid(&region.contig), region.start));

            let initial_start =
                full_regions.first().wrap_err("No regions found").map(|r| r.start)?;
            let chunked = ChunkedRegions {
                full_regions,
                current_region_idx: 0,
                current_start: initial_start,
                max_length: segment_max_length,
                overlap: segment_overlap,
            };
            Ok(chunked)
        }

        /// Fetches the FASTA sequence for `region`, returning a rastair `Segment`.
        /// Mirrors the htslib `Readers::segment()` signature.
        #[instrument(level = "debug", skip_all)]
        pub fn segment(&mut self, region: &ChunkRegion, overfetch: u64) -> Result<Segment> {
            let last_position_to_fetch =
                region.end.wrapping_add(overfetch).min(region.last_position);
            let start = Pos0::try_from(region.start)
                .wrap_err_with(|| format!("region start {} is out of range", region.start))?;
            let end = Pos0::try_from(last_position_to_fetch).wrap_err_with(|| {
                format!("region end {} is out of range", last_position_to_fetch)
            })?;

            let bases =
                self.inner.fetch_base_seq(&region.contig, start, end).wrap_err_with(|| {
                    format!("Failed to get region {} from FASTA file", region.region)
                })?;

            // Base is #[repr(u8)] with ASCII discriminants, so `*b as u8` is safe.
            let sequence: Vec<u8> = bases.iter().map(|b| *b as u8).collect();

            Ok(Segment {
                range: region.clone(),
                sequence,
                overlap_start: region.overlap_start,
                overlap_end: region.overlap_end,
            })
        }

        pub fn header(&self) -> &seqair::bam::BamHeader {
            self.inner.header()
        }

        pub(crate) fn inner_mut(&mut self) -> &mut seqair::Readers<RastairRecordFilter> {
            &mut self.inner
        }
    }

    fn get_selected_region(
        region: &RegionString,
        header: &seqair::bam::BamHeader,
    ) -> Result<SelectedRegion> {
        let chromosome = region.chromosome.as_str();
        let start = region.start.map(|p| p.as_u64()).unwrap_or(1);

        let contig_info = header.resolve_contig(chromosome).wrap_err_with(|| {
            format!("Failed to resolve chromosome {} from header", region.chromosome)
        })?;
        let last_position = contig_info.len;
        let end = region.end.map(|p| p.as_u64()).unwrap_or(last_position);

        ensure!(
            start <= last_position,
            "Specified start position {start} past the end of chromosome {chromosome}"
        );
        ensure!(
            end <= last_position,
            "Specified end position {end} past the end of chromosome {chromosome}"
        );

        Ok(SelectedRegion::UserDefinedSubset {
            region: Region { contig: region.chromosome.clone(), start, end },
            last_position,
        })
    }

    fn get_full_regions(header: &seqair::bam::BamHeader) -> Vec<SelectedRegion> {
        header
            .targets()
            .map(|target| {
                SelectedRegion::EntireContig(Region {
                    contig: target.name,
                    start: 1,
                    end: target.length,
                })
            })
            .collect()
    }
}

#[derive(Debug, clap::Args, Clone)]
pub struct ReaderParams {
    /// Path to sorted and indexed BAM or CRAM file
    #[arg(value_parser=value_parser!(ClioPath).exists().is_file(), value_hint=clap::ValueHint::FilePath)]
    #[arg(help_heading = cli::sections::INPUT)]
    pub bam_file: ClioPath,

    /// Path to sorted and indexed (via samtools faidx) FASTA file. Can be bgzip
    /// compressed, but requires both a gzi index and a fai index
    #[arg(short='r', long, value_parser=value_parser!(ClioPath).exists().is_file(), value_hint=clap::ValueHint::FilePath)]
    #[arg(help_heading = cli::sections::INPUT)]
    pub fasta_file: ClioPath,

    /// Restrict processing to specific genomic regions.
    ///
    /// Accepts either space-separated region strings or a single BED file:
    /// - Region strings: `chr`, `chr:start`, `chr:start-end` (1-based inclusive)
    /// - Multiple regions separated by whitespace: `"chr1 chr2:100-200"`
    /// - BED files with `@` prefix: `@targets.bed`
    #[arg(short = 'l', long = "region", value_parser = clap::value_parser!(CliRegionInput))]
    #[arg(help_heading = cli::sections::INPUT)]
    pub regions: Option<CliRegionInput>,
}

impl ReaderParams {
    pub fn readers(&self) -> Result<Readers> {
        let fasta = open_fasta(&self.fasta_file)?;
        let bam_path = self.bam_file.path();
        let mut bam = bam::IndexedReader::from_path(bam_path)
            .with_suggestion(|| {
                format!(
                    "Ensure the BAM/CRAM file is sorted and indexed with \
                    `samtools sort {bam_path:?}` and `samtools index {bam_path:?}`, respectively."
                )
            })
            .note("If you have a .bai/.crai file, ensure it is in the same directory as the BAM/CRAM file.")?;
        bam.set_reference(self.fasta_file.path())
            .wrap_err("Failed to set FASTA reference for BAM reader")
            .note("Rastair itself already opened the FASTA file successfully, this error is from the BAM reader implementation")?;

        Ok(Readers { fasta, bam, params: self.clone() })
    }

    #[cfg(feature = "experimental-seqair")]
    pub fn open_seqair(&self) -> Result<SeqairReaders> {
        use color_eyre::eyre::WrapErr as _;
        let inner = seqair::Readers::open_customized(
            self.bam_file.path(),
            self.fasta_file.path(),
            seqair_readers::RastairRecordFilter::default(),
        )
        .wrap_err_with(|| {
            format!("Failed to open alignment file {}", self.bam_file.path().display())
        })?;
        Ok(SeqairReaders::new(inner, self.regions.clone()))
    }

    #[cfg(feature = "experimental-seqair")]
    pub fn pileup_readers(&self) -> Result<SeqairReaders> {
        self.open_seqair()
    }

    #[cfg(not(feature = "experimental-seqair"))]
    pub fn pileup_readers(&self) -> Result<Readers> {
        self.readers()
    }
}

pub struct Readers {
    pub fasta: FastaReader,
    pub bam: bam::IndexedReader,
    params: ReaderParams,
}

impl Readers {
    /// Calculate segments based on configuration parameters
    #[instrument(level = "debug", skip_all)]
    pub fn segments(
        &self,
        segment_max_length: u64,
        segment_overlap: u64,
    ) -> Result<impl Iterator<Item = ChunkRegion> + use<>> {
        let mut full_regions = if let Some(input) = &self.params.regions {
            input
                .regions()
                .iter()
                .map(|r| {
                    get_selected_region(r, self.bam.header())
                        .wrap_err("Failed to get selected region from BAM file")
                })
                .collect::<Result<Vec<_>>>()?
        } else {
            debug!("fetching all regions");
            get_full_regions(self.bam.header())
                .wrap_err("Failed to get all regions from BAM file")?
        };

        // Emit records in coordinate order regardless of CLI region order: the VCF
        // index builder requires tids (and positions within a tid) to be monotonically
        // non-decreasing, so `--region chr7 chr1` is reordered to chr1, chr7 here.
        let header = self.bam.header();
        full_regions.sort_by_key(|region| (header.tid(region.contig.as_bytes()), region.start));

        let initial_start = full_regions.first().wrap_err("No regions found")?.start;
        let chunked = ChunkedRegions {
            full_regions,
            current_region_idx: 0,
            current_start: initial_start,
            max_length: segment_max_length,
            overlap: segment_overlap,
        };

        Ok(chunked)
    }

    /// Fetch a segment from the FASTA file, with optional overfetching
    #[instrument(level = "debug", skip_all)]
    pub fn segment(&mut self, region: &ChunkRegion, overfetch: u64) -> Result<Segment> {
        let last_position_to_fetch = region.end.wrapping_add(overfetch).min(region.last_position);

        // Calculate exact capacity needed to avoid reallocations
        let len = usize::try_from(last_position_to_fetch.wrapping_sub(region.start))
            .wrap_err("Failed to convert segment length to usize")?;

        trace!(?region, len, "fetching segment");
        let seq = self
            .fasta
            .fetch_seq(&region.contig, region.start, last_position_to_fetch)
            .wrap_err_with(|| format!("Failed to get region {} from FASTA file", region.region))?;

        Ok(Segment {
            range: region.clone(),
            sequence: seq,
            overlap_start: region.overlap_start,
            overlap_end: region.overlap_end,
        })
    }
}

#[instrument(level = "debug", skip(bam_header))]
fn get_selected_region(region: &RegionString, bam_header: &HeaderView) -> Result<SelectedRegion> {
    let chromosome = region.chromosome.as_str();
    // If no start position is specified, default to beginning of chromosome
    let start = region.start.map(to_u64).unwrap_or(1);

    let target_id = bam_header
        .tid(region.chromosome.as_bytes())
        .wrap_err_with(|| {
            format!("Failed to fetch target ID for chromosome {} from header", region.chromosome)
        })
        .with_note(|| {
            format!(
                "This usually means the specified chromosome {} is not in the input BAM file",
                region.chromosome
            )
        })?;
    let last_position =
        bam_header.target_len(target_id).wrap_err("Failed to fetch header length")?;
    // If no end specified, use chromosome length from BAM header
    let end = region.end.map(to_u64).unwrap_or(last_position);

    ensure!(
        start <= last_position,
        "Specified start position {end} past the end of chromosome {chromosome}"
    );
    ensure!(
        end <= last_position,
        "Specified end position {end} past the end of chromosome {chromosome}"
    );

    // Since the user specified this region, we're only returning that one
    Ok(SelectedRegion::UserDefinedSubset {
        region: Region { contig: region.chromosome.clone(), start, end },
        last_position,
    })
}

fn to_u64(value: seqair_types::Pos1) -> u64 {
    value.as_u64()
}

#[instrument(level = "debug", skip(header))]
fn get_full_regions(header: &bam::HeaderView) -> Result<Vec<SelectedRegion>> {
    header
        .target_names()
        .iter()
        .enumerate()
        .filter(|(_, name)| !name.is_empty())
        .map(|(tid, name)| -> Result<SelectedRegion> {
            let chr = SmolStr::new(
                std::str::from_utf8(name).wrap_err("BAM target name not valid UTF-8")
                    .note("This is against the BAM specification, please check with the tool that created this file")?,
            );
            let length = header
                .target_len(u32::try_from(tid).wrap_err("Failed to get a target ID that was part of the BAM header")
                    .note("The BAM header might be corrupt")?)
                .wrap_err("Failed to get target length")?;

            Ok(SelectedRegion::EntireContig(Region {
                contig: chr,
                start: 1, // 1-based coordinates
                end: length,
            }))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests").join("data")
    }

    fn get_test_bam() -> ClioPath {
        ClioPath::new(test_data_dir().join("test.bam")).expect("test bam path should be valid")
    }

    fn get_test_fasta() -> ClioPath {
        ClioPath::new(test_data_dir().join("test.fasta.gz"))
            .expect("test fasta path should be valid")
    }

    #[test]
    fn test_get_selected_region_variations() -> Result<()> {
        let params =
            ReaderParams { bam_file: get_test_bam(), fasta_file: get_test_fasta(), regions: None };

        let readers = params.readers()?;
        let header = readers.bam.header();

        // Test chromosome-only region
        let region_chr_only: RegionString = "chr19".parse().unwrap();
        let full_region = get_selected_region(&region_chr_only, header)?;

        assert_eq!(full_region.contig, "chr19");
        assert_eq!(full_region.start, 1); // Should default to 1

        // The end should be the chromosome length from the header
        let chr19_tid = header.tid(b"chr19").unwrap();
        let chr19_len = header.target_len(chr19_tid).unwrap();
        assert_eq!(full_region.end, chr19_len);

        // Test chromosome with start but no end
        let region_with_start: RegionString = "chr19:100".parse().unwrap();
        let full_region = get_selected_region(&region_with_start, header)?;

        assert_eq!(full_region.contig, "chr19");
        assert_eq!(full_region.start, 100);
        assert_eq!(full_region.end, chr19_len); // Should default to chromosome length

        Ok(())
    }

    #[test]
    fn test_get_selected_region_errors() -> Result<()> {
        let params =
            ReaderParams { bam_file: get_test_bam(), fasta_file: get_test_fasta(), regions: None };

        let readers = params.readers()?;
        let header = readers.bam.header();

        // Test non-existent chromosome
        let region_invalid_chr: RegionString = "nonexistent".parse().unwrap();
        let result = get_selected_region(&region_invalid_chr, header);
        assert!(result.is_err());

        // Get valid chromosome length
        let chr19_tid = header.tid(b"chr19").unwrap();
        let chr19_len = header.target_len(chr19_tid).unwrap();

        // Test start position beyond chromosome length
        let invalid_start = chr19_len + 100;
        let region_invalid_start: RegionString = format!("chr19:{invalid_start}").parse().unwrap();
        let result = get_selected_region(&region_invalid_start, header);
        assert!(result.is_err());

        // Test end position beyond chromosome length
        let invalid_end = chr19_len + 100;
        let region_invalid_end: RegionString = format!("chr19:100-{invalid_end}").parse().unwrap();
        let result = get_selected_region(&region_invalid_end, header);
        assert!(result.is_err());

        Ok(())
    }

    #[test]
    fn segments_sort_regions_into_tid_order() -> Result<()> {
        use crate::utils::regions::CliRegionInput;

        // Given out of tid order on the CLI: `bacteriophage_lambda_CpG` is tid 2,
        // `chr19` is tid 0. The emitted chunks must come out in tid order so the
        // VCF index builder accepts the (monotonically non-decreasing tid) stream.
        let regions: CliRegionInput = "bacteriophage_lambda_CpG chr19".parse()?;
        let params = ReaderParams {
            bam_file: get_test_bam(),
            fasta_file: get_test_fasta(),
            regions: Some(regions),
        };

        let readers = params.readers()?;
        let chunks: Vec<ChunkRegion> = readers.segments(u64::MAX, 0)?.collect();

        let mut contig_order: Vec<&str> = Vec::new();
        for chunk in &chunks {
            if contig_order.last() != Some(&chunk.contig.as_str()) {
                contig_order.push(chunk.contig.as_str());
            }
        }

        assert_eq!(contig_order, ["chr19", "bacteriophage_lambda_CpG"]);

        Ok(())
    }

    #[test]
    fn test_get_full_regions() -> Result<()> {
        let params =
            ReaderParams { bam_file: get_test_bam(), fasta_file: get_test_fasta(), regions: None };

        let readers = params.readers()?;
        let header = readers.bam.header();

        let full_regions = get_full_regions(header)?;

        // Should have at least one region
        assert!(!full_regions.is_empty());

        // Verify that chromosome names match what's in the BAM header
        for (i, region) in full_regions.iter().enumerate() {
            let target_name = std::str::from_utf8(header.target_names()[i]).unwrap();
            assert_eq!(region.contig, target_name);

            // Start should be 1 (1-based)
            assert_eq!(region.start, 1);

            // End should match the chromosome length
            let tid = u32::try_from(i).unwrap();
            let chr_len = header.target_len(tid).unwrap();
            assert_eq!(region.end, chr_len);
        }

        Ok(())
    }
}
