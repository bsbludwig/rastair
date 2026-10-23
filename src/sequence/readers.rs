//! BAM/CRAM and FASTA readers: [`SeqairReaders`] for `call` and `ml train`, the
//! htslib [`Readers`] for `per-read` and the BAM rewrite.
//!
//! Each record's OT/OB strand is assigned once from its SAM flags, when seqair
//! loads it (`RastairRecordFilter::compute`), and kept in [`RastairReadExtras`].

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

#[cfg(test)]
pub(crate) use seqair_readers::RastairRecordFilter;
pub use seqair_readers::{RastairReadExtras, SeqairReaders};

mod seqair_readers {
    use crate::{
        call::{
            require_tags::TagRequirement,
            variant_calling::{ReadFlags, ReadMaskParams},
        },
        sequence::{
            ChunkRegion, Region, SelectedRegion, chunked::ChunkedRegions, segementation::Segment,
        },
        utils::{CliRegionInput, RegionString},
    };
    use color_eyre::eyre::{ContextCompat as _, Result, WrapErr as _, ensure};
    use seqair::bam::aux::Aux;
    use seqair::bam::record_store::{
        CustomizeRecordStore, FilterRawFields, RecordStore, SlimRecord,
    };
    use seqair_types::{Pos0, Strand};
    use std::ops::Range;
    use tracing::{debug, instrument};

    /// Per-record data computed once at fetch time.
    #[derive(Clone)]
    pub struct RastairReadExtras {
        pub strand: Strand,
        /// Read positions that survive read-end masking. Resolved here from
        /// (strand, reverse, `seq_len`) so the per-column check is one range
        /// comparison — masking is evaluated at every base of every read, and
        /// none of its inputs change between columns.
        pub mask: Range<u32>,
    }

    /// Push-time record filter + per-record extras provider for seqair.
    #[derive(Clone, Default)]
    pub struct RastairRecordFilter {
        pub read_flags: ReadFlags,
        pub unpaired: bool,
        pub min_mapq: u8,
        pub tag_requirement: TagRequirement,
        pub read_masking: ReadMaskParams,
        /// Mapped records the reader fetched, before any of these filters:
        /// the same count a `.bai` keeps, for the phasing depth guard.
        pub fetched_mapped: u64,
    }

    impl CustomizeRecordStore for RastairRecordFilter {
        type Extra = RastairReadExtras;

        /// Every per-read filter runs here, before the depth cap counts the
        /// read, so a deep pile of MAPQ 0 reads cannot fill the cap and evict
        /// the usable ones.
        fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
            if !fields.flags.is_unmapped() {
                self.fetched_mapped = self.fetched_mapped.saturating_add(1);
            }
            self.read_flags.filter_flags(fields.flags.raw(), self.unpaired)
                && fields.mapq >= self.min_mapq
                && self.tag_requirement.allows(Some(Aux::new(fields.aux_bytes)))
        }

        /// The mask needs the strand, which `compute` resolves, so this read
        /// has already taken a depth-cap slot; only an extreme read-end mask
        /// empties a whole read. A read whose extras cannot be read is
        /// rejected too.
        fn filter(&mut self, rec: &SlimRecord, store: &RecordStore<Self::Extra>) -> bool {
            rec.extra(store).is_ok_and(|extra| !extra.mask.is_empty())
        }

        fn compute(
            &mut self,
            rec: &SlimRecord,
            _store: &RecordStore<RastairReadExtras>,
        ) -> RastairReadExtras {
            let strand = Strand::from(rec.flags);

            let mask = self
                .read_masking
                .keep_range(strand, rec.flags.is_reverse(), rec.seq_len)
                .unwrap_or(0..0);

            RastairReadExtras { strand, mask }
        }
    }

    /// Newtype wrapping `seqair::Readers<RastairRecordFilter>` with the
    /// `segments()` / `segment()` interface the htslib `Readers` of `per-read`
    /// and the BAM rewrite also have.
    pub struct SeqairReaders {
        pub(super) inner: seqair::Readers<RastairRecordFilter>,
        regions: Option<CliRegionInput>,
        /// The whole file's mapped reads per covered base, from its index
        /// ([`mapped_reads_per_base`](super::mapped_reads_per_base)); set
        /// only when phasing asks for it.
        pub run_reads_per_base: Option<f64>,
    }

    impl SeqairReaders {
        pub fn new(
            inner: seqair::Readers<RastairRecordFilter>,
            regions: Option<CliRegionInput>,
        ) -> Self {
            Self { inner, regions, run_reads_per_base: None }
        }

        /// A handle for another thread that shares this one's parsed index and
        /// header, and for CRAM its cache of decoded slices, so neighbouring
        /// segments on different threads decode each slice once.
        pub fn fork(&self) -> Result<Self> {
            let inner = self.inner.fork().wrap_err("Failed to fork seqair readers")?;
            Ok(Self {
                inner,
                regions: self.regions.clone(),
                run_reads_per_base: self.run_reads_per_base,
            })
        }

        /// Produces the same `ChunkRegion` iterator as `Readers::segments()`.
        #[instrument(level = "debug", skip_all)]
        pub fn segments(
            &self,
            segment_max_length: u64,
            segment_overlap: u64,
        ) -> Result<impl Iterator<Item = ChunkRegion> + use<>> {
            let header = self.inner.header();
            let full_regions = if let Some(input) = &self.regions {
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

            // Drop BAM contigs the FASTA reference lacks (e.g. decoy/alt contigs);
            // explicitly requested regions instead stay a hard error.
            let fasta_contigs: std::collections::HashSet<seqair_types::SmolStr> =
                self.inner.fasta().index().sequence_names().into_iter().collect();
            let mut full_regions = crate::sequence::regions::retain_fasta_regions(
                full_regions,
                |contig| fasta_contigs.contains(contig),
                self.regions.is_some(),
            )?;

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
        #[instrument(level = "debug", skip_all)]
        pub fn segment(&mut self, region: &ChunkRegion, overfetch: u64) -> Result<Segment> {
            let last_position_to_fetch =
                region.end.wrapping_add(overfetch).min(region.last_position);
            let start = Pos0::try_from(region.start)
                .wrap_err_with(|| format!("region start {} is out of range", region.start))?;
            // `last_position_to_fetch` is an exclusive bound here, as it is for
            // `Readers::segment` (`FastaReader::fetch_seq` takes `[start, stop)`), while
            // seqair's span is closed on both ends — hence the `- 1`. A region
            // whose bound coincides with its start names no base at all, and a
            // closed span cannot say that, so it short-circuits to an empty
            // sequence rather than asking the reader for a reversed span.
            let sequence = match last_position_to_fetch.checked_sub(1) {
                Some(last) if last >= region.start => {
                    let last = Pos0::try_from(last).wrap_err_with(|| {
                        format!("region end {last_position_to_fetch} is out of range")
                    })?;
                    let bases = self
                        .inner
                        .fetch_base_seq(&region.contig, (start..=last).into())
                        .wrap_err_with(|| {
                            format!("Failed to get region {} from FASTA file", region.region)
                        })?;
                    // Base is #[repr(u8)] with ASCII discriminants, so `*b as u8` is safe.
                    bases.iter().map(|b| *b as u8).collect()
                }
                _ => Vec::new(),
            };

            Ok(Segment {
                range: std::sync::Arc::new(region.clone()),
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

        /// The run's mapped reads per base over those of the segment
        /// `get_pileups` loaded last, `segment_len` bases long; `None`
        /// without [`Self::run_reads_per_base`] or reads in the segment.
        pub fn depth_scale(&self, segment_len: u64) -> Option<f64> {
            let segment = self.inner.customize().fetched_mapped as f64 / segment_len as f64;
            let run = self.run_reads_per_base?;
            (segment.is_finite() && segment > 0.0).then(|| run / segment)
        }
    }

    fn get_selected_region(
        region: &RegionString,
        header: &seqair::bam::BamHeader,
    ) -> Result<SelectedRegion> {
        let chromosome = region.chromosome.as_str();
        let start = region.start.map(|p| p.to_zero_based().as_u64()).unwrap_or(0);

        let contig_info = header.resolve_contig(chromosome).wrap_err_with(|| {
            format!("Failed to resolve chromosome {} from header", region.chromosome)
        })?;
        let last_position = contig_info.len;
        // A closed 1-based end is the same number as a half-open 0-based one.
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
                    start: 0,
                    end: target.length,
                })
            })
            .collect()
    }
}

pub use index_depth::mapped_reads_per_base;

/// seqair parses the index but keeps neither its read counts nor its bins, so
/// this reads them itself.
mod index_depth {
    use color_eyre::eyre::{ContextCompat as _, Result, WrapErr as _, ensure};
    use noodles::csi::{BinningIndex as _, binning_index::ReferenceSequence as _};
    use std::{ops::Range, path::Path};

    /// What a BAI or CSI index says about one reference.
    #[derive(Debug, PartialEq, Eq)]
    struct ReferenceCounts {
        /// `None` when the index keeps no pseudo-bin for it.
        mapped: Option<u64>,
        /// Bases in the index's finest bins that hold a read: the part of a
        /// contig a sliced BAM left out is not counted.
        leaf_bases: u64,
    }

    /// The finest bins of a binning scheme, and the bases each spans.
    fn leaves(min_shift: u8, depth: u8) -> (Range<u64>, u64) {
        let first = ((1u64 << (3 * u32::from(depth))) - 1) / 7;
        (first..first + (1u64 << (3 * u32::from(depth))), 1u64 << min_shift)
    }

    /// Mapped reads per covered base over the whole alignment file, from its
    /// `.csi` or `.bai` (looked up as seqair does); `contig_lengths` are the
    /// header's, by tid, and cap a short contig's last bin. `None` without
    /// such an index (a CRAM), or when it holds no read counts.
    pub fn mapped_reads_per_base(alignment: &Path, contig_lengths: &[u64]) -> Result<Option<f64>> {
        let with_suffix = |suffix: &str| {
            let mut path = alignment.as_os_str().to_owned();
            path.push(suffix);
            std::path::PathBuf::from(path)
        };
        let references = if let Some(csi) = [with_suffix(".csi"), alignment.with_extension("csi")]
            .into_iter()
            .find(|path| path.exists())
        {
            let index = noodles::csi::fs::read(&csi)
                .wrap_err_with(|| format!("Failed to read {}", csi.display()))?;
            csi_counts(&index)
        } else if let Some(bai) = [with_suffix(".bai"), alignment.with_extension("bai")]
            .into_iter()
            .find(|path| path.exists())
        {
            let bytes = std::fs::read(&bai)
                .wrap_err_with(|| format!("Failed to read {}", bai.display()))?;
            bai_counts(&bytes).wrap_err_with(|| format!("Failed to parse {}", bai.display()))?
        } else {
            return Ok(None);
        };
        Ok(reads_per_base(&references, contig_lengths))
    }

    fn reads_per_base(references: &[ReferenceCounts], contig_lengths: &[u64]) -> Option<f64> {
        let mut mapped = None;
        let mut covered = 0;
        for (tid, reference) in references.iter().enumerate() {
            if let Some(count) = reference.mapped {
                *mapped.get_or_insert(0) += count;
            }
            let length = contig_lengths.get(tid).copied().unwrap_or(u64::MAX);
            covered += reference.leaf_bases.min(length);
        }
        let mapped = mapped.filter(|&mapped| mapped > 0 && covered > 0)?;
        Some(mapped as f64 / covered as f64)
    }

    fn csi_counts(index: &noodles::csi::Index) -> Vec<ReferenceCounts> {
        let (leaves, leaf_len) = leaves(index.min_shift(), index.depth());
        let references = index.reference_sequences().iter();
        references
            .map(|reference| ReferenceCounts {
                mapped: reference.metadata().map(|metadata| metadata.mapped_record_count()),
                leaf_bases: leaf_len
                    * reference.bins().keys().filter(|&&bin| leaves.contains(&(bin as u64))).count()
                        as u64,
            })
            .collect()
    }

    /// SAM spec §5.2: per reference, bins with their chunks, the pseudo-bin
    /// 37450 among them whose second chunk is (mapped, unmapped), then the
    /// linear index.
    fn bai_counts(bytes: &[u8]) -> Result<Vec<ReferenceCounts>> {
        const META_BIN: u64 = 37450;
        let mut index = LeBytes(bytes);
        ensure!(index.u32()? == u32::from_le_bytes(*b"BAI\x01"), "Not a BAI index");
        let (leaves, leaf_len) = leaves(14, 5);
        let mut references = Vec::new();
        for _ in 0..index.u32()? {
            let mut reference = ReferenceCounts { mapped: None, leaf_bases: 0 };
            for _ in 0..index.u32()? {
                let bin = u64::from(index.u32()?);
                let n_chunk = usize::try_from(index.u32()?)?;
                let chunks = index.take(16 * n_chunk)?;
                if bin == META_BIN {
                    let count = chunks.get(16..24).wrap_err("Pseudo-bin lacks its counts")?;
                    reference.mapped = Some(u64::from_le_bytes(count.try_into()?));
                } else if leaves.contains(&bin) {
                    reference.leaf_bases += leaf_len;
                }
            }
            let n_intv = usize::try_from(index.u32()?)?;
            index.take(8 * n_intv)?;
            references.push(reference);
        }
        Ok(references)
    }

    struct LeBytes<'a>(&'a [u8]);

    impl<'a> LeBytes<'a> {
        fn take(&mut self, n: usize) -> Result<&'a [u8]> {
            let (head, tail) = self.0.split_at_checked(n).wrap_err("Index ends early")?;
            self.0 = tail;
            Ok(head)
        }

        fn u32(&mut self) -> Result<u32> {
            Ok(u32::from_le_bytes(self.take(4)?.try_into()?))
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn bai_counts_match_samtools_idxstats() -> Result<()> {
            let bytes = std::fs::read("tests/data/test.bam.bai")?;
            let counts = |mapped, leaves: u64| ReferenceCounts {
                mapped: Some(mapped),
                leaf_bases: leaves * 16_384,
            };
            // Mapped from `samtools idxstats tests/data/test.bam`; the leaves
            // agree with a CSI that `samtools index -c` built.
            assert_eq!(
                bai_counts(&bytes)?,
                [counts(12_297, 3), counts(2_146, 1), counts(10_500, 3)]
            );
            Ok(())
        }

        #[test]
        fn a_contig_shorter_than_its_bin_counts_its_length() -> Result<()> {
            let bytes = std::fs::read("tests/data/test.bam.bai")?;
            let lengths = [61_431_566, 2_018, 48_502];
            let covered = 3 * 16_384 + 2_018 + 48_502;
            let rate = reads_per_base(&bai_counts(&bytes)?, &lengths);
            assert_eq!(rate, Some(24_943.0 / f64::from(covered)));
            Ok(())
        }
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

    pub fn pileup_readers(&self) -> Result<SeqairReaders> {
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
        let full_regions = if let Some(input) = &self.params.regions {
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

        // Drop BAM contigs the FASTA reference lacks (e.g. decoy/alt contigs);
        // explicitly requested regions instead stay a hard error.
        let fasta_contigs = self.fasta.sequence_names()?;
        let mut full_regions = crate::sequence::regions::retain_fasta_regions(
            full_regions,
            |contig| fasta_contigs.contains(contig),
            self.params.regions.is_some(),
        )?;

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
            range: std::sync::Arc::new(region.clone()),
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
    let start = region.start.map(|p| p.to_zero_based().as_u64()).unwrap_or(0);

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
    // A closed 1-based end is the same number as a half-open 0-based one.
    let end = region.end.map(|p| p.as_u64()).unwrap_or(last_position);

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
                start: 0,
                end: length,
            }))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    /// `--require-tags` is judged on the raw aux bytes, before the depth cap
    /// counts the read, so a read without the tag never takes a slot.
    #[test]
    fn record_filter_rejects_a_missing_tag_before_loading() {
        use crate::call::require_tags::RequireTagsParams;
        use seqair::bam::aux_data::AuxData;
        use seqair::bam::cigar::{CigarOp, CigarOpType};
        use seqair::bam::{BamHeader, MateFields, record_store::RecordStore};
        use seqair::reader::ResolveTid as _;
        use seqair_types::{BamFlags, Base, Pos0};

        let mut filter = RastairRecordFilter {
            min_mapq: 1,
            tag_requirement: RequireTagsParams { require_tags: vec!["RG=keep".parse().unwrap()] }
                .filter(),
            ..Default::default()
        };
        let tid = 0u32.resolve_tid(&BamHeader::synthetic(1, 1_000)).unwrap();
        let mut store = RecordStore::<RastairReadExtras>::new();
        for (name, group) in [(&b"tagged"[..], "keep"), (b"other", "drop")] {
            let mut aux = AuxData::new();
            aux.set_string(*b"RG", group.as_bytes());
            store
                .push_fields(
                    Pos0::new(10).unwrap(),
                    Pos0::new(17).unwrap(),
                    BamFlags::from(99u16),
                    60,
                    8,
                    0,
                    name,
                    &[CigarOp::new(CigarOpType::Match, 8)],
                    &[Base::A; 8],
                    &[40u8; 8],
                    aux.as_bytes(),
                    Some(tid),
                    MateFields { tid: Some(tid), pos: Pos0::new(500), template_len: None },
                    &mut filter,
                )
                .unwrap();
        }
        assert_eq!(store.len(), 1, "only the read carrying RG=keep is loaded");
    }

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

    /// Writes a plain FASTA holding only `contigs` (dummy sequence) with a `.fai`
    /// index, returning the temp dir (kept alive by the caller) and a `ClioPath`
    /// to the FASTA. Used to build BAM/FASTA reference mismatches on the fly.
    fn fasta_with_contigs(contigs: &[&str]) -> (tempfile::TempDir, ClioPath) {
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().expect("tempdir");
        let path = dir.path().join("ref.fa");
        let mut file = std::fs::File::create(&path).expect("create fasta");
        for contig in contigs {
            writeln!(file, ">{contig}\nACGTACGTAC").expect("write fasta");
        }
        drop(file);
        rust_htslib::faidx::build(&path).expect("build fai");
        let clio = ClioPath::new(path).expect("clio path");
        (dir, clio)
    }

    mod htslib_fasta_filtering {
        use super::*;

        #[test]
        fn segments_skip_bam_contig_absent_from_fasta() -> Result<()> {
            // FASTA omits `bacteriophage_lambda_CpG`, which the BAM header carries.
            let (_dir, fasta) = fasta_with_contigs(&["chr19", "2kb_3_Unmodified"]);
            let params =
                ReaderParams { bam_file: get_test_bam(), fasta_file: fasta, regions: None };

            let readers = params.readers()?;
            let chunks: Vec<ChunkRegion> = readers.segments(u64::MAX, 0)?.collect();

            let contigs: std::collections::HashSet<&str> =
                chunks.iter().map(|c| c.contig.as_str()).collect();
            assert!(contigs.contains("chr19"));
            assert!(contigs.contains("2kb_3_Unmodified"));
            assert!(
                !contigs.contains("bacteriophage_lambda_CpG"),
                "contig absent from the FASTA must be skipped, got {contigs:?}"
            );
            Ok(())
        }

        #[test]
        fn segments_bail_when_no_bam_contig_in_fasta() -> Result<()> {
            let (_dir, fasta) = fasta_with_contigs(&["unrelated_contig"]);
            let params =
                ReaderParams { bam_file: get_test_bam(), fasta_file: fasta, regions: None };

            let readers = params.readers()?;
            let result = readers.segments(u64::MAX, 0).map(Iterator::count);
            assert!(result.is_err(), "expected an error when no BAM contig is in the FASTA");
            Ok(())
        }

        #[test]
        fn segments_explicit_region_absent_from_fasta_errors() -> Result<()> {
            use crate::utils::regions::CliRegionInput;

            let (_dir, fasta) = fasta_with_contigs(&["chr19"]);
            let regions: CliRegionInput = "bacteriophage_lambda_CpG".parse()?;
            let params = ReaderParams {
                bam_file: get_test_bam(),
                fasta_file: fasta,
                regions: Some(regions),
            };

            let readers = params.readers()?;
            let result = readers.segments(u64::MAX, 0).map(Iterator::count);
            assert!(result.is_err(), "explicit region absent from FASTA must be a hard error");
            Ok(())
        }
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
        assert_eq!(full_region.start, 0); // Should default to the first base

        // The end should be the chromosome length from the header
        let chr19_tid = header.tid(b"chr19").unwrap();
        let chr19_len = header.target_len(chr19_tid).unwrap();
        assert_eq!(full_region.end, chr19_len);

        // Test chromosome with start but no end
        let region_with_start: RegionString = "chr19:100".parse().unwrap();
        let full_region = get_selected_region(&region_with_start, header)?;

        assert_eq!(full_region.contig, "chr19");
        assert_eq!(full_region.start, 99); // 1-based 100, 0-based 99
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

            // The contig's first base, 0-based
            assert_eq!(region.start, 0);

            // End should match the chromosome length
            let tid = u32::try_from(i).unwrap();
            let chr_len = header.target_len(tid).unwrap();
            assert_eq!(region.end, chr_len);
        }

        Ok(())
    }

    mod seqair_fasta_filtering {
        use super::*;

        #[test]
        fn segments_skip_bam_contig_absent_from_fasta() -> Result<()> {
            // FASTA omits `bacteriophage_lambda_CpG`, which the BAM header carries.
            let (_dir, fasta) = fasta_with_contigs(&["chr19", "2kb_3_Unmodified"]);
            let params =
                ReaderParams { bam_file: get_test_bam(), fasta_file: fasta, regions: None };

            let readers = params.pileup_readers()?;
            let chunks: Vec<ChunkRegion> = readers.segments(u64::MAX, 0)?.collect();

            let contigs: std::collections::HashSet<&str> =
                chunks.iter().map(|c| c.contig.as_str()).collect();
            assert!(contigs.contains("chr19"));
            assert!(contigs.contains("2kb_3_Unmodified"));
            assert!(
                !contigs.contains("bacteriophage_lambda_CpG"),
                "contig absent from the FASTA must be skipped, got {contigs:?}"
            );
            Ok(())
        }

        #[test]
        fn segments_bail_when_no_bam_contig_in_fasta() -> Result<()> {
            // FASTA shares no contig names with the BAM header (wholesale mismatch).
            let (_dir, fasta) = fasta_with_contigs(&["unrelated_contig"]);
            let params =
                ReaderParams { bam_file: get_test_bam(), fasta_file: fasta, regions: None };

            let readers = params.pileup_readers()?;
            let result = readers.segments(u64::MAX, 0).map(Iterator::count);
            assert!(result.is_err(), "expected an error when no BAM contig is in the FASTA");
            Ok(())
        }

        #[test]
        fn segments_explicit_region_absent_from_fasta_errors() -> Result<()> {
            use crate::utils::regions::CliRegionInput;

            let (_dir, fasta) = fasta_with_contigs(&["chr19"]);
            // The contig exists in the BAM header but not in this FASTA: an
            // explicitly requested region must stay a hard error.
            let regions: CliRegionInput = "bacteriophage_lambda_CpG".parse()?;
            let params = ReaderParams {
                bam_file: get_test_bam(),
                fasta_file: fasta,
                regions: Some(regions),
            };

            let readers = params.pileup_readers()?;
            let result = readers.segments(u64::MAX, 0).map(Iterator::count);
            assert!(result.is_err(), "explicit region absent from FASTA must be a hard error");
            Ok(())
        }
    }
}
