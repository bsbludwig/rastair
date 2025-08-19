use crate::{
    bed::per_read::{BedReadsParams, PerRead},
    sequence::{ChunkRegion, ReaderParams, Readers, Region, Segment},
};
use color_eyre::{
    Result,
    eyre::{Context as _, ContextCompat, eyre},
};
use rastair2_types::{Strand, StrandFromRecord};
use rayon::iter::{ParallelBridge as _, ParallelIterator as _};
use rust_htslib::bam::{FetchDefinition, Read, Record, ext::BamRecordExtensions};
use smallvec::SmallVec;
use std::thread::{self, available_parallelism};
use tracing::{debug, instrument, warn};

mod flags;

#[derive(Debug, Clone, clap::Args)]
pub struct PerReadParams {
    // --- Input parameters ---
    #[command(flatten)]
    segments: ReaderParams,
    /// Maximum length of a segment in bases
    ///
    /// Used for splitting work between threads. Tweak this to adjust memory
    /// usage.
    #[arg(long, default_value_t = 100_000)]
    pub segment_max_length: u64,
    /// Number of bases to overlap between segments
    ///
    /// Helpful to avoid missing variants at the edges of segments.
    #[arg(long, default_value_t = 500)]
    pub segment_overlap: u64,

    // --- Calling parameters ---
    #[command(flatten)]
    read_flags: flags::ReadFlags,

    /// expected maximum read length. If set too short, some read positions
    /// might not get counted. Safest to set this a bit higher than the actual
    /// read length, to allow for indels in reads.
    #[arg(short='w', long, default_value_t = 200, value_parser = clap::value_parser!(u32).range(1..))]
    max_read_length: u32,

    /// Minimum mapping quality per aligned read
    #[arg(short = 'q', long, default_value_t = 1)]
    min_mapq: u8,

    /// Report reads with no CpGs in them
    #[arg(short = 'A', long)]
    all_reads: bool,

    /// Exclude reads where the orientation cannot be unambiguously determined
    #[arg(long)]
    exclude_ambiguous: bool,

    // --- Output parameters ---
    #[command(flatten)]
    bed_reads: BedReadsParams,

    // --- Other runtime parameters ---
    /// Number of threads to use for processing the BAM file. Will use all
    /// available threads when not specified.
    ///
    /// Note that VCF writing might use additional threads internally for compression.
    /// This can be overwritten with `--vcf-threads`.
    #[arg(short='@', long = "threads", default_value_t = available_parallelism().map(|n|n.get()).unwrap_or(1).max(1))]
    pub total_threads: usize,
}

#[instrument(level = "debug", skip_all)]
pub fn call_reads(params: &PerReadParams) -> Result<()> {
    let readers = params.segments.readers().wrap_err("Failed to read BAM/FASTA files")?;
    let regions: Vec<ChunkRegion> = readers
        .segments(params.segment_max_length, 0)
        .wrap_err("Could not fetch segments from BAM file")?
        .collect();

    // Process each region and write results to the BED file
    //
    // This is done in parallel to speed up the processing, so here are a few
    // comments on this works. There are two aspects to this: Collecting the
    // variant candidates and writing them to the output file.
    //
    // To use all CPU available, we use rayon to process the regions in
    // parallel. From there, we send ready-made BED records to here and write
    // them in order.
    let writer_threads = 1; // it's this one
    let worker_threads = params.total_threads.saturating_sub(writer_threads).max(1);

    // The connection between the processing threads and the writer is this
    // ordered channel. It buffers `Vec<PerReads>`s, alongside the index from
    // the parallel iterator.
    let (sender, receiver) = {
        // At least 10x buffer for records to account for reordering and processing time
        let buffer_size = worker_threads * 10;
        ordered_channel::bounded(buffer_size)
    };

    // Create a writer for the output
    let writer_thread = thread::Builder::new()
        .name("writer".to_string())
        .spawn({
            let params = params.clone();
            move || -> Result<()> {
                let mut bed_writer =
                    params.bed_reads.writer().wrap_err("Failed to open BED file")?;

                for records in receiver {
                    for row in records {
                        bed_writer.write_record(&row).wrap_err("Failed to write BED record")?;
                    }
                }

                bed_writer.close().wrap_err("Failed to close BED writer")?;
                Ok(())
            }
        })
        .wrap_err("Failed to spawn writer thread")?;

    // Run this in a custom rayon thread pool to control the number of threads
    // and be able to tweak parameters when profiling
    rayon::ThreadPoolBuilder::new()
        .thread_name(|idx| format!("worker-{idx}"))
        .num_threads(worker_threads)
        .build()
        .wrap_err("Failed to create thread pool for rayon")?
        .install(move || {
            regions.iter().enumerate().par_bridge().try_for_each_with(
                (sender, params),
                |(vcf_sender, params), (index, region)| {
                    process_region_wrapper(index, region, vcf_sender, params)
                },
            )
        })?;

    writer_thread
        .join()
        .map_err(|e| eyre!("{e:?}"))
        .wrap_err("Failed to join writer thread")?
        .wrap_err("writer thread error")?;

    Ok(())
}

/// Wrapper function for processing a region in a thread-safe manner.
#[instrument(level = "debug", skip_all, fields(region=%region.region))]
fn process_region_wrapper(
    index: usize,
    region: &ChunkRegion,
    sender: &mut ordered_channel::Sender<Vec<PerRead>>,
    params: &PerReadParams,
) -> Result<()> {
    thread_local! {
        /// Readers for the BAM and FASTA files, initialized per thread to avoid
        /// re-opening files or having a lock
        static READERS: std::cell::RefCell<Option<Readers>> = const { std::cell::RefCell::new(None) };
    }

    // Use thread-local readers to avoid re-opening files in each thread
    let res = READERS.with(|local_readers| -> Result<Vec<PerRead>> {
        let mut local_readers = local_readers.borrow_mut();
        let readers = {
            // Initialize thread-local readers first time the thread accesses them
            if local_readers.is_none() {
                let readers = params
                    .segments
                    .readers()
                    .wrap_err("Failed to open readers in worker thread")?;
                *local_readers = Some(readers);
            }
            local_readers.as_mut().wrap_err("Failed to access thread-local resources")?
        };

        process_region(readers, region, params)
    });

    let records = match res {
        Ok(records) => records,
        Err(e) => {
            warn!(error = format!("{e:#}"), "Failed to process region");
            // We still send an empty vector to the channel to increment the index
            Vec::new()
        }
    };

    if let Err(_err) = sender.send(index, records).wrap_err("Failed to send records to writer") {
        // the channel is closed, because the writer thread has finished
    }

    Ok(())
}

fn process_region(
    readers: &mut Readers,
    region: &ChunkRegion,
    params: &PerReadParams,
) -> Result<Vec<PerRead>> {
    let segment = readers
        .segment(region, params.segment_overlap)
        .wrap_err("Could not fetch segment from BAM file")?;

    FetchDefinition::try_from(&segment.region)
        .wrap_err("Could not convert region string")
        .and_then(|r| readers.bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))
        .wrap_err_with(|| format!("Could not fetch region `{}` from BAM file", region.region))?;

    let size = usize::try_from(segment.range.len())
        .wrap_err("Failed to convert segment range length to usize")?;
    let capacity_est = if params.all_reads {
        size
    } else {
        // Estimate capacity based on CpG count, assuming 1/10 reads have a CpG
        size / 10
    };
    let mut res = Vec::with_capacity(capacity_est);

    let mut record = Record::new();
    while let Some(result) = readers.bam.read(&mut record) {
        if let Err(e) = result {
            return Err(e).wrap_err("Failed to read BAM record");
        }
        if (record.pos() as u64) < segment.range.start {
            continue;
        }
        if !params.read_flags.filter(&record) {
            continue;
        }
        if record.mapq() < params.min_mapq {
            continue;
        }
        if record.seq_len() > params.max_read_length as usize {
            continue;
        }

        record.cache_cigar();
        let row = record_to_row(&record, &segment).wrap_err("Failed to read record")?;

        if params.all_reads || row.cpg_count > 0 {
            res.push(row);
        }
    }

    debug!(records = res.len(), "Processed region with {} records", res.len());

    Ok(res)
}

fn record_to_row(record: &Record, segment: &Segment) -> Result<PerRead> {
    let segment_start_pos =
        usize::try_from(segment.range.start).expect("segment range fits in usize");
    let ref_seq = &segment.sequence;
    let read_seq = record.seq();
    let cigar = record.cigar();

    let mut cpg_count = 0;
    let mut mod_cpgs = SmallVec::new();
    let mut unmod_cpgs = SmallVec::new();
    let mut snp_cpgs = SmallVec::new();

    for [pos_in_read, pos_in_ref] in record.aligned_pairs_full() {
        let Some(pos_in_read) = pos_in_read else {
            continue;
        };
        let Some(pos_in_ref) = pos_in_ref else {
            continue;
        };
        let pos_in_read = usize::try_from(pos_in_read).expect("position fits in usize");
        let pos_in_ref = usize::try_from(pos_in_ref).expect("position fits in usize");
        let idx = pos_in_ref
            .checked_sub(segment_start_pos)
            .wrap_err("Failed to calculate index for position")?;
        let read_base = read_seq[pos_in_read];
        let ref_base = ref_seq.get(idx).copied().unwrap_or(b'N');
        let orientation = StrandFromRecord::strand(record);

        if orientation == Strand::OT && ref_base == b'C' {
            let next_base = ref_seq.get(idx + 1).copied().unwrap_or(b'N');
            if next_base == b'G' {
                cpg_count += 1;
                match read_base {
                    b'C' => unmod_cpgs.push(pos_in_read),
                    b'T' => mod_cpgs.push(pos_in_read),
                    _ => snp_cpgs.push(pos_in_read),
                }
            }
        } else if orientation == Strand::OB && ref_base == b'G' {
            let prev_base =
                idx.checked_sub(1).and_then(|i| ref_seq.get(i)).copied().unwrap_or(b'N');
            if prev_base == b'C' {
                cpg_count += 1;
                match read_base {
                    b'G' => unmod_cpgs.push(pos_in_read),
                    b'A' => mod_cpgs.push(pos_in_read),
                    _ => snp_cpgs.push(pos_in_read),
                }
            }
        }
    }

    Ok(PerRead {
        region: Region {
            contig: segment.range.contig.clone(),
            start: u64::try_from(record.pos()).expect("pos fits in u64"),
            end: u64::try_from(cigar.end_pos()).expect("pos fits in u64"),
        },
        flag: record.flags(),
        mapq: record.mapq(),
        frag_length: record.insert_size().unsigned_abs(),
        read_length: record.seq_len(),
        read_id: String::from_utf8(Vec::from(record.qname())).unwrap_or_default(),
        cpg_count,
        mod_count: mod_cpgs.len(),
        mod_cpgs,
        unmod_cpgs,
        snp_cpgs,
    })
}
