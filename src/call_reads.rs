use crate::{
    bed::per_read::{BedReadsParams, PerRead},
    sequence::{ChunkRegion, ReaderParams, Region, Segment},
};
use color_eyre::{
    Result,
    eyre::{Context as _, ContextCompat},
};
use rastair2_types::{Strand, StrandFromRecord};
use rust_htslib::bam::{FetchDefinition, Read, Record, ext::BamRecordExtensions};
use smallvec::SmallVec;
use tracing::instrument;

mod flags;

#[derive(Debug, clap::Args)]
pub struct PerReadParams {
    // --- Input parameters ---
    #[command(flatten)]
    segments: ReaderParams,

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
}

#[instrument(level = "debug", skip_all)]
pub fn call_reads(params: &PerReadParams) -> Result<()> {
    let mut readers = params.segments.readers().wrap_err("Failed to read BAM/FASTA files")?;
    let regions: Vec<ChunkRegion> =
        readers.segments().wrap_err("Could not fetch segments from BAM file")?.collect();

    let mut bed_writer = params.bed_reads.writer().wrap_err("Failed to open BED file")?;

    for region in regions {
        let segment = readers.segment(&region).wrap_err("Could not fetch segment from BAM file")?;
        FetchDefinition::try_from(&segment)
            .wrap_err("Could not convert region string")
            .and_then(|r| readers.bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))
            .wrap_err_with(|| {
                format!("Could not fetch region `{}` from BAM file", region.region)
            })?;

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
                bed_writer.write_record(&row).wrap_err("Failed to write BED record")?;
            }
        }
    }

    Ok(())
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
