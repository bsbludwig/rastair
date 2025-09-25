use crate::{
    bed::reader::{RastairBedReader, RastairCall, SimpleRastairBedRecord},
    sequence::{ChunkRegion, ReaderParams, Region},
    utils::{cli, logging::ThisIsABug},
};
use clap::{Parser, value_parser};
use clio::ClioPath;
use color_eyre::eyre::{Context, Result};
use rust_htslib::bam::{
    self, FetchDefinition, Header, Read, Record, Writer, ext::BamRecordExtensions as _,
    header::HeaderRecord,
};
use smallvec::SmallVec;

mod base_modification;
pub use base_modification::MethylatedPositions;
use tracing::instrument;

#[derive(Debug, Parser)]
pub struct BamRewriteArgs {
    #[command(flatten)]
    segments: ReaderParams,
    /// Maximum length of a segment in bases
    ///
    /// Used for splitting work between threads. Tweak this to adjust memory
    /// usage.
    #[arg(long, default_value_t = 100_000)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub segment_max_length: u64,

    /// Rastair's calls to determine methylation
    #[arg(value_parser=value_parser!(ClioPath).exists().is_file())]
    #[arg(help_heading = cli::sections::INPUT)]
    calls_file: ClioPath,

    /// Output file
    #[arg(short = 'o', long, default_value = "-")]
    #[arg(help_heading = cli::sections::OUTPUT)]
    output: ClioPath,
}

#[tracing::instrument(level = "info", skip_all)]
pub fn rewrite(params: &BamRewriteArgs) -> Result<()> {
    let mut readers = params.segments.readers().wrap_err("Failed to read BAM/FASTA files")?;
    let regions: Vec<ChunkRegion> = readers
        .segments(params.segment_max_length, 0)
        .wrap_err("Could not fetch segments from BAM file")?
        .collect();

    let mut calls_reader =
        RastairBedReader::new(params.calls_file.path()).wrap_err("failed to open calls file")?;

    let output_file = &params.output;
    let mut writer = {
        let header = {
            let mut header = Header::from_template(readers.bam.header());
            add_rastair_header(&mut header);
            header
        };

        let mut writer = if output_file.is_std() {
            Writer::from_stdout(&header, bam::Format::Bam)
        } else {
            Writer::from_path(output_file.path(), &header, bam::Format::Bam)
        }
        .wrap_err("failed to create writer")?;
        writer
            .set_compression_level(bam::CompressionLevel::Fastest)
            .wrap_err("failed to set compression level")?;
        writer.set_threads(3).wrap_err("failed to set threads")?;
        writer
    };

    for segment in &regions {
        rewrite_region(&mut readers.bam, &mut calls_reader, &mut writer, &segment.region)
            .wrap_err("rewrite")?;
    }

    Ok(())
}

#[instrument(level = "debug", skip_all, fields(region = %region))]
fn rewrite_region(
    bam: &mut bam::IndexedReader,
    calls_reader: &mut RastairBedReader,
    writer: &mut rust_htslib::bam::Writer,
    region: &Region,
) -> Result<()> {
    FetchDefinition::try_from(region)
        .wrap_err("Could not convert region string")
        .and_then(|r| bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))
        .wrap_err_with(|| format!("Could not fetch region `{}` from BAM file", region))?;

    let noodle_region = region
        .clone()
        .try_into()
        .wrap_err("Failed to convert region representation")
        .this_is_a_bug()?;
    let calls = calls_reader.query(&noodle_region).wrap_err("failed to query calls")?;

    let mut record = Record::new();
    while let Some(result) = bam.read(&mut record) {
        if let e @ Err(_) = result {
            return e.wrap_err("Failed to read BAM record");
        }

        rewrite_record(bam, &calls, &mut record)?;

        writer.write(&record).wrap_err("failed to write record")?;
    }

    Ok(())
}

/// Rewrite a single BAM record with methylation information from calls
///
/// 1. Find all calls that are covered by our record
/// 2. Find the positions in the read that correspond to the called positions
/// 3. Build a `MethylatedPositions` struct
/// 4. Rewrite the sequence to un-modify the bases of methylated positions
///    (i.e., change T or A that is methylation evidence back to C or G)
/// 5. Apply the modifications to the record
#[instrument(level = "debug", skip_all, fields(tid = record.tid(), pos = record.pos()))]
fn rewrite_record(
    bam: &mut bam::IndexedReader,
    calls: &[SimpleRastairBedRecord],
    record: &mut Record,
) -> Result<(), color_eyre::eyre::Error> {
    let chr = bam.header().tid2name(record.tid() as u32);
    let called_positions: SmallVec<_, 10> = calls
        .iter()
        .filter(|call| chr == call.chrom.as_bytes())
        .filter(|call| {
            let read_start = record.pos();
            let read_end = read_start + record.insert_size();
            i64::from(call.pos) >= read_start && i64::from(call.pos) < read_end
        })
        .collect();
    let mut methylated_positions = SmallVec::new();
    for [pos_in_read, pos_in_ref] in record.aligned_pairs_full() {
        let Some(pos_in_read) = pos_in_read else {
            continue;
        };
        let Some(pos_in_ref) = pos_in_ref else {
            continue;
        };
        let pos_in_read = u32::try_from(pos_in_read).expect("position fits in u32");

        if let Some(called_pos) =
            called_positions.iter().find(|call| i64::from(call.pos) == pos_in_ref)
            && let RastairCall::Cpg { methylated, .. } = called_pos.call
            && methylated
        {
            methylated_positions.push(pos_in_read);
        }
    }
    let mods = MethylatedPositions::for_cpg_methylation(&*record, methylated_positions);
    mods.apply_to_record(record)?;
    Ok(())
}

/// Add a PG header line for rastair
fn add_rastair_header(header: &mut Header) {
    header.push_record(
        HeaderRecord::new(b"PG")
            .push_tag(b"ID", "rastair.rewrite_bam")
            .push_tag(b"PN", "rastair")
            .push_tag(b"VN", env!("CARGO_PKG_VERSION"))
            .push_tag(b"CL", std::env::args().skip(1).collect::<Vec<_>>().join(" "))
            .push_tag(b"DS", "Rewrote BAM with methylation information"),
    );
}
