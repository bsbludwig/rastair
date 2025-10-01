use crate::{
    sequence::{ChunkRegion, ReaderParams, Region},
    utils::{cli, logging::ThisIsABug},
};
use clap::{Parser, value_parser};
use clio::ClioPath;
use color_eyre::eyre::{Context, Result};
use rastair_types::{Base, Strand, StrandFromRecord};
use rust_htslib::bam::{
    self, FetchDefinition, Header, Read, Record, Writer, ext::BamRecordExtensions as _,
    header::HeaderRecord,
};
use rustc_hash::FxHashMap;
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

#[tracing::instrument(level = "info", skip_all, fields(
    output = %params.output.path().display(),
))]
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
            .wrap_err_with(|| format!("Failed to rewrite region {}", segment.region))?;
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

    // Load all calls in this region into a map of position -> call
    // Since we're going region by region, all calls are from the same contig
    let calls: FxHashMap<u32, RastairCall> = calls_reader
        .query(&noodle_region)
        .wrap_err("failed to query calls file")?
        .iter()
        .map(|call| (call.pos, call.call.clone()))
        .collect();

    let mut record = Record::new();
    while let Some(result) = bam.read(&mut record) {
        if let Err(error) = result {
            warn!(%error, "Failed to read BAM record");
        }

        rewrite_record(&calls, &mut record).wrap_err("failed to rewrite record")?;

        writer.write(&record).wrap_err("failed to write record to new BAM file")?;
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
#[instrument(level = "debug", skip_all, fields(pos = record.pos()))]
fn rewrite_record(
    // All calls in the region of the record (they are all the same contig)
    calls: &FxHashMap<u32, RastairCall>,
    record: &mut Record,
) -> Result<()> {
    use Base::*;

    let strand = StrandFromRecord::strand(record);

    let mut seq = record.seq().as_bytes();
    let mut methylated_positions: SmallVec<u32, 10> = SmallVec::new();

    for [pos_in_read, pos_in_ref] in record.aligned_pairs_full() {
        let Some(pos_in_read) = pos_in_read else {
            continue;
        };
        let Some(pos_in_ref) = pos_in_ref else {
            continue;
        };
        let pos_in_ref = u32::try_from(pos_in_ref).expect("position fits in u32");
        let pos_in_read = u32::try_from(pos_in_read).expect("position fits in u32");

        if let Some(call) = calls.get(&pos_in_ref)
            && let RastairCall::Cpg { methylated, .. } = call
            && *methylated
        {
            let pos = pos_in_read as usize;
            let observed_base = Base::from(seq[pos]);

            // Let's rewrite the sequence to un-modify the bases and collect
            // methylated positions
            match strand {
                // If we are on the top strand, we expect to see T at methylated
                // C positions
                Strand::OT => {
                    if observed_base == T {
                        seq[pos] = *C;
                        methylated_positions.push(pos_in_read);
                    }
                }
                // If we are on the bottom strand, we expect to see A at
                // methylated G
                Strand::OB => {
                    if observed_base == A {
                        seq[pos] = *G;
                        methylated_positions.push(pos_in_read);
                    }
                }
                Strand::Unknown => continue,
            }
        }
    }

    let mods = MethylatedPositions::new(strand, &seq, &methylated_positions);

    record.set_seq(&seq);

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
