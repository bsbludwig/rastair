use std::path::Path;

use clap::{Parser, value_parser};
use clio::ClioPath;
use color_eyre::eyre::{Context, ContextCompat, Result};
use rastair::{
    bed::reader::{RastairBedReader, RastairCall, SimpleRastairBedRecord},
    utils::{Base, MethylatedPositions, logging::setup_logging},
};
use rastair_types::RegionString;
use rust_htslib::bam::{
    self, Header, Read, Record, Writer, ext::BamRecordExtensions as _, header::HeaderRecord,
};
use smallvec::SmallVec;

#[derive(Debug, Parser)]
struct Cli {
    /// Input file, SAM or BAM
    #[arg(value_parser=value_parser!(ClioPath).exists().is_file())]
    bam_file: ClioPath,
    /// Rastair's calls to determine methylation
    #[arg(value_parser=value_parser!(ClioPath).exists().is_file())]
    calls_file: ClioPath,
    /// Output file, BAM
    output_file: ClioPath,

    /// Region to fetch, e.g. "chr19:6105400-6105410"
    #[arg(short = 'l', long)]
    region: Option<RegionString>,

    /// Enable more logging
    ///
    /// You can also use the `RASTAIR_LOG` environment variable to configure
    /// logging in a more precise way. See the documentation of the
    /// `tracing-subscriber` library to learn more.
    #[arg(short, long, global = true)]
    verbose: bool,
}

fn main() -> Result<()> {
    let args = Cli::parse();
    setup_logging(args.verbose);

    let region = "chr19";
    let mut bam =
        bam::IndexedReader::from_path(args.bam_file.path()).wrap_err("failed to open bam file")?;
    bam.fetch(region).wrap_err("error fetching range")?;

    let mut calls_reader =
        RastairBedReader::new(args.calls_file.path()).wrap_err("failed to open calls file")?;
    let calls = calls_reader.query(&region.parse()?).wrap_err("failed to query calls")?;

    rewrite_bam(&mut bam, &args.output_file, &calls).wrap_err("rewrite")?;

    Ok(())
}

#[tracing::instrument(skip_all)]
fn rewrite_bam(
    bam: &mut bam::IndexedReader,
    output_file: &ClioPath,
    calls: &[SimpleRastairBedRecord],
) -> Result<()> {
    let header = {
        let mut header = Header::from_template(bam.header());
        header.push_record(
            HeaderRecord::new(b"PG")
                .push_tag(b"ID", "rastair.rewrite_bam")
                .push_tag(b"PN", "rastair")
                .push_tag(b"VN", env!("CARGO_PKG_VERSION"))
                .push_tag(b"CL", std::env::args().skip(1).collect::<Vec<_>>().join(" "))
                .push_tag(b"DS", "Rewrote BAM with methylation information"),
        );
        header
    };
    let mut writer = {
        if output_file.is_std() {
            Writer::from_stdout(&header, bam::Format::Bam)
        } else {
            Writer::from_path(output_file.path(), &header, bam::Format::Bam)
        }
    }
    .wrap_err("failed to create writer")?;
    writer
        .set_compression_level(bam::CompressionLevel::Fastest)
        .wrap_err("failed to set compression level")?;
    writer.set_threads(3).wrap_err("failed to set threads")?;

    let mut record = Record::new();
    while let Some(result) = bam.read(&mut record) {
        if let Err(e) = result {
            return Err(e).wrap_err("Failed to read BAM record");
        }

        // Find calls for this read
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

        if !methylated_positions.is_empty() {
            let mods = MethylatedPositions { base: Base::C, positions: methylated_positions };
            mods.apply_to_record(&mut record)?;
        }
        writer.write(&record).wrap_err("failed to write record")?;
    }
    drop(writer);

    Ok(())
}
