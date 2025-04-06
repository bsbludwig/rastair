use bitflags::bitflags;
use clap::value_parser;
use clio::ClioPath;
use color_eyre::eyre::{Result, eyre};
use rastair2::{read, utils::RegionString};
use rust_htslib::bam::{self, Read as _, record::Cigar};
use smallvec::SmallVec;
use std::path::Path;
use tracing::{debug, info, instrument, trace, warn};
use tracing_subscriber::layer::SubscriberExt as _;

#[derive(Debug, clap::Parser)]
struct Cli {}

#[derive(Debug, clap::Subcommand)]
enum Subcommand {
    Call(CallParams),
}

#[derive(Debug, clap::Args)]
struct CallParams {
    /// A sorted and indexed bam file
    #[arg(value_name="BAM_FILE", value_parser=value_parser!(ClioPath).exists().is_file())]
    bam_file: ClioPath,

    /// A sorted and indexed (via samtools faidx) fasta file. Can be bgzip compressed, but requires both a gzi index and a fai index
    #[arg(short='r', long, value_name="FASTA_FILE", required=true, value_parser=value_parser!(ClioPath).exists().is_file())]
    fasta_file: ClioPath,

    /// Restrict to a specific chromosome or region of a chromosome. Format is "chr", "chr:start" or "chr:start-end", where start is 1-based and end is inclusive.
    #[arg(short = 'l', long)]
    region: Option<RegionString>,
}

fn main() -> Result<()> {
    color_eyre::install()?;
    let subscriber = tracing_subscriber::Registry::default()
        .with(tracing_error::ErrorLayer::default())
        .with(tracing_subscriber::fmt::Layer::default());

    tracing::subscriber::set_global_default(subscriber)?;

    read("test_data/test.bam".as_ref())?;

    Ok(())
}

// bitflags! {
//     #[derive(Debug, Clone, Copy)]
//     struct Flags: u16 {
//         const IS_PAIRED = 0x1;
//         const IS_PROPERLY_PAIRED = 0x2;
//         const IS_UNMAPPED = 0x4;
//         const MATE_IS_UNMAPPED = 0x8;
//         const IS_REVERSE_STRAND = 0x10;
//         const MATE_IS_REVERSE_STRAND = 0x20;
//         const IS_FIRST_IN_PAIR = 0x40;
//         const IS_SECOND_IN_PAIR = 0x80;
//         const IS_NOT_PRIMARY = 0x100;
//         const IS_FAILED = 0x200;
//         const IS_DUPLICATE = 0x400;
//         const IS_SUPPLEMENTAL = 0x800;
//     }
// }
