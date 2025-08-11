use color_eyre::Result;
use std::io::Write;

pub mod writer;

pub mod per_read;
pub mod rastair1;

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum BedFormat {
    /// BGZIP compressed file, usually `.bed.gz`
    BedGz,
    /// Regular BED file, usually `.bed`
    Bed,
}

pub trait BedRecord {
    const HEADER: &'static str;
    fn write<W: Write>(&self, writer: &mut W) -> Result<()>;
}
