use clap_num::maybe_hex;
use rust_htslib::bam::Record;

#[derive(Debug, Clone, Default, clap::Args)]
pub struct ReadFlags {
    /// Include reads that match all of these bit-flags
    #[arg(short = 'f', long, value_parser=maybe_hex::<u16>)]
    pub include_flags: Option<u16>,

    /// Exclude reads that match any of these bit-flags
    #[arg(short = 'F', long, value_parser=maybe_hex::<u16>)]
    pub exclude_flags: Option<u16>,
}

impl ReadFlags {
    /// Check if the read matches the include flags and does not match the exclude flags
    pub fn filter(&self, record: &Record) -> bool {
        let flags = record.flags();
        let include = self.include_flags.map_or_else(|| true, |f| flags & f == f);
        let exclude = self.exclude_flags.map_or_else(|| false, |f| flags & f != 0);
        include && !exclude
    }
}
