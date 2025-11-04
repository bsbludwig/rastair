#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

pub mod bam;
pub mod bed;
pub mod call;
pub mod call_reads;
pub mod convert;
pub mod mbias;
pub mod io {
    pub mod formats;
    pub mod mpk;
    pub mod vcf_writer;
}
pub mod vcf;
pub mod utils {
    pub use rastair_types::*;

    mod base_counter;
    pub mod file_helpers;
    pub use base_counter::Counter;

    pub mod logging;

    mod surrounding;
    pub use surrounding::surrounding_records;

    mod dedupe_reads;
    pub use dedupe_reads::ReadDeduplicator;

    pub mod cli;

    pub mod conversion;
}
pub mod sequence;
