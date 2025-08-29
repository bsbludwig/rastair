#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

pub mod bed;
pub mod call;
pub mod call_reads;
pub mod convert;
pub mod io {
    pub mod formats;
    pub mod mpk;
    pub mod vcf_writer;
}
pub mod vcf;
pub mod utils {
    pub use rastair2_types::*;

    mod base_modification;
    pub mod file_helpers;
    pub use base_modification::MethylatedPositions;
    mod base_counter;
    pub use base_counter::Counter;

    pub mod logging;

    mod surrounding;
    pub use surrounding::surrounding_records;

    mod dedupe_reads;
    pub use dedupe_reads::ReadDeduplicator;
}
pub mod sequence;
