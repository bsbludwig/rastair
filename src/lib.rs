#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

pub mod bed;
pub mod call;
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
}
pub mod sequence;
